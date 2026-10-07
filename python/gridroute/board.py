"""Grid PCB layout toolkit: explicit placement, raster occupancy, a multi-layer router and its bookkeeping.

Layout scripts place footprints explicitly and route nets one at a time, in an order and on layers they choose
(the equivalent of interactive routing). This module provides the geometry, the obstacle-aware path search (on the
gridroute native kernels, with numpy fallbacks), rip-up, speculation and a persistent route cache; the result is
saved as JSON for a board writer (e.g. a KiCad pcbnew script).

Configure once per process with configure() (signal layers, plane nets, net classes, grid pitch, search options),
then build Board(w, h, comps, pin_net, footprints). See README.md ("Board toolkit").

Conventions: mm, KiCad board coordinates (x right, y down), rotations in degrees counter-clockwise as KiCad shows
them. Grid model: every signal layer is a raster of pitch G. `occ` holds copper grown by the class clearance
surplus (used for clearance checks), `core` holds the true copper (used for connectivity), `thru` marks vias and
plated holes. A route for a net may cross its own copper and must keep `half width + clearance + MARGIN` from
everything else; MARGIN absorbs the rasterisation error of off-grid pad edges.
"""
import heapq
import json
import math
import os
import threading
import weakref
from array import array
from collections import deque
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import gridroute as _native

KEEP = -1                     # occupancy value of keep-outs, holes, the board edge and clearance conflicts
TILE = 64                     # cells per side of a content-hash tile
SPEC_DEP = 2.0                # mm beyond a search window that can still change the search
SPEC_STATS = {'used': 0, 'stale': 0, 'submitted': 0}
_SPECS = weakref.WeakKeyDictionary()      # board -> {key: (ops index, phantom version, window, future)}
_TLS = threading.local()                  # .spec: running on a speculation thread (no search memo writes)
_SPEC_POOL = None


def _env(name, default):
    return os.environ.get('GRIDROUTE_' + name, default)


def configure(layers=('F.Cu', 'B.Cu'), planes=None, pitch=0.05, classes=None, net_class=None, layer_widths=None,
              edge_clear=0.3, hole_clear=0.25, pad_grow=0.025, margin=0.025, search=None, budget=None, field=None,
              accel=None, spec=None, cache=None, cache_name='', salt=b'', hole_to_hole=0.25, layer_clearances=None,
              copper_grow=0.):
    """Set the board rules and router options for this process (call before creating a Board).

    layers: signal layer names, top first. planes: {net: plane layer} for nets carried by planes (their pads get
    via drops, not routes). classes: {class: (width, clearance, via diameter, via drill)} in mm, with 'Default'.
    net_class: net name -> class name. layer_widths: {class: {layer: width}} for classes routed at a per-layer width
    (e.g. impedance-controlled pairs). edge_clear / hole_clear: copper to board edge / hole. pad_grow / margin: raster
    allowances. Router options (default from $GRIDROUTE_<NAME>): search 'hybrid' (default) or 'exact' (plain A*,
    the reference results; also sets REPRODUCE for layout scripts), budget (hybrid: expansions before the field
    takes over, 50000), field (field resolution factor, 2), accel 'on' / 'off' / 'verify' (native kernels; verify
    runs the numpy paths too and asserts equal results), spec (speculation threads, 4; 0 = off), cache (file path,
    or 'off'; default ~/Library/Caches/gridroute/<cache_name hash>-routes.pkl), salt (bytes mixed into cache keys,
    e.g. the caller's rule sources, so changing them starts a fresh cache)."""
    global LAYERS, NL, G, CLASSES, CLEAR, VIA_D, VIA_DRILL, PLANES, EDGE_CLEAR, HOLE_CLEAR, LAYER_WIDTHS, PAD_GROW
    global MARGIN, HOLE_TO_HOLE, LAYER_CLEARANCES, COPPER_GROW, _net_class, _GR, _VERIFY, _HYBRID, REPRODUCE, _BUDGET, _FIELD, _CACHE, _SALT, _SPEC_N, _SPEC_POOL, _CONFIG
    LAYERS = list(layers)
    NL = len(LAYERS)
    G = pitch
    CLASSES = dict(classes or {'Default': (0.15, 0.15, 0.5, 0.3)})
    CLEAR = min(r[1] for r in CLASSES.values())
    VIA_D, VIA_DRILL = CLASSES['Default'][2], CLASSES['Default'][3]
    PLANES = dict(planes or {})
    EDGE_CLEAR, HOLE_CLEAR = edge_clear, hole_clear
    HOLE_TO_HOLE = hole_to_hole
    LAYER_WIDTHS = {k: dict(v) for k, v in (layer_widths or {}).items()}
    LAYER_CLEARANCES = dict(layer_clearances or {})
    PAD_GROW, MARGIN = pad_grow, margin
    COPPER_GROW = copper_grow
    _net_class = net_class or (lambda n: 'Default')
    accel = accel or _env('ACCEL', 'on')
    search = search or _env('SEARCH', 'hybrid')
    _GR = _native if accel != 'off' and _native.available() else None
    _VERIFY = _GR is not None and accel == 'verify'
    _HYBRID = _GR is not None and search == 'hybrid'
    REPRODUCE = search == 'exact'
    _BUDGET = int(budget if budget is not None else _env('BUDGET', '50000'))
    _FIELD = int(field if field is not None else _env('FIELD', '2'))
    if _HYBRID:
        _GR.set_field_factor(_FIELD)
    _CACHE, _SALT = None, ''
    where = cache if cache is not None else _env('CACHE', '')
    if _GR is not None and where != 'off':
        import hashlib
        h = hashlib.blake2b(digest_size=16)
        h.update(open(os.path.abspath(__file__), 'rb').read())
        h.update(salt)
        h.update(repr((LAYERS, PLANES, G, CLASSES, LAYER_WIDTHS, EDGE_CLEAR, HOLE_CLEAR, PAD_GROW, MARGIN, _HYBRID,
                       HOLE_TO_HOLE, LAYER_CLEARANCES, COPPER_GROW, _BUDGET, _FIELD, _GR.__file__, os.path.getmtime(_GR._LIB._name))).encode())
        path = where or os.path.join(os.path.expanduser('~/Library/Caches/gridroute'),
                                     hashlib.sha1(cache_name.encode()).hexdigest()[:12] + '-routes.pkl')
        _CACHE, _SALT = _GR.SearchCache(path), h.hexdigest()
    if _SPEC_POOL is not None:
        _SPEC_POOL.shutdown(wait=True)
    _SPEC_N = int(spec if spec is not None else _env('SPEC', '4')) if _GR is not None else 0
    _SPEC_POOL = ThreadPoolExecutor(_SPEC_N) if _SPEC_N > 0 else None
    _CONFIG = dict(layers=list(LAYERS), planes=dict(PLANES), pitch=G, classes=dict(CLASSES),
                   net_class=_net_class, layer_widths=LAYER_WIDTHS, edge_clear=EDGE_CLEAR, hole_clear=HOLE_CLEAR,
                   pad_grow=PAD_GROW, margin=MARGIN, search=search, budget=_BUDGET, field=_FIELD, accel=accel,
                   spec=_SPEC_N, cache=where, cache_name=cache_name, salt=salt, hole_to_hole=HOLE_TO_HOLE,
                   layer_clearances=LAYER_CLEARANCES,copper_grow=COPPER_GROW)


def net_class(net):
    return _net_class(net)


def short(net):
    return net.rstrip('/').split('/')[-1]


def rot(px, py, deg):
    """Rotate a footprint-local offset by the footprint orientation (KiCad: CCW on screen, y down)."""
    a = math.radians(deg)
    c, s = math.cos(a), math.sin(a)
    return px * c + py * s, -px * s + py * c


class Part:
    def __init__(self, ref, fpdef, net_of):
        self.ref, self.fp = ref, fpdef
        self.net_of = net_of
        self.x = self.y = self.rot = 0.0
        self.placed = False

    def pads(self):
        """World pads: dicts with x, y, w, h (axis-aligned after rotation), shape, layers, net, num. Memoised per
        position and rotation (callers must not modify the dicts)."""
        key = (self.x, self.y, self.rot)
        if getattr(self, '_pads_key', None) != key:
            self._pads, self._pads_key = self._make_pads(), key
        return list(self._pads)

    def _make_pads(self):
        out = []
        for p in self.fp['pads']:
            dx, dy = rot(p['x'], p['y'], self.rot)
            r = (p['rot'] + p.get('angle', 0.) + self.rot) % 180
            orthogonal = min(abs(r), abs(r - 90), abs(r - 180)) < 1e-6
            swap = orthogonal and abs(r - 90) < 1e-6
            w, h = (p['h'], p['w']) if swap else (p['w'], p['h'])
            dw, dh = p.get('drill_w', p['drill']), p.get('drill_h', p['drill'])
            if swap:
                dw, dh = dh, dw
            out.append(dict(num=p['num'], x=self.x + dx, y=self.y + dy, w=w, h=h, shape=p['shape'],
                            layers=p['layers'], npth=p['npth'], drill=p['drill'], ref=self.ref, rr=p.get('rr', 0.0),
                            drill_w=dw, drill_h=dh, angle=0. if orthogonal else r,
                            original_ref=p.get('original_ref',self.ref), original_pin=p.get('original_pin',p['num']),
                            escape_width=p.get('escape_width'), mask_layers=p.get('mask_layers',()),
                            mask_margin=p.get('mask_margin',0.),mask_margins=p.get('mask_margins',{}),
                            clearance=p.get('clearance',0.),
                            polygons=[dict(points=[tuple(a+b for a,b in zip(rot(x,y,self.rot),(self.x+dx,self.y+dy))) for x,y in poly['points']],
                                           holes=[[tuple(a+b for a,b in zip(rot(x,y,self.rot),(self.x+dx,self.y+dy))) for x,y in hole] for hole in poly.get('holes',[])])
                                      for poly in p.get('polygons',[])],
                            net=self.net_of.get((self.ref, p['num'])) if p['num'] else None))
        return out

    def pad(self, num):
        ps = [p for p in self.pads() if p['num'] == num]
        return ps[0] if ps else None

    def courtyard(self):
        key = (self.x, self.y, self.rot)
        if getattr(self, '_cy_key', None) != key:
            self._cy, self._cy_key = self._make_courtyard(), key
        return self._cy

    def _make_courtyard(self):
        x0, y0, x1, y1 = self.fp['courtyard']
        pts = [rot(x, y, self.rot) for x, y in ((x0, y0), (x1, y0), (x1, y1), (x0, y1))]
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        return self.x + min(xs), self.y + min(ys), self.x + max(xs), self.y + max(ys)


def pad_layers(p):
    return list(range(NL)) if p['layers'] == 'all' else [0 if p['layers'] == 'F' else NL - 1]


def offset_polyline(pts, d):
    """Polyline offset by d to the left of its direction as seen on screen (y down), with mitred corners."""
    def nrm(a, b):
        dx, dy = b[0] - a[0], b[1] - a[1]
        l = math.hypot(dx, dy)
        return (dy / l, -dx / l)
    out = []
    n0 = nrm(pts[0], pts[1])
    out.append((pts[0][0] + n0[0] * d, pts[0][1] + n0[1] * d))
    for p0, p1, p2 in zip(pts, pts[1:], pts[2:]):
        na, nb = nrm(p0, p1), nrm(p1, p2)
        bx, by = na[0] + nb[0], na[1] + nb[1]
        c = 1 + na[0] * nb[0] + na[1] * nb[1]
        out.append((p1[0] + d * bx / c, p1[1] + d * by / c))
    n1 = nrm(pts[-2], pts[-1])
    out.append((pts[-1][0] + n1[0] * d, pts[-1][1] + n1[1] * d))
    return out


def _seg_dist(x, y, a, b):
    vx, vy = b[0] - a[0], b[1] - a[1]
    l2 = vx * vx + vy * vy
    t = max(0.0, min(1.0, ((x - a[0]) * vx + (y - a[1]) * vy) / l2)) if l2 > 0 else 0.0
    return math.hypot(x - a[0] - t * vx, y - a[1] - t * vy)


def _pad_distance(x, y, p):
    """Distance from a point to the actual pad copper, not its bounding box."""
    if p['shape']=='polygon':
        def inside(ring):
            value=False
            for (ax,ay),(bx,by) in zip(ring,ring[1:]+ring[:1]):
                if (ay>y)!=(by>y) and x<(bx-ax)*(y-ay)/(by-ay)+ax:value=not value
            return value
        distance=math.inf
        for poly in p['polygons']:
            if inside(poly['points']) and not any(inside(h) for h in poly.get('holes',[])):return 0.
            for ring in [poly['points']]+poly.get('holes',[]):
                distance=min(distance,min(_seg_dist(x,y,a,b) for a,b in zip(ring,ring[1:]+ring[:1])))
        return distance
    if p['shape'] in ('circle', 'oval'):
        rr = min(p['w'], p['h']) / 2
    else:
        rr = min(p.get('rr', 0.), p['w'] / 2, p['h'] / 2)
    dx, dy = rot(x - p['x'], y - p['y'], -p.get('angle', 0.))
    dx = max(abs(dx) - p['w'] / 2 + rr, 0.)
    dy = max(abs(dy) - p['h'] / 2 + rr, 0.)
    return max(0., math.hypot(dx, dy) - rr)


def dedup(pts, eps=1e-6):
    out = [pts[0]]
    for p in pts[1:]:
        if math.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) > eps:
            out.append(p)
    return out


def _collinear_merged(pts, eps=1e-9):
    """Drop interior points that lie on the straight line between their neighbours."""
    out = [pts[0]]
    for q, r in zip(pts[1:], pts[2:] + [None]):
        if r is not None:
            a = out[-1]
            cross = (q[0] - a[0]) * (r[1] - a[1]) - (q[1] - a[1]) * (r[0] - a[0])
            dot = (q[0] - a[0]) * (r[0] - q[0]) + (q[1] - a[1]) * (r[1] - q[1])
            if abs(cross) < eps and dot >= 0:
                continue
        out.append(q)
    return out


class RouteError(Exception):
    pass


class Board:
    def __init__(self, w, h, comps, pin_net, footprints):
        """w, h: board size (mm). comps: {ref: {'footprint': name, ...}}; pin_net: {(ref, pin): net} (see
        gridroute.kicad.load_netlist). footprints: {name: {'pads': [...], 'courtyard': [x0, y0, x1, y1]}}, pads as
        dicts with num, x, y, w, h, rot (multiples of 90 with the part), shape ('rect', 'roundrect', 'circle',
        'oval'), rr (corner radius), layers ('all' through-hole, 'F' top or 'B' bottom), npth, drill, in footprint
        coordinates at 0 degrees (gridroute.kicad.dump_footprints writes this from KiCad libraries)."""
        self.w, self.h = w, h
        self.config = dict(_CONFIG)
        self.pitch = G
        self.nx, self.ny = int(math.ceil(w / G)) + 1, int(math.ceil(h / G)) + 1
        self.occ = np.zeros((NL, self.ny, self.nx), dtype=np.int16)
        self.core = np.zeros((NL, self.ny, self.nx), dtype=np.int16)
        self.thru = np.zeros((self.ny, self.nx), dtype=np.int16)
        self.smd = np.zeros((self.ny, self.nx), dtype=np.int16)      # SMD pad copper (either side): no vias
        self.no_via = np.zeros((self.ny, self.nx), dtype=np.int16)   # explicit via-only rule areas
        self._has_via_areas = False
        self.phantom = np.zeros((self.ny, self.nx), dtype=bool)      # temporary F.Cu / via obstacles (templates)
        self.comps, self.pin_net = comps, pin_net
        self.fpdefs = footprints
        self.nets = sorted(set(self.pin_net.values()))
        self.net_id = {n: i + 1 for i, n in enumerate(self.nets)}
        self.cls = {n: net_class(n) for n in self.nets}
        self._id_clearance = {self.net_id[n]:CLASSES[c][1] for n,c in self.cls.items()}
        self._wide_ids = tuple(i for i, c in self._id_clearance.items() if c > CLEAR)   # nets above the minimum
        self.parts = {r: Part(r, self.fpdefs[c['footprint']], self.pin_net) for r, c in self.comps.items()}
        self.tracks, self.vias = [], []
        self.keepouts, self.rule_keepouts, self.texts, self.rects, self.zones = [], [], [], [], []
        self.labels = []              # functional silkscreen labels placed next to a part by pcb_write
        self.failed = []
        self.overlaps = []
        self.turn = {}                # ref -> extra rotation applied by place() (orientation fixes)
        self.next_dummy = -2
        self._th = np.zeros(((self.ny + TILE - 1) // TILE, (self.nx + TILE - 1) // TILE, 2), dtype=np.uint64)
        self._dirty = np.ones(self._th.shape[:2], dtype=bool)
        # every rasterising operation in order (kind, args, painted bbox, track / via dict): rip() removes nets'
        # copper by clearing a region and replaying the remaining operations there
        self.ops, self._op_win = [], None
        self._edges()

    def _edges(self):
        e = max(1,int(round((EDGE_CLEAR - CLEAR) / G)) + 1)
        for a in (self.occ,):
            a[:, :e, :] = KEEP
            a[:, -e:, :] = KEEP
            a[:, :, :e] = KEEP
            a[:, :, -e:] = KEEP

    def _op(self, kind, args, item=None):
        """Run rasterising operation `kind` and log it with the union of the windows it painted."""
        self._op_win = [self.ny, self.nx, -1, -1]
        getattr(self, '_do_' + kind)(*args)
        w, self._op_win = self._op_win, None
        self.ops.append((kind, args, tuple(w), item))

    # ------------------------------------------------------------------ helpers
    def width(self, net, layer=None):
        lw = LAYER_WIDTHS.get(self.cls[net])
        if layer and lw:
            return lw.get(layer, CLASSES[self.cls[net]][0])
        return CLASSES[self.cls[net]][0]

    def clearance(self, net):
        return CLASSES[self.cls[net]][1]

    def extra(self, net):
        """Clearance surplus of the net's class over the default (grown into `occ` around its copper)."""
        return self.clearance(net) - CLEAR if net else 0.0

    def via_size(self, net):
        return CLASSES[self.cls[net]][2:4]

    def _nid(self, net):
        if net is None:
            self.next_dummy -= 1
            return self.next_dummy
        return self.net_id[net]

    def _win(self, x0, y0, x1, y1):
        i0, j0 = max(int(math.floor(y0 / G)), 0), max(int(math.floor(x0 / G)), 0)
        i1, j1 = min(int(math.ceil(y1 / G)), self.ny - 1), min(int(math.ceil(x1 / G)), self.nx - 1)
        return i0, j0, i1, j1

    def _grid(self, win):
        i0, j0, i1, j1 = win
        return (np.arange(i0, i1 + 1) * G)[:, None], (np.arange(j0, j1 + 1) * G)[None, :]

    def _touch(self, win):
        """Record a painted window: content-hash tiles go dirty, the running operation's bbox grows."""
        i0, j0, i1, j1 = win
        self._dirty[i0 // TILE:i1 // TILE + 1, j0 // TILE:j1 // TILE + 1] = True
        if self._op_win is not None:
            w = self._op_win
            w[0], w[1], w[2], w[3] = min(w[0], i0), min(w[1], j0), max(w[2], i1), max(w[3], j1)

    def _paint_seg(self, arr, L, a, b, r, nid):
        """Paint the capsule a-b of radius r (mm) with nid: _mask_seg + _paint, natively when available."""
        if _GR is None:
            win, m = self._mask_seg(a[0], a[1], b[0], b[1], r)
            self._paint(arr, L, win, m, nid)
            return
        self._touch(_GR.paint_capsule(arr, L or 0, a, b, r, G, nid, keep=KEEP if arr is self.occ else None))

    def _paint(self, arr, L, win, m, nid):
        i0, j0, i1, j1 = win
        self._touch(win)
        w = arr[L, i0:i1 + 1, j0:j1 + 1] if arr.ndim == 3 else arr[i0:i1 + 1, j0:j1 + 1]
        clash = m & (w != 0) & (w != nid)
        w[m & (w == 0)] = nid
        if arr is self.occ:
            w[clash] = KEEP

    def _mask_rect(self, x0, y0, x1, y1, grow=0.0):
        win = self._win(x0 - grow, y0 - grow, x1 + grow, y1 + grow)
        ys, xs = self._grid(win)
        m = (xs >= x0 - grow - 1e-9) & (xs <= x1 + grow + 1e-9) & (ys >= y0 - grow - 1e-9) & (ys <= y1 + grow + 1e-9)
        return win, m

    def _mask_seg(self, x0, y0, x1, y1, r):
        win = self._win(min(x0, x1) - r, min(y0, y1) - r, max(x0, x1) + r, max(y0, y1) + r)
        ys, xs = self._grid(win)
        vx, vy = x1 - x0, y1 - y0
        l2 = vx * vx + vy * vy
        t = np.clip(((xs - x0) * vx + (ys - y0) * vy) / l2, 0, 1) if l2 > 0 else 0
        dx, dy = xs - (x0 + t * vx), ys - (y0 + t * vy)
        return win, dx * dx + dy * dy <= r * r + 1e-9

    def _pad_mask(self, p, grow):
        if p['shape']=='polygon':return self._polygon_mask(p['polygons'],grow)
        if p.get('angle', 0.):
            a = math.radians(p['angle'])
            c, s = math.cos(a), math.sin(a)
            w, h = p['w'] / 2, p['h'] / 2
            bx, by = abs(c) * w + abs(s) * h + grow, abs(s) * w + abs(c) * h + grow
            win = self._win(p['x'] - bx, p['y'] - by, p['x'] + bx, p['y'] + by)
            ys, xs = self._grid(win)
            x, y = xs - p['x'], ys - p['y']
            u, v = x * c - y * s, x * s + y * c
            rr = min(w, h) if p['shape'] in ('circle', 'oval') else min(p.get('rr', 0.), w, h)
            dx, dy = np.maximum(np.abs(u) - w + rr, 0), np.maximum(np.abs(v) - h + rr, 0)
            return win, dx * dx + dy * dy <= (rr + grow) ** 2 + 1e-9
        if p['shape'] in ('circle', 'oval'):
            rr = min(p['w'], p['h']) / 2
            ax, ay = (p['w'] / 2 - rr, 0) if p['w'] >= p['h'] else (0, p['h'] / 2 - rr)
            return self._mask_seg(p['x'] - ax, p['y'] - ay, p['x'] + ax, p['y'] + ay, rr + grow)
        rr = min(p.get('rr', 0.0), p['w'] / 2, p['h'] / 2)
        if rr > 0:
            # rounded rectangle: points within rr + grow of the rectangle shrunk by rr
            win = self._win(p['x'] - p['w'] / 2 - grow, p['y'] - p['h'] / 2 - grow,
                            p['x'] + p['w'] / 2 + grow, p['y'] + p['h'] / 2 + grow)
            ys, xs = self._grid(win)
            dx = np.maximum(np.abs(xs - p['x']) - (p['w'] / 2 - rr), 0)
            dy = np.maximum(np.abs(ys - p['y']) - (p['h'] / 2 - rr), 0)
            return win, dx * dx + dy * dy <= (rr + grow) ** 2 + 1e-9
        return self._mask_rect(p['x'] - p['w'] / 2, p['y'] - p['h'] / 2, p['x'] + p['w'] / 2, p['y'] + p['h'] / 2, grow)

    def _polygon_mask(self, polygons, grow=0.):
        points=[q for poly in polygons for q in poly['points']]
        win=self._win(min(x for x,y in points)-max(grow,0),min(y for x,y in points)-max(grow,0),
                      max(x for x,y in points)+max(grow,0),max(y for x,y in points)+max(grow,0))
        ys,xs=self._grid(win);shape=(win[2]-win[0]+1,win[3]-win[1]+1)
        result=np.zeros(shape,dtype=bool)
        for poly in polygons:
            inside=np.zeros(shape,dtype=bool);distance=np.full(shape,np.inf)
            for ri,ring in enumerate([poly['points']]+list(poly.get('holes',[]))):
                ring_inside=np.zeros(shape,dtype=bool)
                for (ax,ay),(bx,by) in zip(ring,ring[1:]+ring[:1]):
                    if by!=ay:ring_inside^=((ay>ys)!=(by>ys)) & (xs<(bx-ax)*(ys-ay)/(by-ay)+ax)
                    vx,vy=bx-ax,by-ay;den=vx*vx+vy*vy
                    t=np.clip(((xs-ax)*vx+(ys-ay)*vy)/den,0,1) if den else 0
                    distance=np.minimum(distance,(xs-ax-t*vx)**2+(ys-ay-t*vy)**2)
                if ri==0:inside=ring_inside
                else:inside &= ~ring_inside
            result |= (inside|(distance<=grow*grow+1e-12)) if grow>=0 else (inside&(distance>=grow*grow))
        return win,result

    # ------------------------------------------------------------------ placement
    def place(self, ref, x, y, rot=0.0, gap=0.15):
        """Place a footprint; courtyards closer than `gap` to an already placed part are recorded in self.overlaps."""
        p = self.parts[ref]
        assert not p.placed, ref
        rot += self.turn.get(ref, 0)
        p.x, p.y, p.rot, p.placed = x, y, rot % 360, True
        self._pad_index = None
        self._hole_pads = None
        a = p.courtyard()
        for r, q in self.parts.items():
            if q.placed and r != ref:
                b = q.courtyard()
                if a[0] < b[2] + gap and b[0] < a[2] + gap and a[1] < b[3] + gap and b[1] < a[3] + gap:
                    self.overlaps.append((ref, r))
        for pad in p.pads():
            self._op('pad', (pad, None if pad['npth'] else self._nid(pad['net'])))
        return p

    def _do_pad(self, pad, nid):
        required=max(HOLE_CLEAR,pad.get('clearance',0.))
        if pad['npth']:
            hole = dict(pad, shape='oval', w=pad.get('drill_w', pad['drill']), h=pad.get('drill_h', pad['drill']))
            win, m = self._pad_mask(hole, max(0.,required - CLEAR) + PAD_GROW)
            for L in range(NL):
                self._paint(self.occ, L, win, m, KEEP)
            return
        Ls = pad_layers(pad)
        surplus=max(self.extra(pad['net']),pad.get('clearance',0.)-CLEAR,0.)
        win, m = self._pad_mask(pad, PAD_GROW + surplus)
        cwin, cm = self._pad_mask(pad, 0.0)
        for L in Ls:
            self._paint(self.occ, L, win, m, nid)
            self._paint(self.core, L, cwin, cm, nid)
        if pad['layers'] == 'all':
            self._paint(self.thru, None, cwin, cm, nid)
            hole = dict(pad, shape='oval', w=pad.get('drill_w', pad['drill']), h=pad.get('drill_h', pad['drill']))
            hwin, hm = self._pad_mask(hole, max(0., required - CLEAR) + PAD_GROW)
            for L in Ls:
                self._paint(self.occ, L, hwin, hm, nid)
        else:
            self._paint(self.smd, None, cwin, cm, 1)
        # Some footprints expose an aperture on the opposite side to their
        # copper (castellated/edge connectors, or legacy footprint mistakes).
        # A foreign trace through that aperture becomes exposed copper and a
        # solder-mask bridge. Preserve it as an obstacle, never connectivity.
        for layer in pad.get('mask_layers',()):
            L = LAYERS.index(layer)
            margin=pad.get('mask_margins',{}).get(layer,pad.get('mask_margin',0.))
            if L not in Ls or margin>CLEAR+surplus:
                mwin, mm = self._pad_mask(pad, max(0.,margin-CLEAR) + PAD_GROW)
                self._paint(self.occ,L,mwin,mm,nid)

    # ------------------------------------------------------------------ copper
    def add_track(self, net, layer, pts, width=None):
        width = width or self.width(net, layer)
        # CAD units are nanometres. Four-decimal rounding can put an imperial
        # class width below its actual minimum (e.g. .37592 -> .3759 mm).
        width=math.ceil(width*1e6-1e-6)/1e6
        pts=[[round(x,6),round(y,6)] for x,y in pts]
        t = dict(net=net, layer=layer, width=width, pts=pts)
        self._op('track', (net, LAYERS.index(layer), [tuple(q) for q in pts], width), t)
        self.tracks.append(t)

    def _do_track(self, net, L, pts, width):
        nid = self.net_id[net]
        for a, b in zip(pts, pts[1:]):
            self._paint_seg(self.occ, L, a, b, width / 2 + self.extra(net) + COPPER_GROW, nid)
            self._paint_seg(self.core, L, a, b, width / 2, nid)

    def add_via(self, net, x, y, d=None, drill=None):
        d0, dr0 = self.via_size(net)
        d, drill = d or d0, drill or dr0
        v = dict(net=net, x=round(x, 6), y=round(y, 6), d=d, drill=drill)
        self._op('via', (net, v['x'], v['y'], d, drill), v)
        self.vias.append(v)

    def _do_via(self, net, x, y, d, drill):
        nid = self.net_id[net]
        for L in range(NL):
            reach = max(d / 2 + self.extra(net), drill / 2 + HOLE_CLEAR - CLEAR)
            self._paint_seg(self.occ, L, (x, y), (x, y), reach + COPPER_GROW, nid)
            self._paint_seg(self.core, L, (x, y), (x, y), d / 2, nid)
        self._paint_seg(self.thru, None, (x, y), (x, y), d / 2, nid)

    def keepout(self, x0, y0, x1, y1, layers=None, rule=False):
        """Router keep-out (all or the given signal layers); rule=True also writes a KiCad rule area."""
        self._op('keepout', (x0, y0, x1, y1, tuple(layers or LAYERS)))
        self.keepouts.append([x0, y0, x1, y1])
        if rule:
            self.rule_keepouts.append([x0, y0, x1, y1])

    def _do_keepout(self, x0, y0, x1, y1, layers):
        win, m = self._mask_rect(x0, y0, x1, y1)
        for l in layers:
            self._paint(self.occ, LAYERS.index(l), win, m, KEEP)

    def keepout_polygon(self, points, layers=None, holes=(), tracks=True, vias=True, clearance=None):
        """Preserve a CAD copper obstacle or rule area, including polygon holes.

        Track exclusions enter the clearance raster. Via-only areas are kept
        separately, so they do not unnecessarily block legal surface traces.
        """
        self._op('polygon', (tuple(map(tuple, points)), tuple(layers or LAYERS),
                             tuple(tuple(map(tuple, h)) for h in holes), tracks, vias, clearance))

    def _do_polygon(self, points, layers, holes, tracks, vias, clearance=None):
        polygons=[dict(points=points,holes=holes)]
        win,mask=self._polygon_mask(polygons,PAD_GROW+max(0.,(clearance or 0.)-CLEAR))
        if tracks:
            cwin,cm=self._polygon_mask(polygons,0.)
            for layer in layers:
                self._paint(self.occ, LAYERS.index(layer), win, mask, KEEP)
                self._paint(self.core, LAYERS.index(layer), cwin, cm, KEEP)
        if vias:
            self._has_via_areas = True
            self._paint(self.no_via, None, win, mask, 1)

    def reserve(self, net, layer, a, b, r):
        """Reserve the capsule a-b (radius r, mm) for `net` on `layer` (a name, a list of names, or None for every
        layer): other nets keep their clearance from it, as from copper, but it carries no connectivity (e.g. an
        escape corridor or a planned via site of a fine-pitch pin). Not copper: replace_copper() keeps it;
        release_reservations() removes it."""
        Ls = list(range(NL)) if layer is None else [LAYERS.index(l) for l in ([layer] if isinstance(layer, str) else layer)]
        self._op('reserve', (net, tuple(Ls), tuple(a), tuple(b), r), None)

    def _do_reserve(self, net, Ls, a, b, r):
        nid = self.net_id[net]
        win, m = self._mask_seg(a[0], a[1], b[0], b[1], r)
        i0, j0, i1, j1 = win
        self._touch(win)
        if self.__dict__.get('rsv') is None:
            self.rsv = np.zeros_like(self.core)
        for L in ([Ls] if isinstance(Ls, int) else Ls):
            w = self.occ[L, i0:i1 + 1, j0:j1 + 1]
            w[m & (w == 0)] = nid      # never marks a clash: a reservation must not block anything already there
            # like true copper (core), but without connectivity: _blocked enforces surplus class clearances on it
            v = self.rsv[L, i0:i1 + 1, j0:j1 + 1]
            v[m & (v == 0)] = nid

    def reservations(self, nets=None):
        """Logged reservation operations (kind, args, ...) of `nets` (all when None)."""
        return [op for op in self.ops if op[0] == 'reserve' and (nets is None or op[1][0] in nets)]

    def release_reservations(self, nets=None):
        """Remove the reservations of `nets` (all when None); returns how many were removed."""
        return self._rip_ops([k for k, op in enumerate(self.ops) if op[0] == 'reserve' and (nets is None or op[1][0] in nets)])

    def add_neck(self, points, width, nets, layers=None, inset=None):
        """Neck-down region: tracks of `nets` may use `width` (narrower than their class) on `layers` inside the
        polygon `points`, e.g. a CAD rule letting power tracks touching a fine-pitch package's courtyard be thin.
        The raster region is inset (default one pitch) so that every thin segment really touches the polygon."""
        inset = G if inset is None else inset
        win, mask = self._polygon_mask([dict(points=[tuple(q) for q in points])], -inset)
        Ls = [LAYERS.index(l) for l in (layers or LAYERS) if l in LAYERS]
        if mask.any() and Ls:
            self.__dict__.setdefault('necks', []).append(dict(win=win, mask=mask, layers=Ls, width=width,
                                                              nets=frozenset(nets)))

    def _neck_widths(self, net, win):
        """[NL, H, W] neck width per cell of window win for `net` (0 = class width)."""
        out = None
        for nk in self.__dict__.get('necks', ()):
            if net not in nk['nets']:
                continue
            a, b = nk['win'], win
            i0, j0, i1, j1 = max(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), min(a[3], b[3])
            if i0 > i1 or j0 > j1:
                continue
            if out is None:
                out = np.zeros((NL, b[2] - b[0] + 1, b[3] - b[1] + 1), dtype=np.float32)
            reg = nk['mask'][i0 - a[0]:i1 - a[0] + 1, j0 - a[1]:j1 - a[1] + 1]
            for L in nk['layers']:
                view = out[L, i0 - b[0]:i1 - b[0] + 1, j0 - b[1]:j1 - b[1] + 1]
                view[reg & ((view == 0) | (view > nk['width']))] = nk['width']
        return out

    def keepout_ring(self, x, y, radius, width, layers=None, clearance=None):
        """An unfilled circular CAD stroke, with its interior still routable."""
        if radius <= 0 or width < 0:
            raise ValueError('Ring radius must be positive and width nonnegative')
        self._op('ring', (x, y, radius, width, tuple(layers or LAYERS), clearance))

    def _do_ring(self, x, y, radius, width, layers, clearance):
        grow=PAD_GROW+max(0.,(clearance or 0.)-CLEAR)
        outer=radius+width/2+grow
        win=self._win(x-outer,y-outer,x+outer,y+outer)
        ys,xs=self._grid(win)
        distance=np.abs(np.hypot(xs-x,ys-y)-radius)
        mask=distance<=width/2+grow
        core=distance<=width/2
        for layer in layers:
            self._paint(self.occ,LAYERS.index(layer),win,mask,KEEP)
            self._paint(self.core,LAYERS.index(layer),win,core,KEEP)
        self._has_via_areas=True
        self._paint(self.no_via,None,win,mask,1)

    def text(self, s, x, y, layer='F.SilkS', size=1.0):
        self.texts.append(dict(text=s, x=x, y=y, layer=layer, size=size))

    def label(self, ref, text, size=1.0, inside=False):
        self.labels.append(dict(ref=ref, text=text, size=size, inside=inside))

    def rect(self, x0, y0, x1, y1, layer='F.Fab'):
        self.rects.append(dict(x0=x0, y0=y0, x1=x1, y1=y1, layer=layer))

    def zone(self, net, layer, pts, priority=1):
        """Extra copper zone (e.g. a power pour on an inner layer); the router treats its area as that net."""
        self.zones.append(dict(net=net, layer=layer, pts=pts, priority=priority))

    # ------------------------------------------------------------------ search masks
    def _blocked(self, L, nid, win, r):
        """Cells of window win where copper of net nid (or any of a tuple of nets) with reach r (centre to edge +
        clearance) would clash."""
        if _GR is not None:
            out = self._blocked_gr([L], nid, win, r)[0]
            if _VERIFY:
                assert np.array_equal(out, self._blocked_np(L, nid, win, r)), ('_blocked', L, nid, win, r)
            return out
        return self._blocked_np(L, nid, win, r)

    def _blocked_gr(self, Ls, nid, win, r, reduce_or=False):
        """_blocked for several layers in one native call ([len(Ls), H, W], or the OR over them [H, W])."""
        ph = self.phantom.any() and 0 in Ls
        excluded=nid if isinstance(nid, tuple) else (nid,)
        own=max([CLEAR]+[self._id_clearance.get(n,CLEAR) for n in excluded]+
                [LAYER_CLEARANCES.get(LAYERS[L],0.) for L in Ls])
        # Class clearances combine with max(), never by summing both surpluses.
        # The grown raster enforces the obstacle's class; the searching net's
        # class is enforced on the grown raster of minimum-clearance copper
        # (pads keep PAD_GROW, tracks and vias get none, as before classes
        # combined) and on the true copper of wider-class nets, whose grown
        # raster already carries their own surplus. Single-class boards keep one call.
        base_r=r-max(0.,own-CLEAR)
        out = _GR.dilate(self.occ, Ls, win, _GR.disc_span(base_r / G), excl=excluded,
                          extra=self.phantom if ph else None, extra_on=[1] + [0] * (NL - 1) if ph else None,
                          reduce_or=reduce_or)
        if own>CLEAR:
            out |= _GR.dilate(self.occ,Ls,win,_GR.disc_span(r/G),excl=excluded+self._wide_ids,
                              extra=self.phantom if ph else None,extra_on=[1]+[0]*(NL-1) if ph else None,
                              reduce_or=reduce_or)
            out |= _GR.dilate(self.core,Ls,win,_GR.disc_span(r/G),excl=excluded,reduce_or=reduce_or)
            rsv=self.__dict__.get('rsv')
            if rsv is not None and rsv[:,win[0]:win[2]+1,win[1]:win[3]+1].any():
                out |= _GR.dilate(rsv,Ls,win,_GR.disc_span(r/G),excl=excluded,reduce_or=reduce_or)
        return out

    def _blocked_np(self, L, nid, win, r):
        excluded=nid if isinstance(nid,tuple) else (nid,)
        own=max([CLEAR,LAYER_CLEARANCES.get(LAYERS[L],0.)]+[self._id_clearance.get(n,CLEAR) for n in excluded])
        out=self._dilate_np(self.occ,L,nid,win,r-max(0.,own-CLEAR))
        if own>CLEAR:
            out|=self._dilate_np(self.occ,L,excluded+self._wide_ids,win,r)
            out|=self._dilate_np(self.core,L,nid,win,r,phantom=False)
            if self.__dict__.get('rsv') is not None:out|=self._dilate_np(self.rsv,L,nid,win,r)
        return out

    def _dilate_np(self, source, L, nid, win, r, phantom=True):
        i0, j0, i1, j1 = win
        rc = int(math.ceil(r / G))
        pi0, pj0 = max(i0 - rc, 0), max(j0 - rc, 0)
        pi1, pj1 = min(i1 + rc, self.ny - 1), min(j1 + rc, self.nx - 1)
        src = source[L, pi0:pi1 + 1, pj0:pj1 + 1]
        other = src != 0
        for n in (nid if isinstance(nid, tuple) else (nid,)):
            other &= src != n
        if phantom and L == 0 and self.phantom.any():
            other |= self.phantom[pi0:pi1 + 1, pj0:pj1 + 1]
        H, W = other.shape
        cs = np.zeros((H, W + 1), dtype=np.int32)
        np.cumsum(other, axis=1, out=cs[:, 1:])
        dil = np.zeros_like(other)
        rr = (r / G) ** 2
        jj = np.arange(W)
        for di in range(-rc, rc + 1):
            if di * di >= rr - 1e-9:
                continue
            w = int(math.floor(math.sqrt(rr - di * di) - 1e-9))
            lo, hi = np.clip(jj - w, 0, W), np.clip(jj + w + 1, 0, W)
            hd = (cs[:, hi] - cs[:, lo]) > 0
            if di >= 0:
                dil[:H - di] |= hd[di:]
            else:
                dil[-di:] |= hd[:H + di]
        return dil[i0 - pi0:i1 - pi0 + 1, j0 - pj0:j1 - pj0 + 1]

    def _near_smd(self, win, r):
        """Cells of window win within r of SMD pad copper (via centres there would land on a pad)."""
        if _GR is not None:
            out = _GR.dilate(self.smd, [0], win, _GR.disc_span(r / G, strict=False))[0]
            if _VERIFY:
                assert np.array_equal(out, self._near_smd_np(win, r)), ('_near_smd', win, r)
            return out
        return self._near_smd_np(win, r)

    def _near_foreign_copper(self, nid, win, radius):
        """Hole-to-copper distance, independent of either copper netclass.

        Use the true copper raster: _blocked applies class-clearance surplus
        compensation and therefore must not receive a drill-clearance radius.
        """
        radius+=PAD_GROW
        if _GR is not None:
            return _GR.dilate(self.core,list(range(NL)),win,_GR.disc_span(radius/G),
                              excl=(nid,),reduce_or=True)
        out=np.zeros((win[2]-win[0]+1,win[3]-win[1]+1),bool)
        for L in range(NL):out|=self._dilate_np(self.core,L,nid,win,radius)
        return out

    def _near_holes(self, win, radius, exclude_net=None):
        """Exact capsule distance to drilled holes, including same-net PTHs.

        radius is the new object's radius plus its required hole clearance.
        Existing holes are not copper: same-net copper exemptions must never
        permit overlapping drills. Slots use their two physical drill sizes.
        Only intersecting local holes allocate a window-sized distance mask.
        """
        ys, xs = self._grid(win)
        out = np.zeros((win[2] - win[0] + 1, win[3] - win[1] + 1), dtype=bool)
        for x, y, w, h, angle, net in self._holes():
            if exclude_net is not None and net == exclude_net:
                continue
            r = min(w, h) / 2 + radius
            ax, ay = max(0., (w - h) / 2), max(0., (h - w) / 2)
            if abs(w-h) < 1e-9:
                angle = 0.
            c, s = math.cos(math.radians(angle)), math.sin(math.radians(angle))
            bx, by = abs(c) * ax + abs(s) * ay, abs(s) * ax + abs(c) * ay
            if x + bx + r < win[1] * G or x - bx - r > win[3] * G or y + by + r < win[0] * G or y - by - r > win[2] * G:
                continue
            # A drill only affects a small rectangle even when the search covers
            # the entire board. Restrict the distance arithmetic to that box.
            i0=max(win[0],int(math.floor((y-by-r)/G)));i1=min(win[2],int(math.ceil((y+by+r)/G)))
            j0=max(win[1],int(math.floor((x-bx-r)/G)));j1=min(win[3],int(math.ceil((x+bx+r)/G)))
            yy=ys[i0-win[0]:i1-win[0]+1];xx=xs[:,j0-win[1]:j1-win[1]+1]
            # Preserve row/column broadcasting for round and orthogonal holes.
            # Expanding both coordinates to HxW for every ordinary via makes
            # this otherwise small local check needlessly memory-bound.
            if angle:
                u, v = (xx - x) * c - (yy - y) * s, (xx - x) * s + (yy - y) * c
            else:
                u, v = xx - x, yy - y
            dx, dy = np.maximum(np.abs(u) - ax, 0), np.maximum(np.abs(v) - ay, 0)
            out[i0-win[0]:i1-win[0]+1,j0-win[1]:j1-win[1]+1] |= dx * dx + dy * dy < r * r
        return out

    def _holes(self):
        if getattr(self, '_hole_pads', None) is None:
            pads = (p for part in self.parts.values() if part.placed for p in part.pads())
            self._hole_pads = [(p['x'],p['y'],p.get('drill_w',p['drill']),p.get('drill_h',p['drill']),
                                p.get('angle',0.),p['net']) for p in pads if p['drill'] > 0]
        return self._hole_pads + [(v['x'],v['y'],v['drill'],v['drill'],0.,v['net']) for v in self.vias]

    def _via_areas(self, win, r):
        if not self._has_via_areas:
            return np.zeros((win[2]-win[0]+1,win[3]-win[1]+1),dtype=bool)
        if _GR is not None:
            return _GR.dilate(self.no_via, [0], win, _GR.disc_span(r/G, strict=False))[0]
        return self._near_smd_np(win, r, source=self.no_via)

    def _near_smd_np(self, win, r, source=None):
        i0, j0, i1, j1 = win
        rc = int(math.ceil(r / G))
        pi0, pj0 = max(i0 - rc, 0), max(j0 - rc, 0)
        pi1, pj1 = min(i1 + rc, self.ny - 1), min(j1 + rc, self.nx - 1)
        src = (self.smd if source is None else source)[pi0:pi1 + 1, pj0:pj1 + 1] != 0
        H, W = src.shape
        cs = np.zeros((H, W + 1), dtype=np.int32)
        np.cumsum(src, axis=1, out=cs[:, 1:])
        dil = np.zeros_like(src)
        rr = (r / G) ** 2
        jj = np.arange(W)
        for di in range(-rc, rc + 1):
            if di * di > rr:
                continue
            w = int(math.floor(math.sqrt(rr - di * di)))
            lo, hi = np.clip(jj - w, 0, W), np.clip(jj + w + 1, 0, W)
            hd = (cs[:, hi] - cs[:, lo]) > 0
            if di >= 0:
                dil[:H - di] |= hd[di:]
            else:
                dil[-di:] |= hd[:H + di]
        return dil[i0 - pi0:i1 - pi0 + 1, j0 - pj0:j1 - pj0 + 1]

    def _reach(self, net, half, layer=None):
        return half + max(self.clearance(net),LAYER_CLEARANCES.get(layer,0.)) + MARGIN

    def _island(self, nid, seeds, win):
        """Mask (NL, H, W) of window cells of net nid connected to the seed cells through copper and holes."""
        if _GR is not None:
            out = _GR.island(self.core, self.thru, nid, win, seeds)
            if _VERIFY:
                assert np.array_equal(out, self._island_np(nid, seeds, win)), ('_island', nid, win)
            return out
        return self._island_np(nid, seeds, win)

    def _island_np(self, nid, seeds, win):
        i0, j0, i1, j1 = win
        H, W = i1 - i0 + 1, j1 - j0 + 1
        own = (self.core[:, i0:i1 + 1, j0:j1 + 1] == nid)
        thr = (self.thru[i0:i1 + 1, j0:j1 + 1] == nid)
        seen = np.zeros_like(own)
        ownb = bytearray(own.astype(np.uint8).ravel().tobytes())
        thrb = thr.astype(np.uint8).ravel().tobytes()
        seenb = bytearray(len(ownb))
        HW = H * W
        q = deque()
        for (L, i, j) in seeds:
            s = L * HW + i * W + j
            if 0 <= i < H and 0 <= j < W and ownb[s] and not seenb[s]:
                seenb[s] = 1
                q.append(s)
        while q:
            s = q.popleft()
            L, c = divmod(s, HW)
            i, j = divmod(c, W)
            nb = []
            if j > 0:
                nb.append(s - 1)
            if j < W - 1:
                nb.append(s + 1)
            if i > 0:
                nb.append(s - W)
            if i < H - 1:
                nb.append(s + W)
            if thrb[c]:
                nb.extend(L2 * HW + c for L2 in range(NL) if L2 != L)
            for t in nb:
                if ownb[t] and not seenb[t]:
                    seenb[t] = 1
                    q.append(t)
        seen = np.frombuffer(bytes(seenb), dtype=np.uint8).reshape(NL, H, W).astype(bool)
        return seen

    def _pad_seeds(self, p, win):
        i0, j0 = win[0], win[1]
        # Supported pads are filled, axis-aligned convex shapes. One interior
        # cell per copper layer reaches the same island as seeding every pixel.
        # Keep the exhaustive fallback for clipped pads and malformed/overlapping
        # input: these cases must not silently change connectivity semantics.
        i, j = int(round(p['y'] / G)), int(round(p['x'] / G))
        Ls = pad_layers(p)
        nid = self.net_id.get(p.get('net'))
        if (p['shape'] in ('rect', 'roundrect', 'circle', 'oval')
                and min(p['w'], p['h']) >= G
                and win[0] <= i <= win[2] and win[1] <= j <= win[3]
                and 0 <= i < self.ny and 0 <= j < self.nx
                and nid is not None and all(self.core[L, i, j] == nid for L in Ls)):
            return [(L, i - i0, j - j0) for L in Ls]
        cwin, cm = self._pad_mask(p, 0.0)
        ii, jj = np.nonzero(cm)
        return [(L, i + cwin[0] - i0, j + cwin[1] - j0) for L in pad_layers(p) for i, j in zip(ii, jj)]

    # ------------------------------------------------------------------ routing
    def pads_of(self, net):
        if getattr(self, '_pad_index', None) is None:
            idx = {}
            for r, p in self.parts.items():
                if p.placed:
                    for pd in p.pads():
                        idx.setdefault(pd['net'], []).append(pd)
            self._pad_index = idx
        return list(self._pad_index.get(net, ()))

    def pad_components(self, net, pads=None):
        """Canonical component labels at pads, in pad order (0, 1, ...).

        Native flood labelling avoids allocating a full-board mask per net.
        Off-grid/tiny pads without an occupied centre use the exhaustive island
        fallback. An unrepresented pad receives its own disconnected label.
        """
        pads = self.pads_of(net) if pads is None else pads
        cells = [(pad_layers(p)[0], int(round(p['y'] / G)), int(round(p['x'] / G))) for p in pads]
        nid = self.net_id[net]
        if _GR is not None and all(0 <= i < self.ny and 0 <= j < self.nx and self.core[L, i, j] == nid
                                   for L, i, j in cells):
            raw = _GR.label_cells(self.core, self.thru, nid, cells)
            canonical, out = {}, []
            for k, label in enumerate(raw):
                key = int(label) if label >= 0 else ('missing', k)
                out.append(canonical.setdefault(key, len(canonical)))
            return out
        win = (0, 0, self.ny - 1, self.nx - 1)
        seeds = [self._pad_seeds(p, win) for p in pads]
        labels = [-1] * len(pads)
        component = 0
        for k, ss in enumerate(seeds):
            if labels[k] >= 0:
                continue
            labels[k] = component
            island = self._island(nid, ss, win)
            for j in range(k + 1, len(pads)):
                if labels[j] < 0 and any(island[s] for s in seeds[j]):
                    labels[j] = component
            component += 1
        return labels

    def route(self, net, a, b, **kw):
        """Connect the copper island of pad a=(ref, num) to the island of pad b."""
        pa, pb = self.parts[a[0]].pad(a[1]), self.parts[b[0]].pad(b[1])
        assert pa['net'] == net and pb['net'] == net, (net, a, b, pa['net'], pb['net'])
        return self._route(net, [pa], [pb], **kw)

    def connect(self, net, layers=None, only=None, strict=False, **kw):
        """Route every placed pad of the net into one tree (nearest pad first). only: restrict to these refs.
        A connection that fails on `layers` is retried with a wider window on every layer, or on `layers` again
        when strict=True (e.g. a board whose bottom layer is a plane that tracks must not use)."""
        pads = [p for p in self.pads_of(net) if not only or p['ref'] in only]
        if len(pads) < 2:
            return True
        tree = [pads[0]]
        rest = pads[1:]
        ok = True
        # fast mode: label the net's existing copper components at the pads once (whole board); a pad already in
        # the tree's component needs no search. Others go through _route, which still checks connection itself.
        comp, joined = {}, set()
        if _GR is not None and not REPRODUCE:
            lab = _GR.label_cells(self.core, self.thru, self.net_id[net],
                                  [(pad_layers(q)[0], int(round(q['y'] / G)), int(round(q['x'] / G))) for q in pads])
            comp = {id(q): int(c) for q, c in zip(pads, lab) if c >= 0}
            joined = {comp[id(pads[0])]} if id(pads[0]) in comp else set()
        while rest:
            rest.sort(key=lambda p: min(math.hypot(p['x'] - t['x'], p['y'] - t['y']) for t in tree))
            p = rest.pop(0)
            if comp.get(id(p), -1) in joined:
                tree.append(p)
                continue
            if not self._route(net, [p], tree, layers=layers, **kw):
                wide = kw.get('margin', 3.0) * 3
                if not self._route(net, [p], tree, layers=layers if strict else None, margin=wide,
                                   **{k: v for k, v in kw.items() if k != 'margin'}):
                    self.failed.append((net, p['ref'], p['num']))
                    ok = False
                    continue
            tree.append(p)
            if id(p) in comp:
                joined.add(comp[id(p)])     # routed: its whole component now belongs to the tree
        return ok

    AXIAL_MAX_PAD = 2.5       # pads larger than this (mm) may be entered anywhere, like thermal or power pads

    def _axial_entries(self, nid, pads, win, blk):
        """Pad entry rule (Board.pad_entry = 'axial'): a route may leave or enter an SMD pad only on its two
        centre lines and runs straight for a short halo beyond the pad edge. Blocks the other cells of each pad
        and its halo in blk (this search only; the net's other copper there stays usable) and returns the
        pads for snapping the emitted track onto the exact centre line."""
        e = max(G, self.__dict__.get('pad_entry_halo', 0.1))
        snaps, seen = [], set()
        for p in pads:
            key = (p['ref'], p['num'])
            if key in seen or p['layers'] == 'all' or p['shape'] == 'polygon' or p.get('angle', 0.):
                continue
            seen.add(key)
            if max(p['w'], p['h']) > self.AXIAL_MAX_PAD or min(p['w'], p['h']) < G:
                continue
            L = pad_layers(p)[0]
            gw, gm = self._pad_mask(p, e)
            pw, pm = self._pad_mask(p, 0.)
            i0, j0 = max(gw[0], win[0]), max(gw[1], win[1])
            i1, j1 = min(gw[2], win[2]), min(gw[3], win[3])
            if i0 > i1 or j0 > j1:
                continue
            region = gm[i0 - gw[0]:i1 - gw[0] + 1, j0 - gw[1]:j1 - gw[1] + 1].copy()
            ci, cj = int(round(p['y'] / G)), int(round(p['x'] / G))
            if i0 <= ci <= i1:
                region[ci - i0, :] = False
            if j0 <= cj <= j1:
                region[:, cj - j0] = False
            pad = np.zeros_like(region)
            a0, b0 = max(pw[0], i0), max(pw[1], j0)
            a1, b1 = min(pw[2], i1), min(pw[3], j1)
            if a0 <= a1 and b0 <= b1:
                pad[a0 - i0:a1 - i0 + 1, b0 - j0:b1 - j0 + 1] = pm[a0 - pw[0]:a1 - pw[0] + 1, b0 - pw[1]:b1 - pw[1] + 1]
            region &= ~(~pad & (self.core[L, i0:i1 + 1, j0:j1 + 1] == nid))     # own stubs, vias, tracks stay
            blk[L, i0 - win[0]:i1 - win[0] + 1, j0 - win[1]:j1 - win[1] + 1] |= region
            snaps.append((p['x'], p['y'], ci, cj, max(p['w'], p['h']) / 2 + G))
        return snaps

    def _content_key(self, win, *args):
        """Content hash of the board arrays over the tiles covering win, plus args (route cache key)."""
        ti0, tj0, ti1, tj1 = win[0] // TILE, win[1] // TILE, win[2] // TILE, win[3] // TILE
        d = self._dirty[ti0:ti1 + 1, tj0:tj1 + 1]
        if d.any():
            tl = np.argwhere(d) + (ti0, tj0)
            self._th[tl[:, 0], tl[:, 1]] = _GR.tile_hashes([self.occ, self.core, self.thru, self.smd, self.no_via, self.phantom], TILE, tl)
            d[:] = False
        return _GR.hash_arrays([self._th[ti0:ti1 + 1, tj0:tj1 + 1]]) + repr((_SALT, win, args, self._holes()))

    def _route(self, net, src_pads, dst_pads, layers=None, width=None, margin=3.0, via_cost=1.0, layer_cost=None,
               pref=None, turn=0.1, weight=1.2, window=None, max_exp=4_000_000, extra_src=None, extra_dst=None):
        """extra_dst: (x, y) points on existing copper of the net (e.g. a bus track) that belong to the target."""
        args = (layers, width, margin, via_cost, layer_cost, pref, turn, weight, window, max_exp, extra_src, extra_dst)
        key = None
        if _CACHE is not None:
            pts = [(p['x'], p['y']) for p in src_pads + dst_pads] + list(extra_dst or [])
            x0, y0, x1, y1 = window or (min(p[0] for p in pts) - margin, min(p[1] for p in pts) - margin,
                                        max(p[0] for p in pts) + margin, max(p[1] for p in pts) + margin)
            iw = self._win(x0 - 25, y0 - 25, x1 + 25, y1 + 25)
            key = self._content_key(iw, net, [(p['ref'], p['num'], p['x'], p['y']) for p in src_pads],
                                    [(p['ref'], p['num'], p['x'], p['y']) for p in dst_pads], layers, width, margin,
                                    via_cost, sorted((layer_cost or {}).items()), sorted((pref or {}).items()), turn,
                                    weight, window, max_exp, extra_src, extra_dst)
            hit = _CACHE.get(key)
            if hit is not None:
                return self._apply(hit)
        res = self._spec_take(net, src_pads, dst_pads, args)
        axial = self.__dict__.get('pad_entry') == 'axial' and not self.__dict__.get('_relaxed_entry')
        exp0 = _native.stats['astar_exp']
        if res is None:
            res = self._route_search(net, src_pads, dst_pads, *args)
        # relax only a search boxed in at its pads (no legal entry, or a tiny search), not ordinary congestion
        if res[0] == 'fail' and axial and (len(res) > 1 or _native.stats['astar_exp'] - exp0 < 2000):
            # the exceptional case (e.g. a corner pad): any entry angle and point, counted for reporting
            self._relaxed_entry = True
            try:
                res = self._route_search(net, src_pads, dst_pads, *args)
            finally:
                self._relaxed_entry = False
            if res[0] != 'fail':
                self.relaxed_entries = getattr(self, 'relaxed_entries', 0) + 1
        if key is not None:
            _CACHE.put(key, res)
        return self._apply(res)

    def _apply(self, res):
        """Commit a search result: ('path', emit arguments...), ('conn',) (already connected) or ('fail',)."""
        if res[0] == 'path':
            self._emit(*res[1:])
        return res[0] != 'fail'

    @staticmethod
    def _route_args(layers=None, width=None, margin=3.0, via_cost=1.0, layer_cost=None, pref=None, turn=0.1,
                    weight=1.2, window=None, max_exp=4_000_000, extra_src=None, extra_dst=None):
        return (layers, width, margin, via_cost, layer_cost, pref, turn, weight, window, max_exp, extra_src, extra_dst)

    def _spec_key(self, net, src_pads, dst_pads, args):
        return repr((net, [(p['ref'], p['num']) for p in src_pads], [(p['ref'], p['num']) for p in dst_pads], args))

    def speculate(self, net, src_pads, dst_pads, **kw):
        """Start _route(net, src_pads, dst_pads, **kw)'s search on a worker thread against the board as it is now."""
        if _SPEC_POOL is None or (_CACHE is not None and _CACHE.hits > 3 * _CACHE.misses):
            return               # a warm route cache answers faster than a speculation can start
        args = self._route_args(**kw)
        specs = _SPECS.setdefault(self, {})
        key = self._spec_key(net, src_pads, dst_pads, args)
        if key in specs:
            return
        layers, width, margin, window, extra_dst = args[0], args[1], args[2], args[8], args[11]
        pts = [(p['x'], p['y']) for p in src_pads + dst_pads] + list(extra_dst or [])
        x0, y0, x1, y1 = window or (min(p[0] for p in pts) - margin, min(p[1] for p in pts) - margin,
                                    max(p[0] for p in pts) + margin, max(p[1] for p in pts) + margin)
        dep = self._win(x0 - SPEC_DEP, y0 - SPEC_DEP, x1 + SPEC_DEP, y1 + SPEC_DEP)
        fut = _SPEC_POOL.submit(self._spec_run, net, src_pads, dst_pads, args)
        specs[key] = (len(self.ops), getattr(self, '_phantom_ver', 0), dep, fut)
        SPEC_STATS['submitted'] += 1

    def speculate_connect(self, net, **kw):
        """speculate() the first connection connect(net, **kw) will search (its nearest pad to the first pad)."""
        pads = self.pads_of(net)
        if len(pads) < 2:
            return
        rest = sorted(pads[1:], key=lambda p: math.hypot(p['x'] - pads[0]['x'], p['y'] - pads[0]['y']))
        self.speculate(net, [rest[0]], [pads[0]], layers=kw.pop('layers', None), **kw)

    def _spec_run(self, net, src_pads, dst_pads, args):
        _TLS.spec = True
        try:
            return self._route_search(net, src_pads, dst_pads, *args)
        except Exception:          # a speculation never breaks the run: the main thread recomputes
            return None
        finally:
            _TLS.spec = False

    def _spec_take(self, net, src_pads, dst_pads, args):
        """A speculative result for this search if its window is untouched since it started, else None."""
        specs = _SPECS.get(self)
        if not specs:
            return None
        sp = specs.pop(self._spec_key(net, src_pads, dst_pads, args), None)
        if sp is None:
            return None
        k0, pv, dep, fut = sp
        res = fut.result()
        if res is None or pv != getattr(self, '_phantom_ver', 0) or any(
                w[0] <= dep[2] and dep[0] <= w[2] and w[1] <= dep[3] and dep[1] <= w[3] for _, _, w, _ in self.ops[k0:]):
            SPEC_STATS['stale'] += 1
            return None
        SPEC_STATS['used'] += 1
        return res

    def spec_clear(self):
        """Wait for and drop pending speculations (before rip-up or copying the board)."""
        for k0, pv, dep, fut in _SPECS.pop(self, {}).values():
            fut.result()

    def _route_search(self, net, src_pads, dst_pads, layers, width, margin, via_cost, layer_cost, pref, turn, weight,
                      window, max_exp, extra_src, extra_dst):
        explicit = width is not None
        width = width or self.width(net)
        half = width / 2
        nid = self.net_id[net]
        pts = [(p['x'], p['y']) for p in src_pads + dst_pads] + list(extra_dst or [])
        if window:
            x0, y0, x1, y1 = window
        else:
            x0 = min(p[0] for p in pts) - margin
            y0 = min(p[1] for p in pts) - margin
            x1 = max(p[0] for p in pts) + margin
            y1 = max(p[1] for p in pts) + margin
        win = self._win(x0, y0, x1, y1)
        i0, j0, i1, j1 = win
        H, W = i1 - i0 + 1, j1 - j0 + 1
        # islands are traced over a larger window so that copper leaving the search window and coming back (a
        # pair on an inner layer, a trunk) still counts as connected; then cropped to the search window
        iw = self._win(x0 - 25, y0 - 25, x1 + 25, y1 + 25)
        di, dj = i0 - iw[0], j0 - iw[1]
        seeds_s = [s for p in src_pads for s in self._pad_seeds(p, iw)] + \
            [(L, i + di, j + dj) for L, i, j in (extra_src or [])]
        seeds_d = [s for p in dst_pads for s in self._pad_seeds(p, iw)]
        for x, y in extra_dst or []:
            i, j = int(round(y / G)), int(round(x / G))
            if iw[0] <= i <= iw[2] and iw[1] <= j <= iw[3]:
                seeds_d += [(L, i - iw[0], j - iw[1]) for L in range(NL) if self.core[L, i, j] == nid]
        if _GR is not None and not _VERIFY:      # traced over iw, returned for the search window only
            isl_s = _GR.island_crop(self.core, self.thru, nid, iw, seeds_s, win)
            isl_d = _GR.island_crop(self.core, self.thru, nid, iw, seeds_d, win)
        else:
            isl_s = self._island(nid, seeds_s, iw)[:, di:di + H, dj:dj + W]
            isl_d = self._island(nid, seeds_d, iw)[:, di:di + H, dj:dj + W]
        if (isl_s & isl_d).any():
            return ('conn',)
        Ls = [LAYERS.index(l) for l in (layers or LAYERS)]
        r = self._reach(net, half)
        blk = np.ones((NL, H, W), dtype=bool)
        vd, vdr = self.via_size(net)
        vr = self._reach(net, vd / 2)
        per_layer = bool(LAYER_CLEARANCES or (not explicit and self.cls[net] in LAYER_WIDTHS))
        if per_layer:
            for L in Ls:
                wh = half if explicit else self.width(net,LAYERS[L])/2
                blk[L] = self._blocked(L,nid,win,self._reach(net,wh,LAYERS[L]))
            vok = np.ones((H,W),dtype=bool)
            for L in range(NL):
                vok &= ~self._blocked(L,nid,win,self._reach(net,vd/2,LAYERS[L]))
        elif _GR is not None:
            blk[Ls] = self._blocked_gr(Ls, nid, win, r)
            vok = ~self._blocked_gr(list(range(NL)), nid, win, vr, reduce_or=True)
        if not per_layer and (_GR is None or _VERIFY):
            if _VERIFY:
                blk_gr, vok_gr = blk.copy(), vok.copy()
            for L in Ls:
                blk[L] = self._blocked(L, nid, win, r)
            vok = np.ones((H, W), dtype=bool)
            for L in range(NL):
                vok &= ~self._blocked(L, nid, win, vr)
            if _VERIFY:
                assert np.array_equal(blk, blk_gr) and np.array_equal(vok, vok_gr), ('_route masks', net, win)
        snap = None
        if self.__dict__.get('pad_entry') == 'axial' and not self.__dict__.get('_relaxed_entry'):
            snap = self._axial_entries(nid, src_pads + dst_pads, win, blk)
        neck = None if explicit else self._neck_widths(net, win)
        if neck is not None:
            for L in Ls:
                for wn in np.unique(neck[L][neck[L] > 0]):
                    region = neck[L] == wn
                    thin = self._blocked(L, nid, win, self._reach(net, float(wn) / 2, LAYERS[L]))
                    blk[L] &= ~(region & ~thin)
        vok &= ~self._near_smd(win, vd / 2 + 0.1)     # no via-in-pad, not even on the net's own pads
        vok &= ~self._near_holes(win, vdr / 2 + HOLE_TO_HOLE + 0.01)
        vok &= ~self._via_areas(win, vd / 2 + MARGIN)
        # Copper and drill clearances are independent; a small annulus does not
        # make the hole's clearance to foreign copper disappear.
        hr = vdr / 2 + HOLE_CLEAR + MARGIN
        if hr > vr:
            vok &= ~self._near_foreign_copper(nid,win,hr)
        blk[:, 0, :] = blk[:, -1, :] = blk[:, :, 0] = blk[:, :, -1] = True
        # sources as flat state indices (native scan) or (L, i, j) rows (numpy fallback, verify mode)
        src = _GR.mask_scan(isl_s, blk) if _GR is not None and not _VERIFY else np.argwhere(isl_s & ~blk)
        tgt = isl_d & ~blk
        if not len(src) or not tgt.any():
            return ('fail', 'entry') if snap else ('fail',)
        lc = [1.0] * NL
        for l, m in (layer_cost or {}).items():
            lc[LAYERS.index(l)] = m
        # negotiated routing (see router.RoutingController): a hook may make some blocked states passable at a
        # cost; the search then runs on this board's (static) obstacles plus that cost.
        soft = self.__dict__.get('soft_cost')
        cost = None
        if soft is not None:
            cost = soft(self, net, nid, win, Ls, half if explicit else None, blk)
        if cost is not None:
            path = self._astar_gr_run(blk, vok, src, tgt, Ls, lc, via_cost / G, pref or {}, turn / G, weight, max_exp,
                                      cost=cost)
            if path is not None:
                self.contested = getattr(self, 'contested', []) + [(L, i + i0, j + j0) for L, i, j in path
                                                                    if cost[L, i, j] > 0]
        else:
            path = self._astar(blk, vok, src, tgt, Ls, lc, via_cost / G, pref or {}, turn / G, weight, max_exp)
        if path is None:
            return ('fail',)
        return ('path', net, path, win, None if self.cls[net] in LAYER_WIDTHS and not explicit else width, (vd, vdr), H, W,
                None if neck is None else [float(neck[k]) for k in path], snap)
        return True

    def _astar(self, blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp):
        nL, H, W = blk.shape
        HW = H * W
        if _GR is not None:
            path = self._astar_gr(blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp)
            if _VERIFY and not _HYBRID:     # the hybrid search differs by design
                assert path == self._astar_ref(blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp), '_astar'
            return path
        blkb = blk.astype(np.uint8).ravel().tobytes()
        vokb = vok.astype(np.uint8).ravel().tobytes()
        tgtb = tgt.astype(np.uint8).ravel().tobytes()
        ti = np.argwhere(tgt)
        ty0, ty1 = int(ti[:, 1].min()), int(ti[:, 1].max())
        tx0, tx1 = int(ti[:, 2].min()), int(ti[:, 2].max())
        hmul = weight * min(lc[L] for L in Ls)
        S2 = math.sqrt(2)
        moves = {}
        for L in Ls:
            p = pref.get(LAYERS[L])
            mv = []
            for d, (di, dj) in enumerate(((0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1))):
                c = S2 if di and dj else 1.0
                if p == 'h':
                    c *= 1.0 if not di else (1.3 if dj else 2.0)
                elif p == 'v':
                    c *= 1.0 if not dj else (1.3 if di else 2.0)
                mv.append((di * W + dj, c * lc[L], d, di, dj))
            moves[L] = mv
        INF = float('inf')
        N = nL * HW
        dist = array('d', [INF]) * N
        prev = array('i', [-1]) * N
        dirn = bytearray(b'\x08') * N
        heap = []
        for L, i, j in src:
            s = int(L) * HW + int(i) * W + int(j)
            dist[s] = 0.0
            heap.append((0.0, 0.0, s))
        heapq.heapify(heap)
        found = None
        n = 0
        push, pop = heapq.heappush, heapq.heappop
        while heap:
            f, g, s = pop(heap)
            if g > dist[s]:
                continue
            if tgtb[s]:
                found = s
                break
            n += 1
            if n > max_exp:
                break
            L, c = divmod(s, HW)
            i, j = divmod(c, W)
            din = dirn[s]
            base = L * HW
            for off, cost, d, di, dj in moves[L]:
                s2 = s + off
                if blkb[s2]:
                    continue
                if di and dj and (blkb[s + dj] or blkb[s + di * W]):
                    continue
                g2 = g + cost
                if din != 8 and din != d:
                    g2 += turn
                if g2 < dist[s2]:
                    dist[s2] = g2
                    prev[s2] = s
                    dirn[s2] = d
                    i2, j2 = i + di, j + dj
                    dx = tx0 - j2 if j2 < tx0 else (j2 - tx1 if j2 > tx1 else 0)
                    dy = ty0 - i2 if i2 < ty0 else (i2 - ty1 if i2 > ty1 else 0)
                    h = (dx + dy - 0.5858 * min(dx, dy)) * hmul
                    push(heap, (g2 + h, g2, s2))
            if vokb[c]:
                for L2 in Ls:
                    if L2 == L:
                        continue
                    s2 = L2 * HW + c
                    if blkb[s2]:
                        continue
                    g2 = g + vcost
                    if g2 < dist[s2]:
                        dist[s2] = g2
                        prev[s2] = s
                        dirn[s2] = 8
                        push(heap, (f - g + g2, g2, s2))
        if found is None:
            return None
        path = [found]
        while prev[path[-1]] >= 0:
            path.append(prev[path[-1]])
        path.reverse()
        return [(s // HW, (s % HW) // W, s % W) for s in path]

    @staticmethod
    def _move_costs(Ls, lc, pref):
        S2 = math.sqrt(2)
        mc = np.zeros((NL, 8), dtype=np.float32)
        for L in Ls:
            p = pref.get(LAYERS[L])
            for d, (di, dj) in enumerate(((0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1))):
                c = S2 if di and dj else 1.0
                # 'h' / 'v': preferred direction; 'H' / 'V': strict (diagonals and cross moves cost much more)
                diag, cross = {'h': (1.3, 2.0), 'v': (1.3, 2.0), 'H': (1.8, 3.5), 'V': (1.8, 3.5)}.get(p, (1.0, 1.0))
                if p in ('h', 'H'):
                    c *= 1.0 if not di else (diag if dj else cross)
                elif p in ('v', 'V'):
                    c *= 1.0 if not dj else (diag if di else cross)
                mc[L, d] = c * lc[L]
        return mc

    def _astar_gr(self, blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp):
        """Search via gridroute; with the route cache, results are also memoised by the search inputs themselves
        (masks, sources, targets, costs), so a search whose surroundings changed only outside its masks is reused."""
        if _CACHE is None or getattr(_TLS, 'spec', False):
            return self._astar_gr_run(blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp)
        if _HYBRID:
            # keyed by the sources, targets and costs; valid while the masks are unchanged on the tiles the search
            # read (its expanded cells): a far-away edit does not invalidate it. Failed searches need the whole
            # masks unchanged (no path is a property of the whole window).
            src32 = np.asarray(src, dtype=np.int32)
            key = 'R' + _GR.hash_arrays([tgt, src32]) + repr((_SALT, blk.shape, Ls, lc, vcost, sorted(pref.items()),
                                                                turn, weight, max_exp))
            whole = lambda: 'W' + _GR.hash_arrays([blk, vok])
            hit = _CACHE.get(key, valid=lambda e: e[1] == (whole() if e[0] is None else _GR.hash_tiles([blk, vok], e[0])))
            if hit is not None:
                return hit[2]
            path, tiles = self._astar_gr_run(blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp, read_set=True)
            ent = (None, whole(), path) if path is None else (tiles, _GR.hash_tiles([blk, vok], tiles), path)
            _CACHE.put(key, ((_CACHE.peek(key) or [])[-3:]) + [ent])
            return path
        key = 'S' + _GR.hash_arrays([blk, vok, tgt, np.asarray(src, dtype=np.int32)]) + repr(
            (_SALT, Ls, lc, vcost, sorted(pref.items()), turn, weight, max_exp))
        hit = _CACHE.get(key)
        if hit is not None:
            return hit[1]
        path = self._astar_gr_run(blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp)
        _CACHE.put(key, ('S', path))
        return path

    def _astar_gr_run(self, blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp, read_set=False, cost=None):
        nL, H, W = blk.shape
        HW = H * W
        lok = np.zeros(NL, dtype=bool)
        lok[Ls] = True
        flat = src if src.ndim == 1 else src[:, 0] * HW + src[:, 1] * W + src[:, 2]
        s = _GR.astar(blk, vok, flat, tgt, self._move_costs(Ls, lc, pref), lok,
                      vcost, turn, weight * min(lc[L] for L in Ls), max_exp,
                      cost=cost, **(dict(budget=_BUDGET, weight=weight, read_set=read_set) if _HYBRID else {}))
        if read_set:
            s, tiles = s
            return (None if s is None else
                    list(zip((s // HW).tolist(), (s % HW // W).tolist(), (s % W).tolist()))), tiles
        if s is None:
            return None
        return list(zip((s // HW).tolist(), (s % HW // W).tolist(), (s % W).tolist()))

    def _astar_ref(self, blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp):
        """gridroute's reference A* (a transcription of the original C search) for verify mode."""
        nL, H, W = blk.shape
        HW = H * W
        lok = np.zeros(NL, dtype=bool)
        lok[Ls] = True
        flat = src if src.ndim == 1 else src[:, 0] * HW + src[:, 1] * W + src[:, 2]
        s = _GR.astar(blk, vok, flat, tgt, self._move_costs(Ls, lc, pref), lok, vcost, turn,
                      weight * min(lc[L] for L in Ls), max_exp, reference=True)
        if s is None:
            return None
        return list(zip((s // HW).tolist(), (s % HW // W).tolist(), (s % W).tolist()))

    def _snap_entries(self, path, coords, win, snap):
        """Move the straight runs at the ends of a path onto their pad's exact centre line and return the pad
        centres to start / end the track at (None where the end is not on an axial pad)."""
        i0, j0 = win[0], win[1]
        ends = [None, None]
        for end, order in ((0, range(len(path))), (1, range(len(path) - 1, -1, -1))):
            order = list(order)
            L, i, j = path[order[0]]
            gi, gj = i + i0, j + j0
            for x, y, ci, cj, reach in snap:
                if (gi != ci and gj != cj) or math.hypot(gj * G - x, gi * G - y) > reach:
                    continue
                nxt = path[order[1]] if len(order) > 1 else None
                horizontal = gi == ci and (gj != cj or (nxt is not None and nxt[1] + i0 == gi and nxt[0] == L))
                for k in order:
                    s = path[k]
                    if s[0] != L or (s[1] + i0 != gi if horizontal else s[2] + j0 != gj):
                        break
                    coords[k] = (coords[k][0], y) if horizontal else (x, coords[k][1])
                ends[end] = (x, y)
                break
        return ends

    def _emit(self, net, path, win, width, via, H, W, neck=None, snap=None):
        i0, j0 = win[0], win[1]
        coords = [((k[2] + j0) * G, (k[1] + i0) * G) for k in path]
        ends = self._snap_entries(path, coords, win, snap) if snap else [None, None]
        self._emit_ends = (ends[0], ends[1], len(path) - 1)
        path = [(k[0], k[1], k[2], n) for n, k in enumerate(path)]
        xy = lambda k: coords[k[3]]
        neck = neck or [0.] * len(path)
        seg, nw = [path[0]], [neck[0]]
        for k, w in zip(path[1:], neck[1:]):
            if k[0] != seg[-1][0]:
                self._flush_necked(net, seg, nw, xy, width)
                self.add_via(net, *xy(k), d=via[0], drill=via[1])
                seg, nw = [k], [w]
            else:
                seg.append(k)
                nw.append(w)
        self._flush_necked(net, seg, nw, xy, width)

    def _flush_necked(self, net, seg, nw, xy, width):
        """Flush one layer's run, thin where it lies in a neck region: a thin run extends one cell into the
        class-width copper on each side (so it still touches the region), class-width runs cover the rest."""
        if not any(nw):
            return self._flush(net, seg, xy, width)
        k = 0
        while k < len(seg):
            e = k
            while e + 1 < len(seg) and nw[e + 1] == nw[k]:
                e += 1
            if nw[k]:
                self._flush(net, seg[max(k - 1, 0):e + 2], xy, nw[k])
            else:
                self._flush(net, seg[k:e + 1], xy, width)
            k = e + 1

    def _flush(self, net, seg, xy, width):
        if len(seg) < 2:
            return
        pts = [seg[0]]
        for a, b, c in zip(seg, seg[1:], seg[2:]):
            if (b[1] - a[1], b[2] - a[2]) != (c[1] - b[1], c[2] - b[2]):
                pts.append(b)
        pts.append(seg[-1])
        out = [xy(k) for k in pts]
        start, end, last = getattr(self, '_emit_ends', (None, None, -1))
        if len(seg[0]) > 3:
            if start is not None and seg[0][3] == 0 and math.dist(start, out[0]) > 1e-9:
                out.insert(0, start)
            if end is not None and seg[-1][3] == last and math.dist(end, out[-1]) > 1e-9:
                out.append(end)
            out = _collinear_merged(out)
        self.add_track(net, LAYERS[seg[0][0]], out, width)

    # ------------------------------------------------------------------ differential pairs
    def route_pair(self, netp, netn, a, b, layer, width=None, gap=0.2, margin=4.0, turn=0.3, pref=None,
                   max_exp=4_000_000, window=None, a_dir=None, b_dir=None, lead=1.0, avoid=None):
        """Route a coupled pair on one layer. a, b: ((xp, yp), (xn, yn)), the P and N points where the pair starts
        and ends on `layer` (pad ends or via centres already drawn). The centreline is searched as one track of
        width 2w + gap that may cross either net's copper, then offset into the two tracks. Both ends get a straight
        lead of length `lead` square to the end pair (a_dir / b_dir: travel direction there, default towards the
        other end), so the tracks meet their pads or vias head-on. The pair must not need to swap sides between a
        and b. avoid: [(layer, allowed nets)]: also keep the pair clear of other copper on that layer (an adjacent
        signal layer that would otherwise sit under the pair), except the allowed nets, which it may cross.
        Returns True on success."""
        w = width or self.width(netp, layer)
        half = w + gap / 2                          # centreline to outer edge
        off = (w + gap) / 2                         # centreline to each track centre
        ids = (self.net_id[netp], self.net_id[netn])
        L = LAYERS.index(layer)
        ma = ((a[0][0] + a[1][0]) / 2, (a[0][1] + a[1][1]) / 2)
        mb = ((b[0][0] + b[1][0]) / 2, (b[0][1] + b[1][1]) / 2)

        def square(e, m, d, sign):
            ax, ay = e[0][0] - e[1][0], e[0][1] - e[1][1]
            l = math.hypot(ax, ay)
            n = (-ay / l, ax / l)
            if d is None:
                d = n if (n[0] * (mb[0] - ma[0]) + n[1] * (mb[1] - ma[1])) * sign > 0 else (-n[0], -n[1])
            return d
        da = square(a, ma, a_dir, 1)
        db = square(b, mb, b_dir, 1)
        la = (ma[0] + da[0] * lead, ma[1] + da[1] * lead)          # centreline leaves a along da ...
        lb = (mb[0] - db[0] * lead, mb[1] - db[1] * lead)          # ... and arrives at b along db
        xs = [ma[0], mb[0]]
        ys = [ma[1], mb[1]]
        x0, y0, x1, y1 = window or (min(xs) - margin, min(ys) - margin, max(xs) + margin, max(ys) + margin)
        win = self._win(x0, y0, x1, y1)
        i0, j0, i1, j1 = win
        H, W = i1 - i0 + 1, j1 - j0 + 1
        r = half + CLEAR + self.extra(netp) + MARGIN
        blk = np.ones((NL, H, W), dtype=bool)
        blk[L] = self._blocked(L, ids, win, r)
        for av_layer, allow in avoid or ():
            blk[L] |= self._blocked(LAYERS.index(av_layer), ids + tuple(self.net_id[n] for n in allow), win, half + 0.1)
        blk[:, 0, :] = blk[:, -1, :] = blk[:, :, 0] = blk[:, :, -1] = True
        # the pair may only reach its ends along the leads: keep the centreline away from both end pairs (the pair's
        # own pads / vias do not block it otherwise), then free the lead points
        ys_, xs_ = self._grid(win)
        for m, l in ((ma, la), (mb, lb)):
            d2 = (xs_ - m[0]) ** 2 + (ys_ - m[1]) ** 2
            blk[L] |= d2 < (lead - 0.05) ** 2
        for m in (la, lb):
            cw, cm = self._mask_seg(m[0], m[1], m[0], m[1], off + 0.05)
            sub = blk[L, cw[0] - i0:cw[2] - i0 + 1, cw[1] - j0:cw[3] - j0 + 1]
            sub[cm[:sub.shape[0], :sub.shape[1]]] = False
        cell = lambda m: (L, int(round(m[1] / G)) - i0, int(round(m[0] / G)) - j0)
        src = np.array([cell(la)])
        tgt = np.zeros_like(blk)
        tgt[cell(lb)] = True
        lc = [1.0] * NL
        path = self._astar(blk, np.zeros((H, W), dtype=bool), src, tgt, [L], lc, 1e9, pref or {}, turn / G, 1.2,
                           max_exp)
        if path is None:
            self.failed.append(('pair', netp, netn))
            return False
        pts = [path[0]]
        for u, v, q in zip(path, path[1:], path[2:]):
            if (v[1] - u[1], v[2] - u[2]) != (q[1] - v[1], q[2] - v[2]):
                pts.append(v)
        pts.append(path[-1])
        cl = [((k[2] + j0) * G, (k[1] + i0) * G) for k in pts]
        cl[0], cl[-1] = la, lb
        cl = dedup([ma] + cl + [mb])
        # which side is P: sign of the cross product of the first direction with (P - mid)
        dx, dy = da
        side = 1.0 if (a[0][0] - ma[0]) * dy - (a[0][1] - ma[1]) * dx > 0 else -1.0     # P left of travel?
        tp = [a[0]] + offset_polyline(cl, side * off) + [b[0]]
        tn = [a[1]] + offset_polyline(cl, -side * off) + [b[1]]
        self.add_track(netp, layer, dedup(tp), w)
        self.add_track(netn, layer, dedup(tn), w)
        return True

    def keepout_along(self, tracks, layer, r):
        """Router keep-out on `layer` within r of the given tracks (e.g. an adjacent signal layer kept clear under a
        stripline pair, so the pair sees its two planes)."""
        self._op('along', ([[tuple(q) for q in t['pts']] for t in tracks], LAYERS.index(layer), r))

    def _do_along(self, tracks, L, r):
        for t in tracks:
            for a, b in zip(t, t[1:]):
                self._paint_seg(self.occ, L, a, b, r, KEEP)

    def twin_via(self, net, x, y, dist=0.85, width=0.6, layers=('F.Cu', 'B.Cu')):
        """Put a second via of the net next to the via at (x, y), linked to it on each of `layers`, so a power spur
        carries its current through two holes. Tries eight directions; returns the twin's position or None."""
        nid = self.net_id[net]
        vd, vdr = self.via_size(net)
        for k in range(8):
            a = k * math.pi / 4
            cx, cy = x + dist * math.cos(a), y + dist * math.sin(a)
            win = self._win(cx - 1, cy - 1, cx + 1, cy + 1)
            i, j = int(round(cy / G)) - win[0], int(round(cx / G)) - win[1]
            ok = all(not self._blocked(L, nid, win, self._reach(net, vd / 2))[i, j] for L in range(NL))
            ok = ok and not self._near_smd(win, vd / 2 + 0.1)[i, j]
            ok = ok and all(math.hypot(v['x'] - cx, v['y'] - cy) >= (v['drill'] + vdr) / 2 + 0.3
                            for v in self.vias if abs(v['x'] - cx) < 2 and abs(v['y'] - cy) < 2 and (v['x'], v['y']) != (x, y))
            if not ok:
                continue
            n = 8
            for ly in layers:
                L = LAYERS.index(ly)
                sw = self._win(min(x, cx) - 1, min(y, cy) - 1, max(x, cx) + 1, max(y, cy) + 1)
                blk = self._blocked(L, nid, sw, self._reach(net, width / 2))
                pts = [(x + (cx - x) * t / n, y + (cy - y) * t / n) for t in range(n + 1)]
                if any(blk[int(round(py / G)) - sw[0], int(round(px / G)) - sw[1]] for px, py in pts):
                    ok = False
                    break
            if not ok:
                continue
            self.add_via(net, cx, cy)
            for ly in layers:
                self.add_track(net, ly, [(x, y), (cx, cy)], width)
            return cx, cy
        return None

    # ------------------------------------------------------------------ plane drops
    def fanout(self, ref, num, max_r=1.5, dirs=None, width=0.25, in_pad=False):
        """Plane drop for an SMD pad (see _fanout); with the route cache, replayed when the board around the pad is
        unchanged."""
        p = self.parts[ref].pad(num)
        if _CACHE is None or p['layers'] == 'all' or self.has_drop(p):
            return self._fanout(ref, num, max_r, dirs, width, in_pad)
        r = max(4.0, max_r + 2.5) + max(p['w'], p['h'])
        key = self._content_key(self._win(p['x'] - r, p['y'] - r, p['x'] + r, p['y'] + r), 'fanout', ref, num,
                                (p['x'], p['y'], p['w'], p['h'], p['net']), max_r, dirs, width, in_pad)
        hit = _CACHE.get(key)
        if hit is not None:
            ok, ops, failed = hit
            for kind, args, drill in ops:
                if kind == 'track':
                    self.add_track(args[0], LAYERS[args[1]], args[2], args[3])
                else:
                    self.add_via(args[0], args[1], args[2], d=args[3], drill=drill)
            self.failed += failed
            return ok
        k0, nf = len(self.ops), len(self.failed)
        ok = self._fanout(ref, num, max_r, dirs, width, in_pad)
        _CACHE.put(key, (ok, [(kind, args, item.get('drill')) for kind, args, w, item in self.ops[k0:]],
                         self.failed[nf:]))
        return ok

    def _fanout(self, ref, num, max_r=1.5, dirs=None, width=0.25, in_pad=False):
        """Drop a via to the plane next to an SMD pad of a plane net (see configure(planes=...)) with a short stub.

        dirs: preferred stub directions as unit vectors (default: any of the eight 45-degree directions).
        in_pad: allow the via inside the pad (filled via-in-pad).
        """
        p = self.parts[ref].pad(num)
        net = p['net']
        if p['layers'] == 'all' or self.has_drop(p):
            return True
        nid = self.net_id[net]
        L = pad_layers(p)[0]
        vd, vdr = self.via_size(net)
        win = self._win(p['x'] - max_r - 1, p['y'] - max_r - 1, p['x'] + max_r + 1, p['y'] + max_r + 1)
        i0, j0 = win[0], win[1]
        vok = np.ones((win[2] - i0 + 1, win[3] - j0 + 1), dtype=bool)
        for LL in range(NL):
            vok &= ~self._blocked(LL, nid, win, self._reach(net, vd / 2, LAYERS[LL]))
        if not in_pad:
            vok &= ~self._near_smd(win, vd / 2 + 0.1)
        vok &= ~self._near_holes(win, vdr / 2 + HOLE_TO_HOLE + 0.01)
        vok &= ~self._via_areas(win, vd / 2 + MARGIN)
        hr = vdr / 2 + HOLE_CLEAR + MARGIN
        if hr > self._reach(net, vd / 2):
            vok &= ~self._near_foreign_copper(nid,win,hr)
        sok = ~self._blocked(L, nid, win, self._reach(net, width / 2, LAYERS[L]))
        # an existing via of the net within reach: a straight stub to it instead of another hole
        for v in sorted((v for v in self.vias if v['net'] == net),
                        key=lambda v: math.hypot(v['x'] - p['x'], v['y'] - p['y'])):
            d = math.hypot(v['x'] - p['x'], v['y'] - p['y'])
            if d > max_r + max(p['w'], p['h']) / 2:
                break
            n = max(int(d / G), 1)
            if all(sok[int(round((p['y'] + (v['y'] - p['y']) * t / n) / G)) - i0,
                       int(round((p['x'] + (v['x'] - p['x']) * t / n) / G)) - j0] for t in range(n + 1)
                   if 0 <= int(round((p['y'] + (v['y'] - p['y']) * t / n) / G)) - i0 < sok.shape[0]
                   and 0 <= int(round((p['x'] + (v['x'] - p['x']) * t / n) / G)) - j0 < sok.shape[1]):
                self.add_track(net, LAYERS[L], [(p['x'], p['y']), (v['x'], v['y'])], width)
                return True
        # _near_holes above already includes every existing via and pad drill
        # using the configured hole-to-hole rule, including same-net holes.
        dirs = dirs or [(math.cos(a * math.pi / 4), math.sin(a * math.pi / 4)) for a in range(8)]
        best = None
        rmin = 0.0 if in_pad else max(p['w'], p['h']) / 2 + 0.05
        for dx, dy in dirs:
            steps = int(max_r / G)
            for k in range(int(rmin / G), steps + 1):
                x, y = p['x'] + dx * k * G, p['y'] + dy * k * G
                i, j = int(round(y / G)) - i0, int(round(x / G)) - j0
                if not (0 <= i < vok.shape[0] and 0 <= j < vok.shape[1]) or not vok[i, j]:
                    continue
                good = True
                # Test the segment that will actually be emitted: its via end
                # is snapped to the grid and may differ from the nominal ray.
                ex,ey=(j+j0)*G,(i+i0)*G
                for t in range(k + 1):
                    ratio=t/max(k,1)
                    ii = int(round((p['y'] + (ey-p['y'])*ratio) / G)) - i0
                    jj = int(round((p['x'] + (ex-p['x'])*ratio) / G)) - j0
                    if not sok[ii, jj]:
                        good = False
                        break
                if good:
                    d = k * G
                    if best is None or d < best[0]:
                        best = (d, i + i0, j + j0)
                    break
        if not best:
            self.failed.append((net, ref, num))
            return False
        _, i, j = best
        x, y = j * G, i * G
        if math.hypot(x - p['x'], y - p['y']) > 1e-6:
            self.add_track(net, LAYERS[L], [(p['x'], p['y']), (x, y)], width)
        self.add_via(net, x, y)
        return True

    def via_ok(self, net, win):
        """[H, W] cells of grid window win where a via of `net` fits against the current board (the same tests as
        a plane drop: copper on every layer, SMD pads, holes, via rule areas, drill-to-copper clearance)."""
        nid = self.net_id[net]
        vd, vdr = self.via_size(net)
        vok = np.ones((win[2] - win[0] + 1, win[3] - win[1] + 1), dtype=bool)
        for L in range(NL):
            vok &= ~self._blocked(L, nid, win, self._reach(net, vd / 2, LAYERS[L]))
        vok &= ~self._near_smd(win, vd / 2 + 0.1)
        vok &= ~self._near_holes(win, vdr / 2 + HOLE_TO_HOLE + 0.01)
        vok &= ~self._via_areas(win, vd / 2 + MARGIN)
        hr = vdr / 2 + HOLE_CLEAR + MARGIN
        if hr > self._reach(net, vd / 2):
            vok &= ~self._near_foreign_copper(nid, win, hr)
        return vok

    def track_ok(self, net, L, width, win):
        """[H, W] cells of window win where a track centre of `net` with `width` fits on layer index L."""
        return ~self._blocked(L, self.net_id[net], win, self._reach(net, width / 2, LAYERS[L]))

    def plane_reached(self, p):
        """A plane net's SMD pad reaches a via or plated hole of its net through copper: its own drop (3 mm), or
        a trace tied to a neighbour's drop (8 mm)."""
        return self.has_drop(p, 3.) or self.has_drop(p, 8.)

    def has_drop(self, p, reach=4.0):
        """True if pad p already reaches a via or plated hole of its own net through copper (remembered: copper is
        only removed by rip(), which forgets)."""
        known = self.__dict__.setdefault('_dropped', set())
        k = (p['ref'], p['num'], reach)
        if k in known:
            return True
        if self._has_drop(p, reach):
            known.add(k)
            return True
        return False

    def _has_drop(self, p, reach):
        nid = self.net_id[p['net']]
        w = self._win(p['x'] - reach, p['y'] - reach, p['x'] + reach, p['y'] + reach)
        isl = self._island(nid, self._pad_seeds(p, w), w)
        thr = self.thru[w[0]:w[2] + 1, w[1]:w[3] + 1] == nid
        return bool((isl.any(axis=0) & thr).any())

    def fanout_all(self, refs, nets=None, **kw):
        nets = tuple(PLANES) if nets is None else nets
        for r in refs:
            for p in self.parts[r].pads():
                if p['net'] in nets and p['layers'] != 'all':
                    if not self.fanout(r, p['num'], **kw) and 'max_r' not in kw:
                        self.fanout(r, p['num'], max_r=3.0, **kw)     # crowded spot: a longer stub

    # ------------------------------------------------------------------ rip-up
    def nets_in(self, win, since=0):
        """Nets with tracks or vias logged since op index `since` painted inside grid window win."""
        out = set()
        for kind, args, w, item in self.ops[since:]:
            if item is not None and w[0] <= win[2] and win[0] <= w[2] and w[1] <= win[3] and win[1] <= w[3]:
                out.add(item['net'])
        return out

    def rip(self, nets, since=0):
        """Remove the tracks and vias of `nets` logged since op index `since`: clear the region they painted and
        replay the remaining operations that touch it, in their original order (so the rasters are exactly what
        they would be had the copper never been added). Returns the number of items removed."""
        nets = set(nets)
        gone = [k for k in range(since, len(self.ops)) if self.ops[k][3] is not None and self.ops[k][3]['net'] in nets]
        return self._rip_ops(gone)

    def _rip_ops(self, gone):
        """Remove the logged operations with indices `gone` (tracks / vias) and repaint their region."""
        self.spec_clear()
        self._dropped = set()
        if not gone:
            return 0
        i0 = min(self.ops[k][2][0] for k in gone)
        j0 = min(self.ops[k][2][1] for k in gone)
        i1 = max(self.ops[k][2][2] for k in gone)
        j1 = max(self.ops[k][2][3] for k in gone)
        ids = {id(self.ops[k][3]) for k in gone}
        drop = set(gone)
        self.ops = [op for k, op in enumerate(self.ops) if k not in drop]
        self.tracks = [t for t in self.tracks if id(t) not in ids]
        self.vias = [v for v in self.vias if id(v) not in ids]
        for a in (self.occ, self.core) + ((self.rsv,) if self.__dict__.get('rsv') is not None else ()):
            a[:, i0:i1 + 1, j0:j1 + 1] = 0
        for a in (self.thru, self.smd, self.no_via):
            a[i0:i1 + 1, j0:j1 + 1] = 0
        self._dirty[i0 // TILE:i1 // TILE + 1, j0 // TILE:j1 // TILE + 1] = True
        # replay into a scratch board region: each operation repaints in full (outside the region that is idempotent)
        self._edges()
        for kind, args, w, item in self.ops:
            if w[0] <= i1 and i0 <= w[2] and w[1] <= j1 and j0 <= w[3]:
                getattr(self, '_do_' + kind)(*args)
        return len(gone)

    # ------------------------------------------------------------------ relaxation
    @staticmethod
    def _cost(items, via_cost):
        """Track length (mm) plus via_cost mm per via of a list of track / via dicts."""
        return sum(math.hypot(b[0] - a[0], b[1] - a[1]) for t in items if 'pts' in t for a, b in zip(t['pts'], t['pts'][1:])) + \
            via_cost * sum(1 for v in items if 'pts' not in v)

    def _net_complete(self, net):
        pads = self.pads_of(net)
        if len(pads) < 2:
            return True
        return len(set(self.pad_components(net, pads))) <= 1

    def replace_copper(self, tracks, vias):
        """Replace routed copper and synchronize all rasters and operation logs.

        Footprints, keepouts and other static geometry remain in place. Inputs
        may be the board's current lists; new dictionaries are made on replay.
        This is also the rollback primitive for an external routing controller.
        """
        tracks, vias = list(tracks), list(vias)
        self._rip_ops([i for i, op in enumerate(self.ops) if op[3] is not None])
        self.tracks, self.vias = [], []
        for t in tracks:
            self.add_track(t['net'], t['layer'], t['pts'], t['width'])
        for v in vias:
            self.add_via(v['net'], v['x'], v['y'], d=v['d'], drill=v['drill'])

    def relax(self, nets=None, since=0, passes=3, via_cost=1.0, turn=0.05, margin=3.0, gain=0.2, layers=None):
        """Rip-up-and-reroute relaxation, run once the whole board is routed.

        Routes found early are shaped by copper that later moved or never came, by layer direction preferences and
        by the weighted / field-guided search; the result is detours. Here every net with copper logged since op
        index `since` is ripped (that copper only: hand-drawn breakouts and anything before `since` stay) and
        re-routed against the finished board with an unweighted search and no direction preferences. The new route
        is kept if the net is complete and it is at least `gain` mm cheaper (track length plus `via_cost` mm per
        via); otherwise the old copper goes back. Plane nets are skipped. Passes repeat, in net order, until
        one improves nothing. Returns [(net, mm saved)]."""
        if nets is None:
            nets = sorted({op[3]['net'] for op in self.ops[since:] if op[3] is not None})
        nets = [n for n in nets if n not in PLANES]
        saved = []
        for _ in range(passes):
            better = 0
            for n in nets:
                old = [op[3] for op in self.ops[since:] if op[3] is not None and op[3]['net'] == n]
                if not old:
                    continue
                before = self._cost(old, via_cost)
                nf = len(self.failed)
                self.rip([n], since)
                k0 = len(self.ops)
                ok = self.connect(n, layers=layers, strict=layers is not None, via_cost=via_cost, turn=turn,
                                  margin=margin, weight=1.0)
                new = [op[3] for op in self.ops[k0:] if op[3] is not None and op[3]['net'] == n]
                after = self._cost(new, via_cost)
                del self.failed[nf:]
                if ok and after < before - gain and self._net_complete(n):
                    saved.append((n, before - after))
                    better += 1
                    continue
                self.rip([n], k0)                               # put the old copper back
                for it in old:
                    if 'pts' in it:
                        self.add_track(it['net'], it['layer'], [tuple(p) for p in it['pts']], it['width'])
                    else:
                        self.add_via(it['net'], it['x'], it['y'], d=it['d'], drill=it['drill'])
            if not better:
                break
        return saved

    def pull_tight(self, since=0, nets=None, gain=0.05):
        """Octilinear string pulling of the tracks logged since op index `since`.

        Within each track polyline, a run of vertices between two of its vertices is replaced by a straight
        octilinear link or a one-bend link (diagonal then straight, or straight then diagonal) when that link is
        clear of other nets' copper and keep-outs at the track's width and clearance, at least `gain` mm shorter,
        and no other copper of the net (track ends, vias, pads) attaches to the replaced run. Repeats until nothing
        shortens. Returns the mm saved."""
        total = 0.0
        while True:
            changed = False
            for k in range(since, len(self.ops)):
                if k >= len(self.ops):
                    break
                kind, args, w, t = self.ops[k]
                if kind != 'track' or t is None or len(t['pts']) < 3 or (nets and t['net'] not in nets):
                    continue
                pts = self._pull(t, gain)
                if pts is None:
                    continue
                old = self._cost([t], 0.0)
                self._rip_ops([k])
                self.add_track(t['net'], t['layer'], pts, t['width'])
                total += old - self._cost([{'pts': pts}], 0.0)
                changed = True
            if not changed:
                return total

    def _anchors(self, t):
        """Points where other copper of t's net may attach: other tracks' ends and every vertex, vias, pads."""
        net, out = t['net'], []
        for o in self.tracks:
            if o['net'] == net and o is not t:
                out += [(p[0], p[1], o['width'] / 2) for p in o['pts']]
        out += [(v['x'], v['y'], v['d'] / 2) for v in self.vias if v['net'] == net]
        out += [(p['x'], p['y'], max(p['w'], p['h']) / 2) for p in self.pads_of(net)]
        return out

    def _pull(self, t, gain=0.05):
        """A shorter vertex list for track t (see pull_tight), or None."""
        pts = [tuple(p) for p in t['pts']]
        half = t['width'] / 2
        anchors = self._anchors(t)
        # which segments each anchor touches; a run may only be replaced if every anchor touching it sits at one
        # of the run's two kept vertices (the pad or via the track starts / ends on, a branch at a kept vertex)
        touch = [(x, y, r, {s for s, (a, b) in enumerate(zip(pts, pts[1:])) if _seg_dist(x, y, a, b) < r + half + 0.02})
                 for x, y, r in anchors]
        touch = [a for a in touch if a[3]]
        L = LAYERS.index(t['layer'])
        nid = self.net_id[t['net']]
        reach = self._reach(t['net'], half, t['layer'])

        def free(i, j):
            for x, y, r, segs in touch:
                if any(i <= s < j for s in segs) and not any(math.hypot(x - pts[k][0], y - pts[k][1]) < r + half + 0.02
                                                              for k in (i, j)):
                    return False
            return True
        for i in range(len(pts) - 2):
            for j in range(len(pts) - 1, i + 1, -1):
                if not free(i, j):
                    continue
                run = sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(pts[i:j], pts[i + 1:j + 1]))
                for link in self._links(pts[i], pts[j]):
                    ln = sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(link, link[1:]))
                    if ln < run - gain and self._clear(link, L, nid, reach):
                        return pts[:i] + link + pts[j + 1:]
        return None

    @staticmethod
    def _links(p, q):
        """Octilinear links from p to q: straight if p-q is at 0/45/90 degrees, else the two one-bend links."""
        dx, dy = q[0] - p[0], q[1] - p[1]
        ax, ay = abs(dx), abs(dy)
        if ax < 1e-6 or ay < 1e-6 or abs(ax - ay) < 1e-6:
            return [[p, q]]
        m = min(ax, ay)
        sx, sy = math.copysign(1, dx), math.copysign(1, dy)
        c1 = (p[0] + sx * m, p[1] + sy * m)
        c2 = (q[0] - sx * m, q[1] - sy * m)
        return [[p, c1, q], [p, c2, q]]

    def _clear(self, link, L, nid, reach):
        xs = [p[0] for p in link]
        ys = [p[1] for p in link]
        win = self._win(min(xs) - 1, min(ys) - 1, max(xs) + 1, max(ys) + 1)
        blk = self._blocked(L, nid, win, reach)
        for a, b in zip(link, link[1:]):
            n = max(int(math.hypot(b[0] - a[0], b[1] - a[1]) / (G / 2)), 1)
            for s in range(n + 1):
                x, y = a[0] + (b[0] - a[0]) * s / n, a[1] + (b[1] - a[1]) * s / n
                i, j = int(round(y / G)) - win[0], int(round(x / G)) - win[1]
                if not (0 <= i < blk.shape[0] and 0 <= j < blk.shape[1]) or blk[i, j]:
                    return False
        return True

    # ------------------------------------------------------------------ replication
    def mark(self):
        return len(self.tracks), len(self.vias)

    def since(self, mark):
        return self.tracks[mark[0]:], self.vias[mark[1]:]

    def net_map(self, ref_map):
        """Nets of the source refs -> nets of the mapped refs, from pad correspondences. A source net whose pads map
        to different nets (a strap tied to +3V3 in one copy and GND in another) takes the majority; copper touching the
        odd pad must be kept out of the template."""
        votes = {}
        for a, b in ref_map.items():
            for (r, pin), n in self.pin_net.items():
                if r == a:
                    t = self.pin_net.get((b, pin))
                    votes.setdefault(n, {}).setdefault(t, 0)
                    votes[n][t] += 1
        return {n: max(v, key=lambda t: (v[t], t == n)) for n, v in votes.items()}

    def replicate(self, ref_map, xf, drot, items):
        """Place the mapped parts and copy tracks/vias through the rigid transform xf (x, y) -> (x', y')."""
        nm = self.net_map(ref_map)
        for a, b in ref_map.items():
            pa = self.parts[a]
            if pa.placed and not self.parts[b].placed:
                self.place(b, *xf(pa.x, pa.y), pa.rot + drot)
        tracks, vias = items
        for t in list(tracks):
            self.add_track(nm[t['net']], t['layer'], [xf(*p) for p in t['pts']], t['width'])
        for v in list(vias):
            self.add_via(nm[v['net']], *xf(v['x'], v['y']), d=v['d'], drill=v['drill'])

    def prune_vias(self):
        """Remove one-layer signal vias (see _prune_vias); memoised on the copper and placement."""
        def run():
            n = self._prune_vias()
            return self.tracks, self.vias, n
        before_tracks, before_vias = list(self.tracks), list(self.vias)
        before_missing = set(self.unrouted())
        self.tracks, self.vias, n = self._memo('prune', (self.tracks, self.vias, self._placement()), run)
        self.replace_copper(self.tracks, self.vias)
        # A via's copper disc can bridge disjoint track ends even on one layer.
        # Do not silently destroy that connection during a cosmetic cleanup.
        if set(self.unrouted()) != before_missing:
            self.replace_copper(before_tracks, before_vias)
            return 0
        return n

    def via_bridges_copper(self, v):
        """Whether deleting a via can break a physical copper joint.

        Raster flood-fill joins adjacent pixels even when physical track ends
        are only tangent. Require positive overlap between all local branches.
        Pad/via contacts are conservatively retained by this cosmetic cleanup.
        """
        eps=0.002
        ts=[t for t in self.tracks if t['net']==v['net'] and any(
            _seg_dist(v['x'],v['y'],a,b)<(v['d']+t['width'])/2+eps
            for a,b in zip(t['pts'],t['pts'][1:]))]
        if any(_pad_distance(v['x'],v['y'],p)<v['d']/2+eps for p in self.pads_of(v['net'])):
            return True
        if any(u is not v and u['net']==v['net'] and math.hypot(u['x']-v['x'],u['y']-v['y'])<(u['d']+v['d'])/2+eps for u in self.vias):
            return True
        if len(ts)<2:return False
        def distance(a,b,c,d):
            def orient(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
            if orient(a,b,c)*orient(a,b,d)<0 and orient(c,d,a)*orient(c,d,b)<0:return 0.
            return min(_seg_dist(*a,c,d),_seg_dist(*b,c,d),_seg_dist(*c,a,b),_seg_dist(*d,a,b))
        seen={0};todo=[0]
        while todo:
            i=todo.pop();t=ts[i]
            for j,u in enumerate(ts):
                if j in seen or u['layer']!=t['layer']:continue
                if any(distance(a,b,c,d)<(t['width']+u['width'])/2-eps
                       for a,b in zip(t['pts'],t['pts'][1:]) for c,d in zip(u['pts'],u['pts'][1:])):
                    seen.add(j);todo.append(j)
        return len(seen)!=len(ts)

    def via_replacement_tracks(self,v):
        """Replace a one-layer via joint inside the existing copper union.

        Returns None when any contact cannot be joined without adding copper
        outside the old via disc and touching tracks. The caller retains those vias.
        """
        if any(_pad_distance(v['x'],v['y'],p)<v['d']/2+.002 for p in self.pads_of(v['net'])):return None
        if any(u is not v and u['net']==v['net'] and math.hypot(u['x']-v['x'],u['y']-v['y'])<(u['d']+v['d'])/2+.002 for u in self.vias):return None
        stars=[];layers=set()
        for t in self.tracks:
            if t['net']!=v['net']:continue
            closest=None
            for a,b in zip(t['pts'],t['pts'][1:]):
                dx,dy=b[0]-a[0],b[1]-a[1];den=dx*dx+dy*dy
                k=max(0.,min(1.,((v['x']-a[0])*dx+(v['y']-a[1])*dy)/den)) if den else 0.
                q=(a[0]+k*dx,a[1]+k*dy);d=math.hypot(q[0]-v['x'],q[1]-v['y'])
                if closest is None or d<closest[0]:closest=(d,q)
            if closest is None or closest[0]>=(v['d']+t['width'])/2+.002:continue
            layers.add(t['layer']);width=self.width(v['net'],t['layer'])
            # The cylinder up to the old track's centre must lie in the via;
            # the final round cap is already covered by that existing track.
            if width>t['width']+1e-9 or math.hypot(closest[0],width/2)>v['d']/2-.001:return None
            if closest[0]>.0001:stars.append(dict(net=v['net'],layer=t['layer'],width=width,pts=[(v['x'],v['y']),closest[1]]))
        return stars if len(layers)==1 else None

    def _prune_vias(self):
        """Remove signal vias that touch copper on one layer only (e.g. a template's branch via a copy did not use).
        Tracks meeting at such a via stay connected to each other (same layer, overlapping ends); a track that only
        led to it is a dead-end spur and goes too. Plane-net vias always reach their plane and are kept. Repeats
        until nothing changes; returns the number of vias removed."""
        removed = 0
        while True:
            by_net = {}
            for ti, t in enumerate(self.tracks):
                by_net.setdefault(t['net'], []).append(ti)
            pads = {}
            for q in self.parts.values():
                if q.placed:
                    for p in q.pads():
                        if p['net']:
                            pads.setdefault(p['net'], []).append(p)
            drop_v, drop_t = set(), set()
            for vi, v in enumerate(self.vias):
                if v['net'] in PLANES:
                    continue
                touch = []
                for ti in by_net.get(v['net'], []):
                    t = self.tracks[ti]
                    d = min(_seg_dist(v['x'], v['y'], a, b) for a, b in zip(t['pts'], t['pts'][1:]))
                    if d < v['d'] / 2 + t['width'] / 2 - 1e-3:
                        touch.append(ti)
                layers = {self.tracks[ti]['layer'] for ti in touch}
                at_pad = False
                for p in pads.get(v['net'], []):
                    if _pad_distance(v['x'], v['y'], p) < v['d'] / 2 - 1e-3:
                        layers |= {LAYERS[L] for L in pad_layers(p)}
                        at_pad = True
                if len(layers) >= 2 or self.via_bridges_copper(v):
                    continue
                drop_v.add(vi)
                if len(touch) == 1 and not at_pad:
                    t = self.tracks[touch[0]]
                    ends = [e for e in (t['pts'][0], t['pts'][-1]) if math.hypot(e[0] - v['x'], e[1] - v['y']) < v['d'] / 2]
                    others = [tj for tj in by_net[v['net']] if tj != touch[0] and ends and any(
                        _seg_dist(ends[0][0], ends[0][1], a, b) < 1e-3 for a, b in zip(self.tracks[tj]['pts'], self.tracks[tj]['pts'][1:]))]
                    if ends and not others and not layers - {t['layer']}:
                        drop_t.add(touch[0])
            if not drop_v:
                return removed
            removed += len(drop_v)
            self.vias = [v for i, v in enumerate(self.vias) if i not in drop_v]
            self.tracks = [t for i, t in enumerate(self.tracks) if i not in drop_t]

    def trim_dangling(self):
        """Cut back free track ends (see _trim_dangling); memoised on the copper and placement."""
        def run():
            n = self._trim_dangling()
            return self.tracks, n
        self.tracks, n = self._memo('trim', (self.tracks, self.vias, self._placement()), run)
        self.replace_copper(self.tracks, self.vias)
        return n

    def _trim_dangling(self):
        """Cut track polylines back from free ends (an end touching no pad, via or other track of its net) to the
        last point where other copper of the net touches them; tracks touching nothing else go entirely.
        Returns the length removed (mm)."""
        removed = 0.0
        pads = {}
        for q in self.parts.values():
            if q.placed:
                for p in q.pads():
                    if p['net']:
                        pads.setdefault(p['net'], []).append(p)
        changed = True
        while changed:
            changed = False
            by_net = {}
            for ti, t in enumerate(self.tracks):
                by_net.setdefault(t['net'], []).append(ti)
            vias = {}
            for v in self.vias:
                vias.setdefault(v['net'], []).append(v)

            def touched(ti, x, y, layer, w):
                n = self.tracks[ti]['net']
                for p in pads.get(n, []):
                    if layer in [LAYERS[L] for L in pad_layers(p)] and _pad_distance(x, y, p) < 1e-6:
                        return True
                for v in vias.get(n, []):
                    if math.hypot(v['x'] - x, v['y'] - y) < v['d'] / 2 - 1e-3:
                        return True
                for tj in by_net[n]:
                    if tj == ti:
                        continue
                    u = self.tracks[tj]
                    if u is not None and u['layer'] == layer and any(_seg_dist(x, y, a, b) < u['width'] / 2 - 1e-3 for a, b in
                                                   zip(u['pts'], u['pts'][1:])):
                        return True
                return False
            for ti in range(len(self.tracks)):
                t = self.tracks[ti]
                if t is None:
                    continue
                pts, w, L = t['pts'], t['width'], t['layer']
                for rev in (False, True):
                    seq = pts[::-1] if rev else pts
                    if touched(ti, seq[0][0], seq[0][1], L, w):
                        continue
                    # walk inwards in 0.05 mm steps to the first point touched by other copper
                    cut = None
                    for k, (a, b) in enumerate(zip(seq, seq[1:])):
                        l = math.hypot(b[0] - a[0], b[1] - a[1])
                        steps = max(int(l / 0.05), 1)
                        for i in range(1, steps + 1):
                            x = a[0] + (b[0] - a[0]) * i / steps
                            y = a[1] + (b[1] - a[1]) * i / steps
                            if touched(ti, x, y, L, w):
                                cut = (k, (x, y))
                                break
                        if cut:
                            break
                    if cut is None:
                        removed += sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(pts, pts[1:]))
                        self.tracks[ti] = None
                        changed = True
                        break
                    k, pt = cut
                    new = [pt] + seq[k + 1:]
                    removed += sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(seq[:k + 1] + [pt], (seq[:k + 1] + [pt])[1:]))
                    if len(new) < 2 or math.hypot(new[-1][0] - new[0][0], new[-1][1] - new[0][1]) < 1e-6 and len(new) == 2:
                        self.tracks[ti] = None
                    else:
                        pts = new[::-1] if rev else new
                        self.tracks[ti] = dict(t, pts=[[round(x, 6), round(y, 6)] for x, y in pts])
                        t = self.tracks[ti]
                    changed = True
                    if self.tracks[ti] is None:
                        break
            self.tracks = [t for t in self.tracks if t is not None]
        return removed

    def set_phantoms(self, rects):
        """Temporary F.Cu (and so via) obstacles, e.g. other blocks' parts mapped into a template cell; [] clears."""
        self.spec_clear()
        self.phantom[:] = False
        self._phantom_ver = getattr(self, '_phantom_ver', 0) + 1
        self._dirty[:] = True
        for x0, y0, x1, y1 in rects:
            win, m = self._mask_rect(x0, y0, x1, y1)
            self.phantom[win[0]:win[2] + 1, win[1]:win[3] + 1] |= m

    # ------------------------------------------------------------------ checks and output
    def _memo(self, kind, key, fn):
        """fn() memoised in the route cache under a content key (pure functions of the board only)."""
        if _CACHE is None:
            return fn()
        import copy
        k = kind + _GR.hash_arrays([np.frombuffer(repr(key).encode(), dtype=np.uint8)]) + _SALT
        hit = _CACHE.get(k)
        if hit is None:
            hit = fn()
            _CACHE.put(k, copy.deepcopy(hit))     # the caller goes on to modify what fn() returned
            return hit
        return copy.deepcopy(hit)

    def _placement(self):
        return [(r, p.x, p.y, p.rot) for r, p in self.parts.items() if p.placed]

    def check(self, since=(0, 0)):
        """Clearance check (see _check); memoised on the new copper and the copper raster."""
        return self._memo('check', (self.tracks[since[0]:], self.vias[since[1]:],
                                    _GR.hash_arrays([self.core]) if _CACHE is not None else None),
                          lambda: self._check(since))

    def _check(self, since=(0, 0)):
        """Clearance check of tracks and vias (from mark `since`) against other nets' copper, using the class
        clearance of either net (the larger wins). Rasterised, so it can miss by a few hundredths of a mm; KiCad
        DRC is the final word. Returns [(description, other net)]."""
        out = []

        def test(net, L, win, m):
            i0, j0, i1, j1 = win
            c = self.core[L, i0:i1 + 1, j0:j1 + 1]
            hit = m & (c > 0) & (c != self.net_id[net])
            return sorted({self.nets[v - 1] for v in np.unique(c[hit])})

        def near(net, L, geom, r_own):
            r_own -= 0.02                  # raster tolerance: gaps of exactly the rule value must pass
            found = []
            for extra in sorted({0.0, 0.05}):
                win, m = geom(r_own + extra)
                for o in test(net, L, win, m):
                    if extra == 0.0 or self.clearance(o) > self.clearance(net) + 1e-9:
                        found.append(o)
            return sorted(set(found))
        for t in self.tracks[since[0]:]:
            L = LAYERS.index(t['layer'])
            for a, b in zip(t['pts'], t['pts'][1:]):
                geom = lambda r: self._mask_seg(a[0], a[1], b[0], b[1], r)
                for o in near(t['net'], L, geom, t['width'] / 2 + self.clearance(t['net'])):
                    out.append(('%s %s %s-%s' % (short(t['net']), t['layer'], a, b), o))
        for v in self.vias[since[1]:]:
            geom = lambda r: self._mask_seg(v['x'], v['y'], v['x'], v['y'], r)
            for L in range(NL):
                for o in near(v['net'], L, geom, v['d'] / 2 + self.clearance(v['net'])):
                    out.append(('via %s (%s, %s) %s' % (short(v['net']), v['x'], v['y'], LAYERS[L]), o))
        return out

    def unrouted(self, refs=None):
        """Nets whose placed pads are not in one copper island (plane nets: pads without a via or hole)."""
        out = []
        byn = {}
        for r, p in self.parts.items():
            if not p.placed or (refs and r not in refs):
                continue
            for pd in p.pads():
                if pd['net'] and not pd['net'].startswith('unconnected-'):
                    byn.setdefault(pd['net'], []).append(pd)
        full = (0, 0, self.ny - 1, self.nx - 1)
        for net, pads in sorted(byn.items()):
            nid = self.net_id[net]
            if net in PLANES:
                for pd in pads:
                    if pd['layers'] == 'all':
                        continue
                    if not self.plane_reached(pd):
                        out.append((net, pd['ref'], pd['num']))
                continue
            if len(pads) < 2:
                continue
            labels = self.pad_components(net, pads)
            for pd, label in zip(pads[1:], labels[1:]):
                if label != labels[0]:
                    out.append((net, pd['ref'], pd['num']))
        return out

    def save(self, path, outline):
        data = dict(outline=outline, parts={r: dict(x=round(p.x, 4), y=round(p.y, 4), rot=p.rot) for r, p in self.parts.items() if p.placed},
                    tracks=self.tracks, vias=self.vias, zones=self.zones, keepouts=self.keepouts,
                    keepouts_rule=self.rule_keepouts, texts=self.texts, rects=self.rects, labels=self.labels)
        json.dump(data, open(path, 'w'), indent=0)
        if _CACHE is not None:
            _CACHE.save()
            print('route cache: %d hits, %d misses' % (_CACHE.hits, _CACHE.misses))


configure()
