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
              accel=None, spec=None, cache=None, cache_name='', salt=b''):
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
    global MARGIN, _net_class, _GR, _VERIFY, _HYBRID, REPRODUCE, _BUDGET, _FIELD, _CACHE, _SALT, _SPEC_N, _SPEC_POOL
    LAYERS = list(layers)
    NL = len(LAYERS)
    G = pitch
    CLASSES = dict(classes or {'Default': (0.15, 0.15, 0.5, 0.3)})
    CLEAR = CLASSES['Default'][1]
    VIA_D, VIA_DRILL = CLASSES['Default'][2], CLASSES['Default'][3]
    PLANES = dict(planes or {})
    EDGE_CLEAR, HOLE_CLEAR = edge_clear, hole_clear
    LAYER_WIDTHS = {k: dict(v) for k, v in (layer_widths or {}).items()}
    PAD_GROW, MARGIN = pad_grow, margin
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
                       _BUDGET, _FIELD, _GR.__file__, os.path.getmtime(_GR._LIB._name))).encode())
        path = where or os.path.join(os.path.expanduser('~/Library/Caches/gridroute'),
                                     hashlib.sha1(cache_name.encode()).hexdigest()[:12] + '-routes.pkl')
        _CACHE, _SALT = _GR.SearchCache(path), h.hexdigest()
    if _SPEC_POOL is not None:
        _SPEC_POOL.shutdown(wait=True)
    _SPEC_N = int(spec if spec is not None else _env('SPEC', '4')) if _GR is not None else 0
    _SPEC_POOL = ThreadPoolExecutor(_SPEC_N) if _SPEC_N > 0 else None


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
            r = (p['rot'] + self.rot) % 180
            assert min(abs(r), abs(r - 90), abs(r - 180)) < 1, (self.ref, p['num'], r)
            w, h = (p['h'], p['w']) if abs(r - 90) < 1 else (p['w'], p['h'])
            out.append(dict(num=p['num'], x=self.x + dx, y=self.y + dy, w=w, h=h, shape=p['shape'],
                            layers=p['layers'], npth=p['npth'], drill=p['drill'], ref=self.ref, rr=p.get('rr', 0.0),
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


def dedup(pts, eps=1e-6):
    out = [pts[0]]
    for p in pts[1:]:
        if math.hypot(p[0] - out[-1][0], p[1] - out[-1][1]) > eps:
            out.append(p)
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
        self.nx, self.ny = int(math.ceil(w / G)) + 1, int(math.ceil(h / G)) + 1
        self.occ = np.zeros((NL, self.ny, self.nx), dtype=np.int16)
        self.core = np.zeros((NL, self.ny, self.nx), dtype=np.int16)
        self.thru = np.zeros((self.ny, self.nx), dtype=np.int16)
        self.smd = np.zeros((self.ny, self.nx), dtype=np.int16)      # SMD pad copper (either side): no vias
        self.phantom = np.zeros((self.ny, self.nx), dtype=bool)      # temporary F.Cu / via obstacles (templates)
        self.comps, self.pin_net = comps, pin_net
        self.fpdefs = footprints
        self.nets = sorted(set(self.pin_net.values()))
        self.net_id = {n: i + 1 for i, n in enumerate(self.nets)}
        self.cls = {n: net_class(n) for n in self.nets}
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
        e = int(round((EDGE_CLEAR - CLEAR) / G)) + 1
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

    # ------------------------------------------------------------------ placement
    def place(self, ref, x, y, rot=0.0, gap=0.15):
        """Place a footprint; courtyards closer than `gap` to an already placed part are recorded in self.overlaps."""
        p = self.parts[ref]
        assert not p.placed, ref
        rot += self.turn.get(ref, 0)
        p.x, p.y, p.rot, p.placed = x, y, rot % 360, True
        self._pad_index = None
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
        if pad['npth']:
            win, m = self._mask_seg(pad['x'], pad['y'], pad['x'], pad['y'], pad['drill'] / 2 + HOLE_CLEAR - CLEAR + PAD_GROW)
            for L in range(NL):
                self._paint(self.occ, L, win, m, KEEP)
            return
        Ls = pad_layers(pad)
        win, m = self._pad_mask(pad, PAD_GROW + self.extra(pad['net']))
        cwin, cm = self._pad_mask(pad, 0.0)
        for L in Ls:
            self._paint(self.occ, L, win, m, nid)
            self._paint(self.core, L, cwin, cm, nid)
        if pad['layers'] == 'all':
            self._paint(self.thru, None, cwin, cm, nid)
        else:
            self._paint(self.smd, None, cwin, cm, 1)

    # ------------------------------------------------------------------ copper
    def add_track(self, net, layer, pts, width=None):
        width = width or self.width(net, layer)
        t = dict(net=net, layer=layer, width=round(width, 4), pts=[[round(x, 4), round(y, 4)] for x, y in pts])
        self._op('track', (net, LAYERS.index(layer), [tuple(q) for q in pts], width), t)
        self.tracks.append(t)

    def _do_track(self, net, L, pts, width):
        nid = self.net_id[net]
        for a, b in zip(pts, pts[1:]):
            self._paint_seg(self.occ, L, a, b, width / 2 + self.extra(net), nid)
            self._paint_seg(self.core, L, a, b, width / 2, nid)

    def add_via(self, net, x, y, d=None, drill=None):
        d0, dr0 = self.via_size(net)
        d, drill = d or d0, drill or dr0
        v = dict(net=net, x=round(x, 4), y=round(y, 4), d=d, drill=drill)
        self._op('via', (net, x, y, d), v)
        self.vias.append(v)

    def _do_via(self, net, x, y, d):
        nid = self.net_id[net]
        for L in range(NL):
            self._paint_seg(self.occ, L, (x, y), (x, y), d / 2 + self.extra(net), nid)
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
        return _GR.dilate(self.occ, Ls, win, _GR.disc_span(r / G), excl=nid if isinstance(nid, tuple) else (nid,),
                          extra=self.phantom if ph else None, extra_on=[1] + [0] * (NL - 1) if ph else None,
                          reduce_or=reduce_or)

    def _blocked_np(self, L, nid, win, r):
        i0, j0, i1, j1 = win
        rc = int(math.ceil(r / G))
        pi0, pj0 = max(i0 - rc, 0), max(j0 - rc, 0)
        pi1, pj1 = min(i1 + rc, self.ny - 1), min(j1 + rc, self.nx - 1)
        src = self.occ[L, pi0:pi1 + 1, pj0:pj1 + 1]
        other = src != 0
        for n in (nid if isinstance(nid, tuple) else (nid,)):
            other &= src != n
        if L == 0 and self.phantom.any():
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

    def _near_smd_np(self, win, r):
        i0, j0, i1, j1 = win
        rc = int(math.ceil(r / G))
        pi0, pj0 = max(i0 - rc, 0), max(j0 - rc, 0)
        pi1, pj1 = min(i1 + rc, self.ny - 1), min(j1 + rc, self.nx - 1)
        src = self.smd[pi0:pi1 + 1, pj0:pj1 + 1] != 0
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

    def _reach(self, net, half):
        return half + CLEAR + self.extra(net) + MARGIN

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

    def route(self, net, a, b, **kw):
        """Connect the copper island of pad a=(ref, num) to the island of pad b."""
        pa, pb = self.parts[a[0]].pad(a[1]), self.parts[b[0]].pad(b[1])
        assert pa['net'] == net and pb['net'] == net, (net, a, b, pa['net'], pb['net'])
        return self._route(net, [pa], [pb], **kw)

    def connect(self, net, layers=None, only=None, **kw):
        """Route every placed pad of the net into one tree (nearest pad first). only: restrict to these refs."""
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
                if not self._route(net, [p], tree, layers=None, margin=kw.get('margin', 3.0) * 3,
                                   **{k: v for k, v in kw.items() if k != 'margin'}):
                    self.failed.append((net, p['ref'], p['num']))
                    ok = False
                    continue
            tree.append(p)
            if id(p) in comp:
                joined.add(comp[id(p)])     # routed: its whole component now belongs to the tree
        return ok

    def _content_key(self, win, *args):
        """Content hash of the board arrays over the tiles covering win, plus args (route cache key)."""
        ti0, tj0, ti1, tj1 = win[0] // TILE, win[1] // TILE, win[2] // TILE, win[3] // TILE
        d = self._dirty[ti0:ti1 + 1, tj0:tj1 + 1]
        if d.any():
            tl = np.argwhere(d) + (ti0, tj0)
            self._th[tl[:, 0], tl[:, 1]] = _GR.tile_hashes([self.occ, self.core, self.thru, self.smd, self.phantom], TILE, tl)
            d[:] = False
        return _GR.hash_arrays([self._th[ti0:ti1 + 1, tj0:tj1 + 1]]) + repr((_SALT, win, args))

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
        if res is None:
            res = self._route_search(net, src_pads, dst_pads, *args)
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
        if not explicit and self.cls[net] in LAYER_WIDTHS:   # search at the widest width the emitted layers may use
            width = max(self.width(net, l) for l in (layers or LAYERS))
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
        if _GR is not None:
            blk[Ls] = self._blocked_gr(Ls, nid, win, r)
            vok = ~self._blocked_gr(list(range(NL)), nid, win, vr, reduce_or=True)
        if _GR is None or _VERIFY:
            if _VERIFY:
                blk_gr, vok_gr = blk.copy(), vok.copy()
            for L in Ls:
                blk[L] = self._blocked(L, nid, win, r)
            vok = np.ones((H, W), dtype=bool)
            for L in range(NL):
                vok &= ~self._blocked(L, nid, win, vr)
            if _VERIFY:
                assert np.array_equal(blk, blk_gr) and np.array_equal(vok, vok_gr), ('_route masks', net, win)
        vok &= ~self._near_smd(win, vd / 2 + 0.1)     # no via-in-pad, not even on the net's own pads
        # the hole needs HOLE_CLEAR from other copper as well (approximated by the pad reach above)
        blk[:, 0, :] = blk[:, -1, :] = blk[:, :, 0] = blk[:, :, -1] = True
        # sources as flat state indices (native scan) or (L, i, j) rows (numpy fallback, verify mode)
        src = _GR.mask_scan(isl_s, blk) if _GR is not None and not _VERIFY else np.argwhere(isl_s & ~blk)
        tgt = isl_d & ~blk
        if not len(src) or not tgt.any():
            return ('fail',)
        lc = [1.0] * NL
        for l, m in (layer_cost or {}).items():
            lc[LAYERS.index(l)] = m
        path = self._astar(blk, vok, src, tgt, Ls, lc, via_cost / G, pref or {}, turn / G, weight, max_exp)
        if path is None:
            return ('fail',)
        return ('path', net, path, win, None if self.cls[net] in LAYER_WIDTHS and not explicit else width, (vd, vdr), H, W)
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

    def _astar_gr_run(self, blk, vok, src, tgt, Ls, lc, vcost, pref, turn, weight, max_exp, read_set=False):
        nL, H, W = blk.shape
        HW = H * W
        lok = np.zeros(NL, dtype=bool)
        lok[Ls] = True
        flat = src if src.ndim == 1 else src[:, 0] * HW + src[:, 1] * W + src[:, 2]
        s = _GR.astar(blk, vok, flat, tgt, self._move_costs(Ls, lc, pref), lok,
                      vcost, turn, weight * min(lc[L] for L in Ls), max_exp,
                      **(dict(budget=_BUDGET, weight=weight, read_set=read_set) if _HYBRID else {}))
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

    def _emit(self, net, path, win, width, via, H, W):
        i0, j0 = win[0], win[1]
        xy = lambda k: ((k[2] + j0) * G, (k[1] + i0) * G)
        seg = [path[0]]
        for k in path[1:]:
            if k[0] != seg[-1][0]:
                self._flush(net, seg, xy, width)
                self.add_via(net, *xy(k), d=via[0], drill=via[1])
                seg = [k]
            else:
                seg.append(k)
        self._flush(net, seg, xy, width)

    def _flush(self, net, seg, xy, width):
        if len(seg) < 2:
            return
        pts = [seg[0]]
        for a, b, c in zip(seg, seg[1:], seg[2:]):
            if (b[1] - a[1], b[2] - a[2]) != (c[1] - b[1], c[2] - b[2]):
                pts.append(b)
        pts.append(seg[-1])
        self.add_track(net, LAYERS[seg[0][0]], [xy(k) for k in pts], width)

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
            vok &= ~self._blocked(LL, nid, win, self._reach(net, vd / 2))
        if not in_pad:
            vok &= ~self._near_smd(win, vd / 2 + 0.1)
        sok = ~self._blocked(L, nid, win, self._reach(net, width / 2))
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
        # keep drilled holes apart (hole-to-hole 0.25 mm), whatever their net
        ys, xs = self._grid(win)
        for v in self.vias:
            if abs(v['x'] - p['x']) < max_r + 2 and abs(v['y'] - p['y']) < max_r + 2:
                vok &= (xs - v['x']) ** 2 + (ys - v['y']) ** 2 >= ((v['drill'] + vdr) / 2 + 0.3) ** 2
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
                for t in range(k + 1):
                    ii = int(round((p['y'] + dy * t * G) / G)) - i0
                    jj = int(round((p['x'] + dx * t * G) / G)) - j0
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
        self.spec_clear()
        self._dropped = set()
        nets = set(nets)
        gone = [k for k in range(since, len(self.ops)) if self.ops[k][3] is not None and self.ops[k][3]['net'] in nets]
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
        for a in (self.occ, self.core):
            a[:, i0:i1 + 1, j0:j1 + 1] = 0
        for a in (self.thru, self.smd):
            a[i0:i1 + 1, j0:j1 + 1] = 0
        self._dirty[i0 // TILE:i1 // TILE + 1, j0 // TILE:j1 // TILE + 1] = True
        # replay into a scratch board region: each operation repaints in full (outside the region that is idempotent)
        self._edges()
        for kind, args, w, item in self.ops:
            if w[0] <= i1 and i0 <= w[2] and w[1] <= j1 and j0 <= w[3]:
                getattr(self, '_do_' + kind)(*args)
        return len(gone)

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
        self.tracks, self.vias, n = self._memo('prune', (self.tracks, self.vias, self._placement()), run)
        return n

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
                    if abs(p['x'] - v['x']) < p['w'] / 2 + v['d'] / 2 and abs(p['y'] - v['y']) < p['h'] / 2 + v['d'] / 2:
                        layers |= {LAYERS[L] for L in pad_layers(p)}
                        at_pad = True
                if len(layers) >= 2:
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
                    if layer in [LAYERS[L] for L in pad_layers(p)] and abs(p['x'] - x) <= p['w'] / 2 + w / 2 and \
                            abs(p['y'] - y) <= p['h'] / 2 + w / 2:
                        return True
                for v in vias.get(n, []):
                    if math.hypot(v['x'] - x, v['y'] - y) <= v['d'] / 2 + w / 2:
                        return True
                for tj in by_net[n]:
                    if tj == ti:
                        continue
                    u = self.tracks[tj]
                    if u['layer'] == layer and any(_seg_dist(x, y, a, b) <= (u['width'] + w) / 2 for a, b in
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
                        self.tracks[ti] = dict(t, pts=[[round(x, 4), round(y, 4)] for x, y in pts])
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
                    w = self._win(pd['x'] - 3, pd['y'] - 3, pd['x'] + 3, pd['y'] + 3)
                    isl = self._island(nid, self._pad_seeds(pd, w), w)
                    thr = self.thru[w[0]:w[2] + 1, w[1]:w[3] + 1] == nid
                    if not (isl.any(axis=0) & thr).any():
                        out.append((net, pd['ref'], pd['num']))
                continue
            if len(pads) < 2:
                continue
            xs, ys = [p['x'] for p in pads], [p['y'] for p in pads]
            w = self._win(min(xs) - 15, min(ys) - 15, max(xs) + 15, max(ys) + 15)
            isl = self._island(nid, self._pad_seeds(pads[0], w), w)
            for pd in pads[1:]:
                if not any(isl[s] for s in self._pad_seeds(pd, w)):
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
