"""Escape planning for dense (fine-pitch) packages, done before general routing.

Routing pin by pin lets early nets drop vias into the channels later pins need. A human
escapes a dense part as a whole:

1. Local buses: nets that join a dense package to nearby pads on the same copper side are
   routed first, on that side only (e.g. a level shifter's row straight into its connector).
   Only nets whose straight connections cross no other bus net are routed this way; a net
   that would cross (pin order not preserved) changes layer instead.
2. Dog-bones: pins that must change layer (plane nets, crossing bus nets, or nets with no
   other pad on the pin's side) get a stub straight out of the pad end and a via. The sites
   of one package are chosen together: candidate sites are tested against the board's own
   raster rules, mutual conflicts are geometric (with the raster margin), and a beam search
   keeps the assignment that places the most vias with the shortest stubs. Final native
   validation is still required.
"""
import math
from . import board as geometry

GAPS = (.1, .3, .55, .8, 1.1)                      # pad end to via edge (mm)
LATERAL = (0., .5, -.5, 1., -1., 1.5, -1.5, 2., -2., 2.5, -2.5)   # sideways via offset, in pin pitches
BEAM = 24


def dense_packages(bd, fine_pitch=.65, min_pads=6):
    """{original_ref: (smd pads, pitch)} for packages whose nearest pad spacing is at most fine_pitch."""
    groups = {}
    for part in bd.parts.values():
        if part.placed:
            for q in part.pads():
                if q['layers'] != 'all':
                    groups.setdefault(q['original_ref'], []).append(q)
    out = {}
    for ref, ps in groups.items():
        if len(ps) < min_pads:
            continue
        pitch = min(min(math.hypot(a['x'] - b['x'], a['y'] - b['y']) for b in ps if b is not a) for a in ps)
        if 0 < pitch <= fine_pitch:
            out[ref] = (ps, pitch)
    return out


def _clusters(values, tol=.3):
    """Centres of groups of values closer than tol (pin rows / columns)."""
    out = []
    for v in sorted(values):
        if out and v - out[-1][-1] < tol:
            out[-1].append(v)
        else:
            out.append([v])
    return [sum(c) / len(c) for c in out]


def _side(p):
    return geometry.pad_layers(p)[0]


def _key(p):
    return (p['ref'], p['num'])


def _outward(p, cx, cy):
    """(unit direction away from the package centre along the pad's axis, pad half-length along it, short side,
    whether the pad lies on the package edge rather than in its centre)."""
    vx, vy = p['x'] - cx, p['y'] - cy
    if p['w'] > p['h'] * 1.2 or (p['h'] <= p['w'] * 1.2 and abs(vx) >= abs(vy)):
        return (math.copysign(1., vx) if vx else 1., 0.), p['w'] / 2, p['h'], abs(vx) >= p['w'] / 2
    return (0., math.copysign(1., vy) if vy else 1.), p['h'] / 2, p['w'], abs(vy) >= p['h'] / 2


def _seg_seg(a, b, c, d):
    def orient(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])
    if orient(a, b, c) * orient(a, b, d) < 0 and orient(c, d, a) * orient(c, d, b) < 0:
        return 0.
    return min(geometry._seg_dist(*a, c, d), geometry._seg_dist(*b, c, d),
               geometry._seg_dist(*c, a, b), geometry._seg_dist(*d, a, b))


def _mst_edges(ps):
    """Prim's tree over pad centres: the straight connections a net would make."""
    pts = [(q['x'], q['y']) for q in ps]
    seen, edges = [0], []
    while len(seen) < len(pts):
        a, b = min(((i, j) for i in seen for j in range(len(pts)) if j not in seen),
                   key=lambda e: math.dist(pts[e[0]], pts[e[1]]))
        seen.append(b)
        edges.append((pts[a], pts[b]))
    return edges


def local_bus_nets(bd, packages, local_mm=6.):
    """(planar, crossing): single-side nets touching a dense package within local_mm. Planar nets' straight
    connections cross no other candidate's; crossing nets would have to change layer."""
    dense = {_key(q) for ps, _ in packages.values() for q in ps}
    cands = {}
    for net in bd.nets:
        if net in bd.config['planes'] or net.startswith('unconnected-'):
            continue
        ps = bd.pads_of(net)
        if len(ps) < 2 or any(q['layers'] == 'all' for q in ps) or len({_side(q) for q in ps}) > 1:
            continue
        if not any(_key(q) in dense for q in ps):
            continue
        span = math.hypot(max(q['x'] for q in ps) - min(q['x'] for q in ps),
                          max(q['y'] for q in ps) - min(q['y'] for q in ps))
        if span <= local_mm:
            cands[net] = (span, _side(ps[0]), _mst_edges(ps))
    crossing = set()
    nets = sorted(cands)
    for i, a in enumerate(nets):
        for b in nets[i + 1:]:
            if cands[a][1] == cands[b][1] and any(_seg_seg(*e, *f) == 0. for e in cands[a][2] for f in cands[b][2]):
                crossing.update((a, b))
    planar = sorted((n for n in cands if n not in crossing), key=lambda n: (cands[n][0], n))
    return planar, sorted(crossing)


def route_local_buses(bd, connect, packages, local_mm=6., time_left=lambda: True):
    """Route the planar local bus nets on their own side, shortest first. Returns the nets completed."""
    planar, crossing = local_bus_nets(bd, packages, local_mm)
    done = []
    for net in planar:
        if not time_left():
            break
        if connect(net, [geometry.LAYERS[_side(bd.pads_of(net)[0])]]):
            done.append(net)
    bd.__dict__['_bus_crossing'] = set(crossing)
    return done


def _needs_via(bd, p):
    net = p['net']
    if not net or net.startswith('unconnected-'):
        return False
    if net in bd.config['planes']:
        return not bd.plane_reached(p)
    ps = bd.pads_of(net)
    if len(ps) < 2 or bd._net_complete(net):
        return False
    if net in bd.__dict__.get('_bus_crossing', ()):
        return True
    others = [q for q in ps if _key(q) != _key(p)]
    return bool(others) and all(q['layers'] != 'all' and _side(q) != _side(p) for q in others)


class _Site:
    __slots__ = ('pad', 'via', 'pts', 'width', 'vd', 'drill', 'clear', 'length', 'layer')

    def __init__(self, pad, pts, width, vd, drill, clear, layer):
        self.pad, self.pts, self.width, self.vd, self.drill, self.clear = pad, pts, width, vd, drill, clear
        self.layer = layer
        self.via = pts[-1]
        self.length = sum(math.dist(a, b) for a, b in zip(pts, pts[1:]))

    def segs(self):
        return list(zip(self.pts, self.pts[1:]))

    def conflicts(self, o, rules, slack):
        """Copper clearance (class, and per-layer minimums wherever both objects exist), drill to copper and
        drill to drill, as the native rules apply them; slack covers the raster margin."""
        cls = max(self.clear, o.clear)
        lc = rules['layer_clearances']
        via_via = max([cls] + list(lc.values()))
        d = math.dist(self.via, o.via)
        if d < max((self.vd + o.vd) / 2 + via_via, self.drill / 2 + o.vd / 2 + rules['hole_clear'],
                   o.drill / 2 + self.vd / 2 + rules['hole_clear'],
                   (self.drill + o.drill) / 2 + rules['hole_to_hole']) + slack:
            return True
        for a, b in ((self, o), (o, self)):      # a's via against b's stub (on b's layer)
            need = max(a.vd / 2 + b.width / 2 + max(cls, lc.get(b.layer, 0.)),
                       a.drill / 2 + b.width / 2 + rules['hole_clear']) + slack
            if any(geometry._seg_dist(*a.via, *s) < need for s in b.segs()):
                return True
        if self.layer != o.layer:
            return False
        need = (self.width + o.width) / 2 + max(cls, lc.get(self.layer, 0.)) + slack
        return any(_seg_seg(*s, *t) < need for s in self.segs() for t in o.segs())


def _sites(bd, p, pitch, cx, cy, inward, masks):
    """Legal dog-bone sites of pad p against the current board: a stub straight out of the middle of the pad
    end, then an optional angled leg to a via in one of several staggered rows."""
    net = p['net']
    G = bd.pitch
    L = _side(p)
    d, half, short, edge = _outward(p, cx, cy)
    if not edge:
        return []                       # a central (exposed) pad: the plane drop handles it in the pad
    vd, drill = bd.via_size(net)
    width = min(p.get('escape_width') or bd.width(net, geometry.LAYERS[L]), short)
    key = (net, width)
    if key not in masks:
        win = masks['win']
        masks[key] = (bd.via_ok(net, win), bd.track_ok(net, L, width, win))
    vok, tok = masks[key]
    i0, j0 = masks['win'][0], masks['win'][1]

    def free(mask, x, y):
        i, j = int(round(y / G)) - i0, int(round(x / G)) - j0
        return 0 <= i < mask.shape[0] and 0 <= j < mask.shape[1] and bool(mask[i, j])

    def clear_line(a, b):
        n = max(int(math.dist(a, b) / (G / 2)), 1)
        return all(free(tok, a[0] + (b[0] - a[0]) * t / n, a[1] + (b[1] - a[1]) * t / n) for t in range(n + 1))

    out = []
    start = (p['x'], p['y'])
    for sign in ((1., -1.) if inward else (1.,)):
        ux, uy = d[0] * sign, d[1] * sign
        for gap in GAPS:
            along = half + gap + vd / 2
            for lat in LATERAL:
                x = round((p['x'] + ux * along - uy * lat * pitch) / G) * G
                y = round((p['y'] + uy * along + ux * lat * pitch) / G) * G
                if not free(vok, x, y):
                    continue
                if lat == 0.:
                    # straight on the pad axis: keep the via on the pad's centre line, off-grid if need be
                    x, y = (p['x'], y) if ux == 0 else (x, p['y'])
                    pts = [start, (x, y)]
                else:
                    # straight out of the pad end first, so the angled leg cannot clip the neighbouring pads
                    k = half + min(gap, .15)
                    knee = (p['x'] + ux * k, p['y'] + uy * k)
                    pts = [start, knee, (x, y)]
                if all(clear_line(a, b) for a, b in zip(pts, pts[1:])):
                    out.append(_Site(p, pts, width, vd, drill, bd.clearance(net), geometry.LAYERS[L]))
    return out


def _assign(pins, options, rules, slack):
    """Beam search: at most one site per pin, no two chosen sites in conflict; most vias, then shortest stubs."""
    order = sorted((p for p in pins if options[_key(p)]), key=lambda p: (len(options[_key(p)]), _key(p)))
    beam = [((), 0., ())]                       # (chosen sites, total length, skipped pins)
    for p in order:
        nxt = []
        for chosen, length, skipped in beam:
            for s in options[_key(p)]:
                if not any(s.conflicts(c, rules, slack) for c in chosen):
                    nxt.append((chosen + (s,), length + s.length, skipped))
            nxt.append((chosen, length, skipped + (p,)))
        nxt.sort(key=lambda b: (-len(b[0]), b[1]))
        beam = nxt[:BEAM]
    best = beam[0]
    return list(best[0]), list(best[2]) + [p for p in pins if not options[_key(p)]]


def plan_dogbones(bd, packages, time_left=lambda: True, inward_gap=1.2):
    """Commit stub + via for every pin of the dense packages that must change layer. Returns
    (number placed, [(ref, pin, net) without a site]). Two-row packages with at least `inward_gap` mm
    between the rows also get sites between the rows."""
    placed = 0
    unplaced = []
    rules = dict(hole_to_hole=geometry.HOLE_TO_HOLE, hole_clear=geometry.HOLE_CLEAR,
                 layer_clearances=dict(bd.config['layer_clearances']))
    slack = geometry.MARGIN + bd.pitch / 2
    for ref, (ps, pitch) in sorted(packages.items()):
        if not time_left():
            break
        pins = [p for p in ps if _needs_via(bd, p)]
        if not pins:
            continue
        cx, cy = sum(q['x'] for q in ps) / len(ps), sum(q['y'] for q in ps) / len(ps)
        xs, ys = [q['x'] for q in ps], [q['y'] for q in ps]
        inward = any(len(c) == 2 and c[1] - c[0] >= inward_gap for c in (_clusters(ys), _clusters(xs)))
        r = 3.5
        masks = {'win': bd._win(min(xs) - r, min(ys) - r, max(xs) + r, max(ys) + r)}
        options = {_key(p): _sites(bd, p, pitch, cx, cy, inward, masks) for p in pins}
        chosen, skipped = _assign(pins, options, rules, slack)
        for s in chosen:
            p = s.pad
            bd.add_track(p['net'], geometry.LAYERS[_side(p)], s.pts, s.width)
            bd.add_via(p['net'], *s.via)
            placed += 1
        unplaced += [(ref, p['original_pin'], p['net']) for p in skipped]
    return placed, unplaced
