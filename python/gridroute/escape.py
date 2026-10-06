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
        # coincident pads (e.g. a USB-C connector's paired A/B pins) are one position
        gaps = [math.hypot(a['x'] - b['x'], a['y'] - b['y']) for a in ps for b in ps if b is not a]
        gaps = [g for g in gaps if g > 1e-6]
        pitch = min(gaps) if gaps else 0.
        if 0 < pitch <= fine_pitch:
            out[ref] = (ps, pitch)
    return out


def package_centre(bd, ref):
    """Centre of all of a package's pads, plated ones included: a single row of SMD pins (a connector) then
    still has an inside (towards its shell / body) and an outside."""
    ps = [q for part in bd.parts.values() if part.placed for q in part.pads() if q['original_ref'] == ref]
    xs = sorted(q['x'] for q in ps);ys = sorted(q['y'] for q in ps)
    return (xs[0] + xs[-1]) / 2, (ys[0] + ys[-1]) / 2


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
    """One planned escape of a pin: copper capsules (a, b, radius, layer) and optionally a via at the end.
    kind 'via': a dog-bone (stub + via); kind 'neck': a neck-down stub to where the class width fits,
    with a class-width landing capsule at its end."""
    __slots__ = ('pad', 'kind', 'capsules', 'via', 'vd', 'drill', 'clear', 'length', 'stub', 'width')

    def __init__(self, pad, kind, pts, width, clear, layer, via=None, vd=0., drill=0., landing=0.):
        self.pad, self.kind, self.clear, self.width, self.stub = pad, kind, clear, width, pts
        self.capsules = [(a, b, width / 2, layer) for a, b in zip(pts, pts[1:])]
        if landing:
            self.capsules.append((pts[-1], pts[-1], landing / 2, layer))
        self.via, self.vd, self.drill = via, vd, drill
        self.length = sum(math.dist(a, b) for a, b in zip(pts, pts[1:]))

    def conflicts(self, o, rules, slack):
        """Copper clearance (class, and per-layer minimums wherever both objects exist), drill to copper and
        drill to drill, as the native rules apply them; slack covers the raster margin."""
        cls = max(self.clear, o.clear)
        lc = rules['layer_clearances']
        for a, b, r, L in self.capsules:
            for c, d, q, M in o.capsules:
                if L == M and _seg_seg(a, b, c, d) < r + q + max(cls, lc.get(L, 0.)) + slack:
                    return True
        for x, y in ((self, o), (o, self)):
            if x.via is None:
                continue
            for c, d, q, M in y.capsules:       # x's via against y's copper (on y's layer)
                need = max(x.vd / 2 + q + max(cls, lc.get(M, 0.)), x.drill / 2 + q + rules['hole_clear']) + slack
                if geometry._seg_dist(*x.via, c, d) < need:
                    return True
        if self.via is not None and o.via is not None:
            via_via = max([cls] + list(lc.values()))
            held = lambda v: max(v.vd / 2, v.drill / 2 + rules['hole_clear'] - v.clear)   # reserved disk radius
            if math.dist(self.via, o.via) < max(held(self) + o.vd / 2 + via_via, held(o) + self.vd / 2 + via_via,
                                                self.drill / 2 + o.vd / 2 + rules['hole_clear'],
                                                o.drill / 2 + self.vd / 2 + rules['hole_clear'],
                                                (self.drill + o.drill) / 2 + rules['hole_to_hole']) + slack:
                return True
        return False


def _masks(bd, masks, net, L, width, via=False):
    key = (net, L, width, via)
    if key not in masks:
        masks[key] = bd.via_ok(net, masks['win']) if via else bd.track_ok(net, L, width, masks['win'])
    return masks[key]


def _free(bd, masks, mask, x, y):
    G = bd.pitch
    i, j = int(round(y / G)) - masks['win'][0], int(round(x / G)) - masks['win'][1]
    return 0 <= i < mask.shape[0] and 0 <= j < mask.shape[1] and bool(mask[i, j])


def _clear_line(bd, masks, mask, a, b):
    n = max(int(math.dist(a, b) / (bd.pitch / 2)), 1)
    return all(_free(bd, masks, mask, a[0] + (b[0] - a[0]) * t / n, a[1] + (b[1] - a[1]) * t / n) for t in range(n + 1))


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
    vok = _masks(bd, masks, net, L, width, via=True)
    tok = _masks(bd, masks, net, L, width)
    out = []
    start = (p['x'], p['y'])
    for sign in ((1., -1.) if inward else (1.,)):
        ux, uy = d[0] * sign, d[1] * sign
        for gap in GAPS:
            along = half + gap + vd / 2
            for lat in LATERAL:
                x = round((p['x'] + ux * along - uy * lat * pitch) / G) * G
                y = round((p['y'] + uy * along + ux * lat * pitch) / G) * G
                if not _free(bd, masks, vok, x, y):
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
                if all(_clear_line(bd, masks, tok, a, b) for a, b in zip(pts, pts[1:])):
                    out.append(_Site(p, 'via', pts, width, bd.clearance(net), geometry.LAYERS[L],
                                     via=(x, y), vd=vd, drill=drill))
    return out


def _needs_neck(bd, p):
    """A pad whose net may only leave it at a neck-down width (Board.add_neck rules) and has no stub yet."""
    net = p['net']
    if not p.get('escape_width') or not bd.__dict__.get('necks') or not net or net in bd.config['planes']:
        return False
    if len(bd.pads_of(net)) < 2 or bd._net_complete(net):
        return False
    return p['escape_width'] < bd.width(net, geometry.LAYERS[_side(p)])


def _neck_sites(bd, p, cx, cy, masks, reach=1.5):
    """Straight neck-down stubs from the pad centre (outwards, then 45 degrees either side) to the first
    grid point where the class width fits, each with a class-width landing for the continuation."""
    net = p['net']
    G = bd.pitch
    L = _side(p)
    layer = geometry.LAYERS[L]
    ew, full = p['escape_width'], bd.width(net, layer)
    thin = _masks(bd, masks, net, L, ew)
    fat = _masks(bd, masks, net, L, full)
    d, half, short, edge = _outward(p, cx, cy)
    c, s_ = math.cos(math.pi / 4), math.sin(math.pi / 4)
    dirs = [d, (d[0] * c - d[1] * s_, d[0] * s_ + d[1] * c), (d[0] * c + d[1] * s_, -d[0] * s_ + d[1] * c)]
    out = []
    for ux, uy in dirs:
        for k in range(1, int((reach + half) / (G / 2)) + 1):
            x, y = p['x'] + ux * k * G / 2, p['y'] + uy * k * G / 2
            if not _free(bd, masks, thin, x, y):
                break
            ex, ey = round(x / G) * G, round(y / G) * G
            if not _free(bd, masks, fat, ex, ey):
                continue
            if _clear_line(bd, masks, thin, (p['x'], p['y']), (ex, ey)):
                out.append(_Site(p, 'neck', [(p['x'], p['y']), (ex, ey)], ew, bd.clearance(net), layer, landing=full))
            break
    return out


def _assign(pins, options, rules, slack):
    """Beam search: at most one site per pin, no two chosen sites in conflict; most escapes, then shortest."""
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


def plan_dogbones(bd, packages, time_left=lambda: True, inward_gap=1.2, commit=True):
    """Plan the escapes of every dense package's pins together: dog-bones (stub + via) for pins that must
    change layer and, on boards with neck-down rules, the neck-down stubs of power pins. Returns
    (number placed, [(ref, pin, net) without an escape]). Two-row packages with at least `inward_gap` mm
    between the rows also get via sites between the rows.

    Neck-down stubs and plane pins' dog-bones are copper at once (they are final). Signal dog-bones are
    copper when commit=True; with commit=False their stub and a via-sized disk on every layer are reserved
    instead (Board.reserve), so the net keeps its way out until it is routed, and later stages (power nets,
    plane drops, other signals) must route around them. A neck stub's landing is always reserved."""
    placed = 0
    unplaced = []
    rules = dict(hole_to_hole=geometry.HOLE_TO_HOLE, hole_clear=geometry.HOLE_CLEAR,
                 layer_clearances=dict(bd.config['layer_clearances']))
    slack = geometry.MARGIN + bd.pitch / 2
    for ref, (ps, pitch) in sorted(packages.items()):
        if not time_left():
            break
        necks = [p for p in ps if _needs_neck(bd, p)]
        vias = [p for p in ps if _needs_via(bd, p) and not any(_key(q) == _key(p) for q in necks)]
        pins = necks + vias
        if not pins:
            continue
        cx, cy = package_centre(bd, ref)
        xs, ys = [q['x'] for q in ps], [q['y'] for q in ps]
        inward = any(len(c) == 2 and c[1] - c[0] >= inward_gap for c in (_clusters(ys), _clusters(xs)))
        r = 3.5
        masks = {'win': bd._win(min(xs) - r, min(ys) - r, max(xs) + r, max(ys) + r)}
        options = {_key(p): _neck_sites(bd, p, cx, cy, masks) for p in necks}
        options.update({_key(p): _sites(bd, p, pitch, cx, cy, inward, masks) for p in vias})
        chosen, skipped = _assign(pins, options, rules, slack)
        for s in chosen:
            p = s.pad
            layer = geometry.LAYERS[_side(p)]
            if s.kind == 'neck':
                bd.add_track(p['net'], layer, s.stub, s.width)
                a, b, rad, L = s.capsules[-1]
                bd.reserve(p['net'], layer, a, b, rad)
            elif commit or p['net'] in bd.config['planes']:
                bd.add_track(p['net'], layer, s.stub, s.width)
                bd.add_via(p['net'], *s.via)
            else:
                for a, b, rad, L in s.capsules:
                    bd.reserve(p['net'], L, a, b, rad)
                # big enough that others also keep the drill-to-copper clearance of the future via
                rad = max(s.vd / 2, s.drill / 2 + geometry.HOLE_CLEAR - s.clear)
                bd.reserve(p['net'], None, s.via, s.via, rad)
            placed += 1
        unplaced += [(ref, p['original_pin'], p['net']) for p in skipped]
    return placed, unplaced


def _flood(free, seed):
    """Cells of `free` connected to `seed` with the A* move set (8 neighbours; a diagonal needs both orthogonal
    neighbours free)."""
    import numpy as np
    reach = seed & free
    while True:
        n = reach.copy()
        n[1:, :] |= reach[:-1, :]
        n[:-1, :] |= reach[1:, :]
        n[:, 1:] |= reach[:, :-1]
        n[:, :-1] |= reach[:, 1:]
        orth = n & free
        d = np.zeros_like(reach)
        d[1:, 1:] |= reach[:-1, :-1] & orth[:-1, 1:] & orth[1:, :-1]
        d[1:, :-1] |= reach[:-1, 1:] & orth[:-1, :-1] & orth[1:, 1:]
        d[:-1, 1:] |= reach[1:, :-1] & orth[1:, 1:] & orth[:-1, :-1]
        d[:-1, :-1] |= reach[1:, 1:] & orth[1:, :-1] & orth[:-1, 1:]
        n = orth | (d & free)
        if (n == reach).all():
            return reach
        reach = n


def audit_escapes(bd, refs=None, margin=.3, reach_mm=3.):
    """Check that every connected SMD pin can escape its package on the current board.

    For each pin of a routed net: 'exit' if a track of the net's width (its permitted escape width at a
    neck-down pad) reaches beyond the package outline + margin on the pin's layer within reach_mm;
    'via' if, where the net must change layer, such a track also reaches a legal via site. Then the
    dense packages' via-needing pins are assigned sites jointly (plan_dogbones, nothing committed) to
    find pins whose individual escapes exist but cannot all be used at once. Returns a report dict:
    {'pins': n, 'no_exit': [...], 'no_via': [...], 'joint_unplaced': [...], 'packages': {...}}."""
    import numpy as np
    G = bd.pitch
    groups = {}
    for part in bd.parts.values():
        if part.placed:
            for q in part.pads():
                groups.setdefault(q['original_ref'], []).append(q)
    report = dict(pins=0, no_exit=[], no_via=[], joint_unplaced=[], packages={})
    for ref, ps in sorted(groups.items()):
        if refs and ref not in refs:
            continue
        x0, x1 = min(q['x'] - q['w'] / 2 for q in ps) - margin, max(q['x'] + q['w'] / 2 for q in ps) + margin
        y0, y1 = min(q['y'] - q['h'] / 2 for q in ps) - margin, max(q['y'] + q['h'] / 2 for q in ps) + margin
        stats = dict(pins=0, exit=0, via_needed=0, via=0)
        for p in ps:
            net = p['net']
            if p['layers'] == 'all' or not net or net.startswith('unconnected-'):
                continue
            if len(bd.pads_of(net)) < 2 and net not in bd.config['planes']:
                continue
            if net in bd.config['planes'] and (bd.plane_reached(p) or not _outward(p, *package_centre(bd, ref))[3]):
                continue        # already on its plane, or a central (exposed) pad: a via in the pad serves it
            stats['pins'] += 1
            report['pins'] += 1
            L = _side(p)
            width = p.get('escape_width') or bd.width(net, geometry.LAYERS[L])
            win = bd._win(p['x'] - reach_mm, p['y'] - reach_mm, p['x'] + reach_mm, p['y'] + reach_mm)
            free = bd.track_ok(net, L, width, win)
            pw, pm = bd._pad_mask(p, 0.)
            seed = np.zeros_like(free)
            a0, b0, a1, b1 = max(pw[0], win[0]), max(pw[1], win[1]), min(pw[2], win[2]), min(pw[3], win[3])
            seed[a0 - win[0]:a1 - win[0] + 1, b0 - win[1]:b1 - win[1] + 1] = \
                pm[a0 - pw[0]:a1 - pw[0] + 1, b0 - pw[1]:b1 - pw[1] + 1]
            free |= seed
            reach = _flood(free, seed)
            ii, jj = np.nonzero(reach)
            xs, ys = (jj + win[1]) * G, (ii + win[0]) * G
            outside = ((xs < x0) | (xs > x1) | (ys < y0) | (ys > y1)).any()
            pin = (ref, p['original_pin'], net)
            if outside:
                stats['exit'] += 1
            else:
                report['no_exit'].append(pin)
            needs = (net in bd.config['planes'] and not bd.plane_reached(p)) or _needs_via(bd, p)
            if needs:
                stats['via_needed'] += 1
                if (reach & bd.via_ok(net, win)).any():
                    stats['via'] += 1
                else:
                    report['no_via'].append(pin)
        report['packages'][ref] = stats
    dense = dense_packages(bd)
    if refs:
        dense = {k: v for k, v in dense.items() if k in refs}
    saved = (list(bd.tracks), list(bd.vias), list(bd.ops))
    import copy
    trial = copy.copy(bd)          # plan on a shallow copy whose copper ops are discarded afterwards
    trial.__dict__ = dict(bd.__dict__)
    for k in ('occ', 'core', 'thru', 'smd', 'no_via', '_th', '_dirty'):
        trial.__dict__[k] = bd.__dict__[k].copy()
    trial.tracks, trial.vias, trial.ops = list(bd.tracks), list(bd.vias), list(bd.ops)
    _, unplaced = plan_dogbones(trial, dense, commit=False)
    report['joint_unplaced'] = unplaced
    assert (bd.tracks, bd.vias, bd.ops) == saved
    return report
