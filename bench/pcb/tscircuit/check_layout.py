#!/usr/bin/env python3
"""Router-independent check of a layout.json against its model.json: connectivity and clearance.

Used for both engines alike, so a board without a KiCad source (dut-pi) is judged the same way as the others;
where a KiCad board exists, its DRC (kicad_io.py export, kicad-cli, accept.py) remains the reference.

Copper is a core (point, segment or polygon) grown by a radius: pads as topor reads them (KiCad orientation),
track segments, vias on every layer, small copper zones of the model (pieces) as polygons. Planes: a via or
through pad of a plane net joins its plane when it is inside one of the net's pours on that layer (the whole
outline when the net has no regions); an SMD pad joins one on its own layer.

Reports per net whether all its pads are joined, the open pad groups, different-net clearance violations
(the larger of the two class clearances and pad local clearances; cold_routing.class_clearance between two nets
of one class; hole clearance for drills), copper too close
to the board edge, via drills too close to other drills, track width below the net's class width and signal
copper laid on plane layers. Overlaps (no gap at all: a short, or copper through an obstacle of no net) are also
counted and listed on their own. Pad-to-pad gaps of the placed board are not counted (they are not routing).
"""
import argparse
import collections
import json
import math
import pathlib

TOL = 0.002  # mm


def seg_dist(a, b, c, d):
    def pt(p, a, b):
        dx, dy = b[0] - a[0], b[1] - a[1]
        L = dx * dx + dy * dy
        t = 0.0 if L == 0 else max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L))
        return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)

    def cross(o, p, q):
        return (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0])
    d1, d2, d3, d4 = cross(a, b, c), cross(a, b, d), cross(c, d, a), cross(c, d, b)
    if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and d1 and d2 and d3 and d4:
        return 0.0
    return min(pt(c, a, b), pt(d, a, b), pt(a, c, d), pt(b, c, d))


def in_poly(poly, p):
    x, y = p
    inside = False
    for (x0, y0), (x1, y1) in zip(poly, poly[1:] + poly[:1]):
        if (y0 > y) != (y1 > y) and x < x0 + (y - y0) * (x1 - x0) / (y1 - y0):
            inside = not inside
    return inside


class Shape:
    __slots__ = ('core', 'r', 'net', 'kind', 'clear', 'bbox', 'id', 'drill', 'layer', 'width', 'name')

    def __init__(self, core, r, net, kind, clear, layer, drill=0.0, width=0.0, name=''):
        self.core, self.r, self.net, self.kind, self.clear, self.layer = core, r, net, kind, clear, layer
        self.drill, self.width, self.name = drill, width, name
        xs, ys = [p[0] for p in core], [p[1] for p in core]
        self.bbox = (min(xs) - r, min(ys) - r, max(xs) + r, max(ys) + r)


def core_dist(a, b):
    A, B = a.core, b.core
    if len(A) >= 3 and in_poly(A, B[0]) or len(B) >= 3 and in_poly(B, A[0]):
        return 0.0
    ea = list(zip(A, A[1:] + A[:1])) if len(A) >= 3 else [(A[0], A[-1])]
    eb = list(zip(B, B[1:] + B[:1])) if len(B) >= 3 else [(B[0], B[-1])]
    return min(seg_dist(p, q, r, s) for p, q in ea for r, s in eb)


def gap(a, b):
    return core_dist(a, b) - a.r - b.r


def closest_on(s, o):
    """About where s comes closest to o: the nearest of 21 samples along a segment core (else its points)."""
    A = s.core
    if len(A) == 2:
        A = [(A[0][0] + (A[1][0] - A[0][0]) * k / 20, A[0][1] + (A[1][1] - A[0][1]) * k / 20) for k in range(21)]
    return min(A, key=lambda p: core_dist(Shape([p], 0.0, None, 'point', 0.0, s.layer), o))


def pad_shape(p, layer, clear):
    x, y, w, h = p['x'], p['y'], p['w'], p['h']
    ang = math.radians(p.get('angle', 0))
    shape = p.get('shape', 'rect')
    if shape == 'circle':
        core, r = [(0.0, 0.0)], w / 2
    elif shape == 'oval':
        if abs(w - h) < 1e-6:
            core, r = [(0.0, 0.0)], w / 2
        elif w > h:
            core, r = [(-(w - h) / 2, 0.0), ((w - h) / 2, 0.0)], h / 2
        else:
            core, r = [(0.0, -(h - w) / 2), (0.0, (h - w) / 2)], w / 2
    elif shape == 'roundrect':
        r = min(p.get('rr', 0), w / 2, h / 2)
        ax, ay = w / 2 - r, h / 2 - r
        core = [(-ax, -ay), (ax, -ay), (ax, ay), (-ax, ay)]
    else:
        r = 0.0
        core = [(-w / 2, -h / 2), (w / 2, -h / 2), (w / 2, h / 2), (-w / 2, h / 2)]
    s, c = math.sin(ang), math.cos(ang)
    core = [(x + qx * c + qy * s, y - qx * s + qy * c) for qx, qy in core]
    drill = max(p.get('drill_w', 0), p.get('drill_h', 0), p.get('drill', 0))
    return Shape(core, r, p.get('net') or None, 'pad', clear, layer, drill=drill,
                 name=f"{p.get('original_ref')}.{p.get('original_pin')}")


class Grid:
    def __init__(self, cell=1.0):
        self.cell, self.cells = cell, collections.defaultdict(list)

    def keys(self, bb):
        c = self.cell
        for i in range(math.floor(bb[0] / c), math.floor(bb[2] / c) + 1):
            for j in range(math.floor(bb[1] / c), math.floor(bb[3] / c) + 1):
                yield i, j

    def add(self, s):
        for k in self.keys(s.bbox):
            self.cells[k].append(s)

    def near(self, s, margin):
        bb = (s.bbox[0] - margin, s.bbox[1] - margin, s.bbox[2] + margin, s.bbox[3] + margin)
        seen = set()
        for k in self.keys(bb):
            for o in self.cells.get(k, ()):
                if o.id not in seen and o.bbox[0] <= bb[2] and o.bbox[2] >= bb[0] and o.bbox[1] <= bb[3] and o.bbox[3] >= bb[1]:
                    seen.add(o.id)
                    yield o


def check(m, layout):
    layers = m['layers']
    cold = m.get('cold_routing', {})
    planes = cold.get('planes', {})
    plane_layers = set(cold.get('plane_layers', [])) | set(planes.values())
    classes = m['classes']
    default = classes.get('Default', m['rules'])
    net_class = m.get('net_class', {})
    cls = lambda net: classes.get(net_class.get(net, 'Default'), default)
    class_clear = lambda net: cls(net)[1] if net else default[1]
    pair_clear = cold.get('class_clearance', {})  # between two nets of the same class
    cname = lambda net: net_class.get(net, 'Default') if net else None
    gr = m.get('global_rules', {})
    hole_clear, edge_clear = gr.get('min_hole_clearance', 0.25), gr.get('min_copper_edge_clearance', 0.5)
    hole_to_hole = gr.get('min_hole_to_hole', 0.25)

    shapes = []
    by_layer = {l: Grid() for l in layers}

    def add(s):
        s.id = len(shapes)
        shapes.append(s)
        by_layer[s.layer].add(s)
        return s

    # copper of the placed board
    pad_ids = collections.defaultdict(list)  # (net) -> [(pad index, shape ids)]
    for i, p in enumerate(m['pads']):
        if p.get('npth'):
            continue
        lay = [layers[0]] if p['layers'] == 'F' else [layers[-1]] if p['layers'] == 'B' else layers
        ids = []
        for l in lay:
            ids.append(add(pad_shape(p, l, max(p.get('clearance', 0), class_clear(p.get('net'))))).id)
        if p.get('net'):
            pad_ids[p['net']].append((i, ids))
    zone_pads = []  # a zone joins the pads the model lists for it (its pieces stop short of them; the fill reaches)
    for z in m.get('zones', []):
        if z.get('layer') not in by_layer or not z.get('net'):
            continue
        for pc in z.get('pieces') or [z['points']]:
            zone_pads.append((add(Shape([tuple(q) for q in pc], 0.0, z['net'], 'zone', class_clear(z['net']), z['layer'])).id,
                              z.get('pads', [])))
    fixed = len(shapes)

    # routed copper
    narrow = 0
    track_len = collections.Counter()
    for t in layout['tracks']:
        net = t['net']
        if t['width'] < cls(net)[0] - 1e-6:
            narrow += 1
        for a, b in zip(t['pts'], t['pts'][1:]):
            if a == b:
                continue
            add(Shape([tuple(a), tuple(b)], t['width'] / 2, net, 'track', class_clear(net), t['layer'], width=t['width']))
            track_len[(net in planes, t['layer'] in plane_layers)] += math.dist(a, b)
    vias = []
    for v in layout['vias']:
        ids = [add(Shape([(v['x'], v['y'])], v['d'] / 2, v['net'], 'via', class_clear(v['net']), l, drill=v['drill'])).id for l in layers]
        vias.append((v, ids))

    # connectivity
    parent = list(range(len(shapes)))
    plane_node = {}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        a, b = find(i), find(j)
        if a != b:
            parent[a] = b

    for s in shapes:
        for o in by_layer[s.layer].near(s, TOL):
            if o.id > s.id and o.net == s.net and s.net and gap(s, o) <= TOL:
                union(s.id, o.id)
    for _, ids in [x for v in pad_ids.values() for x in v] + vias:
        for a in ids[1:]:
            union(ids[0], a)
    pad_by_name = {f"{p.get('original_ref')}.{p.get('original_pin')}": ids for v in pad_ids.values() for i, ids in v
                   for p in [m['pads'][i]]}
    for zid, names in zone_pads:
        for n in names:
            if n in pad_by_name:
                union(zid, pad_by_name[n][0])

    regions = collections.defaultdict(list)
    for r in cold.get('plane_regions', []):
        regions[(r['net'], r['layer'])].append(([tuple(q) for q in r['points']], [[tuple(q) for q in h] for h in r.get('holes', [])]))
    outline = [tuple(q) for q in m['outline']]

    def joins_plane(net, layer, xy):
        regs = regions.get((net, layer))
        if not regs:
            return in_poly(outline, xy)
        return any(in_poly(o, xy) and not any(in_poly(h, xy) for h in hs) for o, hs in regs)

    for s in shapes:
        if s.net in planes and (s.kind == 'via' or (s.kind == 'pad' and s.drill > 0)) and s.layer == planes[s.net]:
            xy = s.core[0] if len(s.core) == 1 else (sum(p[0] for p in s.core) / len(s.core), sum(p[1] for p in s.core) / len(s.core))
            if joins_plane(s.net, s.layer, xy):
                union(s.id, plane_node.setdefault(s.net, s.id))

    nets_total = nets_complete = open_groups = 0
    incomplete = []
    for net, pads in pad_ids.items():
        if len(pads) < 2 and net not in planes:
            continue
        groups = {find(ids[0]) for _, ids in pads}
        nets_total += 1
        if len(groups) == 1:
            nets_complete += 1
        else:
            open_groups += len(groups) - 1
            incomplete.append(dict(net=net, groups=len(groups), pads=len(pads)))

    # clearance (routed copper against everything of another net), edge, drills
    violations = collections.Counter()
    examples, overlaps = [], []
    for s in shapes[fixed:]:
        for o in by_layer[s.layer].near(s, 0.6):
            if o.net == s.net and s.net is not None:
                continue
            if o.id >= fixed and o.id < s.id:
                continue
            if s.kind == 'via' and o.kind == 'via' and s.core == o.core:
                continue
            need = max(s.clear, o.clear)
            if s.net and o.net and cname(s.net) == cname(o.net) and cname(s.net) in pair_clear:
                need = max(need, pair_clear[cname(s.net)])
            g = gap(s, o)
            if g < need - TOL:
                k = f'{s.kind}-{o.kind}'
                violations[k] += 1
                if len(examples) < 20:
                    examples.append(dict(kind=k, layer=s.layer, nets=[s.net, o.net], gap=round(g, 4), need=need,
                                         at=[round(c, 3) for c in s.core[0]], other=o.name or None))
                if g < -TOL:  # touching: a short, or copper through an obstacle of no net
                    overlaps.append(dict(kind=k, layer=s.layer, nets=[s.net, o.net], gap=round(g, 4),
                                         at=[round(c, 3) for c in closest_on(s, o)], other=o.name or None))
    # drills: via to other-net copper (hole clearance from the drill edge) and to other drills (any net)
    holes = Grid()
    for p in m['pads']:
        d = max(p.get('drill_w', 0), p.get('drill_h', 0), p.get('drill', 0))
        if d > 0:
            s = Shape([(p['x'], p['y'])], d / 2, p.get('net') or None, 'hole', 0, layers[0])
            s.id = -1 - len(holes.cells)
            holes.add(s)
    hole_list = []
    for v, ids in vias:
        s = Shape([(v['x'], v['y'])], v['drill'] / 2, v['net'], 'hole', 0, layers[0])
        s.id = 10 ** 9 + len(hole_list)
        for o in holes.near(s, hole_to_hole):
            if o.core != s.core and gap(s, o) < hole_to_hole - TOL:
                violations['drill-drill'] += 1
        holes.add(s)
        hole_list.append(s)
        for l in layers:
            for o in by_layer[l].near(s, hole_clear):
                if o.net != v['net'] and o.kind != 'via' and gap(s, o) < hole_clear - TOL:
                    violations['drill-copper'] += 1
    # board edge
    edges = list(zip(outline, outline[1:] + outline[:1]))
    edge_viol = 0
    for s in shapes[fixed:]:
        if s.kind == 'via' and s.layer != layers[0]:
            continue
        for a, b in edges:
            if seg_dist(s.core[0], s.core[-1], a, b) - s.r < edge_clear - TOL:
                edge_viol += 1
                break
    violations['edge'] = edge_viol

    return dict(nets=nets_total, complete_nets=nets_complete, open_pad_groups=open_groups,
                incomplete=sorted(incomplete, key=lambda x: x['net']),
                clearance_violations=dict(violations), clearance_total=sum(violations.values()),
                overlaps=len(overlaps), overlap_list=overlaps,
                examples=examples, tracks=len(layout['tracks']), vias=len(layout['vias']),
                track_length_mm=round(sum(track_len.values()), 1),
                signal_track_on_plane_layers_mm=round(track_len[(False, True)], 1),
                tracks_below_class_width=narrow)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('model', type=pathlib.Path)
    p.add_argument('layout', type=pathlib.Path)
    p.add_argument('--out', type=pathlib.Path)
    a = p.parse_args()
    r = check(json.loads(a.model.read_text()), json.loads(a.layout.read_text()))
    if a.out:
        a.out.write_text(json.dumps(r, indent=1))
    print(json.dumps({k: v for k, v in r.items() if k not in ('incomplete', 'examples', 'overlap_list')}))


if __name__ == '__main__':
    main()
