#!/usr/bin/env python3
"""Convert a prepared topor case (model.json) to tscircuit's SimpleRouteJson.

The board is given to tscircuit as it is given to topor: every copper layer, every pad as an obstacle (rotated
rectangles; circles as circles), rule areas, small copper zones and plane pours as obstacles, the board outline
with its edge clearance. Coordinates are flipped to y up, so KiCad's orientations (counter-clockwise as seen)
are tscircuit's ccwRotationDegrees.

Planes: SimpleRouteJson has no plane drop, so each SMD pad of a plane net becomes a two-point connection from
the pad to the same place on its plane layer (inside its pour); the router has to bring a via down, which is
what topor's plane drop does. Through pads of plane nets join their plane in the KiCad fill and need nothing.
Plane pours are copper-pour obstacles of their net; the rest of a plane layer is kept free of signal tracks
by filler pours of no net.

Not representable, and listed in meta.json: per-class clearance (the default class clearance is global; a pad
with a larger clearance is grown by the difference), via-only keepouts, per-layer widths, neck-downs, RF
single-layer nets and differential pairs (the models list none).
"""
import argparse
import json
import math
import pathlib

STRIP = 0.25  # mm, height of the rectangles a polygon is filled with


def layer_names(layers):
    n = len(layers)
    return {l: ('top' if i == 0 else 'bottom' if i == n - 1 else f'inner{i}') for i, l in enumerate(layers)}


def intervals_at(rings, y):
    """x intervals inside the even-odd union of rings on the horizontal line y."""
    xs = []
    for ring in rings:
        for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]):
            if (y0 <= y < y1) or (y1 <= y < y0):
                xs.append(x0 + (y - y0) * (x1 - x0) / (y1 - y0))
    xs.sort()
    return [(xs[i], xs[i + 1]) for i in range(0, len(xs) - 1, 2)]


def merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1] + 1e-9:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def subtract(iv, cut):
    out = []
    for a, b in iv:
        segs = [(a, b)]
        for c, d in cut:
            segs = [s for x, y in segs for s in ((x, min(y, c)), (max(x, d), y)) if s[1] - s[0] > 1e-6]
        out += segs
    return out


def strips(rings, cover=True):
    """Rows of (y0, y1, intervals) filling the rings. cover=True over-covers (obstacles, keepouts): each row is
    the union of the scanlines at its edges and middle; otherwise the intersection (pours stay inside)."""
    ys = [p[1] for r in rings for p in r]
    if not ys:
        return []
    lo, hi = min(ys), max(ys)
    rows = []
    n = max(1, math.ceil((hi - lo) / STRIP))
    h = (hi - lo) / n
    for k in range(n):
        y0, y1 = lo + k * h, lo + (k + 1) * h
        samples = [intervals_at(rings, y) for y in (y0 + 1e-6, (y0 + y1) / 2, y1 - 1e-6)]
        if cover:
            iv = merge([i for s in samples for i in s])
        else:
            iv = samples[0]
            for s in samples[1:]:
                iv = [(max(a, c), min(b, d)) for a, b in iv for c, d in s if min(b, d) - max(a, c) > 1e-6]
        rows.append((y0, y1, iv))
    return rows


def rects(rows, layers, connected, **extra):
    out = []
    for y0, y1, iv in rows:
        for a, b in iv:
            if b - a < 1e-3:
                continue
            out.append(dict(type='rect', layers=layers, center=dict(x=(a + b) / 2, y=-(y0 + y1) / 2),
                            width=b - a, height=y1 - y0, connectedTo=connected, **extra))
    return out


def ring(pts):
    return [(float(p[0]), float(p[1])) for p in pts]


def inside(rings, x, y):
    return any(a <= x <= b for a, b in intervals_at(rings, y))


def convert(m):
    layers = m['layers']
    lname = layer_names(layers)
    all_layers = [lname[l] for l in layers]
    cold = m.get('cold_routing', {})
    planes = cold.get('planes', {})
    plane_layers = sorted(set(cold.get('plane_layers', [])) | set(planes.values()), key=layers.index)
    signal = [l for l in layers if l not in plane_layers]
    gr = m.get('global_rules', {})
    classes = m['classes']
    default = classes.get('Default', m['rules'])
    c0 = default[1]
    net_class = m.get('net_class', {})
    cls = lambda net: classes.get(net_class.get(net, 'Default'), default)
    hole_clear = gr.get('min_hole_clearance', 0.25)
    edge_clear = gr.get('min_copper_edge_clearance', 0.5)
    unsupported = {}

    def note(k, n=1):
        unsupported[k] = unsupported.get(k, 0) + n

    # nets, connections
    pads_by_net = {}
    for p in m['pads']:
        if p.get('net') and not p.get('npth'):
            pads_by_net.setdefault(p['net'], []).append(p)
    # plane layers: each net's pour (higher priority first, without fill keepouts and rule areas), the rest of
    # the layer left to pours of no net so no signal runs there
    outline = ring(m['outline'])
    pours = {}  # (net, layer) -> rows; net None: the rest
    for layer in plane_layers:
        regs = [r for r in cold.get('plane_regions', []) if r['layer'] == layer]
        if not regs:
            regs = [dict(net=n, points=m['outline'], holes=[], priority=0) for n, l in planes.items() if l == layer][:1]
        regs.sort(key=lambda r: -r.get('priority', 0))
        keep = [ring(k['points']) for k in cold.get('fill_keepouts', []) if k.get('layer') == layer]
        keep += [ring(pg) for a in m.get('rule_areas', []) if a.get('tracks', True)
                 and (layer in a.get('layers', []) or not a.get('layers')) for pg in a.get('polygons', [])]
        rows_outline = strips([outline], cover=False)
        taken = {k: [] for k in range(len(rows_outline))}
        for r in regs:
            rr = [ring(r['points'])] + [ring(h) for h in r.get('holes', [])]
            rows = []
            for k, (y0, y1, _) in enumerate(rows_outline):
                iv = intervals_at(rr, (y0 + y1) / 2)
                lo, hi = intervals_at(rr, y0 + 1e-6), intervals_at(rr, y1 - 1e-6)
                iv = [(max(a, c, e), min(b, d, f)) for a, b in iv for c, d in lo for e, f in hi if min(b, d, f) - max(a, c, e) > 1e-6]
                for kr in keep:
                    iv = subtract(iv, merge(intervals_at([kr], y0 + 1e-6) + intervals_at([kr], (y0 + y1) / 2) + intervals_at([kr], y1 - 1e-6)))
                iv = subtract(merge(iv), merge(taken[k]))
                taken[k] += iv
                rows.append((y0, y1, iv))
            pours.setdefault((r['net'], layer), []).extend(rows)
        pours[(None, layer)] = [(y0, y1, subtract(iv, merge(taken[k]))) for k, (y0, y1, iv) in enumerate(rows_outline)]

    def point(p):
        x, y = p['x'], -p['y']
        if p['layers'] == 'F':
            return dict(x=x, y=y, layer=lname[layers[0]])
        if p['layers'] == 'B':
            return dict(x=x, y=y, layer=lname[layers[-1]])
        return dict(x=x, y=y, layers=[lname[l] for l in signal])

    def drop_target(net, layer, x, y, margin=0.4):
        """Where on the plane layer a pad's drop lands: the nearest point at least margin inside its net's pour
        (under the pad when it can)."""
        best = None
        for y0, y1, iv in pours.get((net, layer), []):
            for a, b in iv:
                if b - a < 2 * margin:
                    continue
                cx, cy = min(max(x, a + margin), b - margin), (y0 + y1) / 2
                d = math.hypot(cx - x, cy - y)
                if best is None or d < best[0]:
                    best = (d, cx, cy)
        return (best[1], best[2]) if best else (x, y)

    connections, conn_net, net_conns = [], {}, {}
    for net, pads in sorted(pads_by_net.items()):
        width = cls(net)[0]
        if net in planes:
            for p in pads:
                if p['layers'] == 'all':
                    continue
                tx, ty = drop_target(net, planes[net], p['x'], p['y'])
                name = f'c{len(connections)}'
                connections.append(dict(name=name, nominalTraceWidth=width,
                                        pointsToConnect=[point(p), dict(x=tx, y=-ty, layer=lname[planes[net]])]))
                conn_net[name] = dict(net=net, kind='plane_drop', pad=f"{p['original_ref']}.{p['original_pin']}")
                net_conns.setdefault(net, []).append(name)
        elif len(pads) >= 2:
            name = f'c{len(connections)}'
            connections.append(dict(name=name, nominalTraceWidth=width, pointsToConnect=[point(p) for p in pads]))
            conn_net[name] = dict(net=net, kind='signal')
            net_conns.setdefault(net, []).append(name)

    def tied(net):
        return [net] + net_conns.get(net, []) if net else []

    obstacles = []
    for p in m['pads']:
        x, y, w, h = p['x'], -p['y'], p['w'], p['h']
        lay = [lname[layers[0]]] if p['layers'] == 'F' else [lname[layers[-1]]] if p['layers'] == 'B' else all_layers
        if p.get('npth'):
            d = max(p.get('drill_w', 0), p.get('drill_h', 0), p.get('drill', 0))
            obstacles.append(dict(type='rect', layers=all_layers, center=dict(x=x, y=y), width=d, height=d, shape='circle',
                                  isNonPlatedHole=True, connectedTo=[]))
            continue
        drill = max(p.get('drill_w', 0), p.get('drill_h', 0), p.get('drill', 0))
        clear = p.get('clearance', c0)
        if drill > 0:
            ring_w = max(0.0, min(w - max(p.get('drill_w', 0), p.get('drill', 0)), h - max(p.get('drill_h', 0), p.get('drill', 0))) / 2)
            clear = max(clear, hole_clear - ring_w)
        grow = max(0.0, clear - c0)
        o = dict(type='rect', layers=lay, center=dict(x=x, y=y), width=w + 2 * grow, height=h + 2 * grow,
                 connectedTo=tied(p.get('net')))
        if p['shape'] == 'circle':
            o['shape'] = 'circle'
        elif p.get('angle'):
            o['ccwRotationDegrees'] = p['angle']
        obstacles.append(o)

    for z in m.get('zones', []):
        if z.get('layer') not in lname:
            continue
        pieces = z.get('pieces') or [z['points']]
        for pc in pieces:
            obstacles += rects(strips([ring(pc)]), [lname[z['layer']]], tied(z.get('net')))

    # rule areas: 'tracks' and 'vias' are KiCad's do-not-allow flags; an area that forbids neither only keeps
    # pours out (already in the plane pours above)
    for a in m.get('rule_areas', []):
        lay = [lname[l] for l in a.get('layers', []) if l in lname] or all_layers
        if a.get('tracks', True):
            for pg in a.get('polygons', []):
                obstacles += rects(strips([ring(pg)]), lay, [])
        elif a.get('vias', True):
            note('via-only keepout areas (tscircuit has no via keepout)')

    for g in m.get('copper_graphics', []):
        if g.get('layer') in lname:
            obstacles += rects(strips([ring(g['points'])] + [ring(h) for h in g.get('holes', [])]), [lname[g['layer']]], tied(g.get('net')))
    for g in m.get('edge_graphics', []) + [dict(points=h) for h in m.get('outline_holes', [])]:
        obstacles += rects(strips([ring(g['points'])]), all_layers, [])

    for (net, layer), rows in pours.items():
        obstacles += rects(rows, [lname[layer]], tied(net), isCopperPour=True)

    xs, ys = [p[0] for p in outline], [-p[1] for p in outline]
    vd, vdrill = default[2], default[3]
    srj = dict(
        layerCount=len(layers),
        minTraceWidth=min(c[0] for c in classes.values()),
        nominalTraceWidth=default[0],
        minViaDiameter=vd, minViaPadDiameter=vd, min_via_pad_diameter=vd,
        minViaHoleDiameter=vdrill, min_via_hole_diameter=vdrill,
        defaultObstacleMargin=c0,
        minTraceToPadEdgeClearance=c0,
        minViaEdgeToPadEdgeClearance=c0,
        minTraceToHoleEdgeClearance=max(hole_clear, edge_clear),
        minViaHoleEdgeToViaHoleEdgeClearance=gr.get('min_hole_to_hole', 0.25),
        minPlatedHoleDrillEdgeToDrillEdgeClearance=gr.get('min_hole_to_hole', 0.25),
        minBoardEdgeClearance=edge_clear,
        allowViaInPad=False,
        bounds=dict(minX=min(xs), maxX=max(xs), minY=min(ys), maxY=max(ys)),
        outline=[dict(x=x, y=y) for x, y in zip(xs, ys)],
        obstacles=obstacles,
        connections=connections,
    )
    if len(classes) > 1:
        note('per-class clearance (global default used; larger-clearance pads grown)')
    for k in ('layer_widths', 'layer_min_widths', 'disallow', 'track_clearance'):
        if cold.get(k):
            note(f'cold_routing.{k}')
    meta = dict(layers=layers, layer_names=lname, plane_layers=plane_layers, planes=planes, connections=conn_net,
                classes=classes, net_class=net_class, origin=m.get('origin'), unsupported=unsupported,
                counts=dict(obstacles=len(obstacles), connections=len(connections),
                            signal_connections=sum(1 for c in conn_net.values() if c['kind'] == 'signal'),
                            plane_drops=sum(1 for c in conn_net.values() if c['kind'] == 'plane_drop')))
    return srj, meta


def main():
    p = argparse.ArgumentParser()
    p.add_argument('input', type=pathlib.Path, help='case directory with model.json')
    p.add_argument('--out', type=pathlib.Path, required=True)
    a = p.parse_args()
    srj, meta = convert(json.loads((a.input / 'model.json').read_text()))
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / 'srj.json').write_text(json.dumps(srj))
    (a.out / 'srj-meta.json').write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta['counts']), json.dumps(meta['unsupported']))


if __name__ == '__main__':
    main()
