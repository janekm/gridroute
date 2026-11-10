#!/usr/bin/env python3
"""Turn one of tscircuit's own boards (a SimpleRouteJson) into a topor case, the reverse of to_srj.py.

    srj_to_model.py SRJ.json OUT_DIR [--clearance 0.15] [--width W] [--via 0.3/0.15] [--edge 0.2]

OUT_DIR gets model.json for topor, the SRJ itself as srj.json (compare.py gives tscircuit that file, not a
conversion back) and srj-meta.json (connections to nets, the coordinate shift) for from_srj.py.

The board is the SRJ bounds, y flipped to KiCad's y down with the origin at the bounds' top left. Each obstacle
becomes a pad: on one layer (F/B), or on every layer ('all', no drill: tscircuit gives holes as rectangular
copper obstacles, and they stay that way). Connections that share an obstacle (connectedTo) are one net. An
obstacle is a pad of the net of the connections it lists, else of its connectivity net (a pad nothing is routed
to), else of no net.

The SRJ gives a minimum trace width and nothing else, so the rules are set here: tscircuit's own via (0.3 mm,
0.15 mm drill) and a clearance (0.15 mm by default; tscircuit's own DRC asks for 0.1 mm).
"""
import argparse
import collections
import json
import pathlib
import shutil


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('srj', type=pathlib.Path)
    ap.add_argument('out', type=pathlib.Path)
    ap.add_argument('--clearance', type=float, default=0.15)
    ap.add_argument('--width', type=float, help='track width (default: the SRJ minTraceWidth)')
    ap.add_argument('--via', default='0.3/0.15', help='via diameter/drill')
    ap.add_argument('--edge', type=float, default=0.2, help='copper to board edge')
    a = ap.parse_args()
    srj = json.loads(a.srj.read_text())
    if srj.get('layerCount', 2) != 2:
        ap.error('two-layer boards only')
    b = srj['bounds']
    x0, y1 = b['minX'], b['maxY']
    width = a.width or srj.get('minTraceWidth', 0.15)
    vd, vdrill = (float(v) for v in a.via.split('/'))

    conns = {c['name']: c for c in srj['connections']}
    par = {n: n for n in conns}

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    obstacles = srj['obstacles']
    for o in obstacles:
        cs = [c for c in o.get('connectedTo', []) if c in conns]
        for c in cs[1:]:
            par[find(c)] = find(cs[0])
    # a net is named after the connectivity net its obstacles carry, else its first connection
    named = {}
    for o in obstacles:
        cs = [c for c in o.get('connectedTo', []) if c in conns]
        cn = [c for c in o.get('connectedTo', []) if c.startswith('connectivity_net')]
        if cs and cn:
            named.setdefault(find(cs[0]), cn[0])
    roots = {find(n) for n in conns}
    net_name = {r: named.get(r, r) for r in roots}
    if len(set(net_name.values())) != len(net_name):  # two nets on one connectivity net: keep them apart
        net_name = {r: f'{named.get(r, r)}~{r}' for r in roots}

    pads, notes = [], collections.Counter()
    for i, o in enumerate(obstacles):
        L = set(o['layers'])
        lay = 'all' if {'top', 'bottom'} <= L else 'F' if L == {'top'} else 'B' if L == {'bottom'} else None
        if o['type'] != 'rect' or lay is None:
            notes[f"skipped: {o['type']} on {sorted(L)}"] += 1
            continue
        ids = o.get('connectedTo', [])
        cs = [c for c in ids if c in conns]
        cn = [c for c in ids if c.startswith('connectivity_net')]
        net = net_name[find(cs[0])] if cs else cn[0] if cn else ''
        notes['pads of a routed net' if cs else 'pads of a net nothing is routed to' if cn else 'obstacles of no net'] += 1
        ref = next((c for c in ids if c.startswith('pcb_smtpad') or c.startswith('pcb_plated')), f'obstacle{i}')
        pads.append(dict(ref=f'pad_{i}', original_ref=ref, original_pin='1', num='1', x=o['center']['x'] - x0,
                         y=y1 - o['center']['y'], w=o['width'], h=o['height'], rot=0,
                         angle=o.get('ccwRotationDegrees', 0) or 0, shape='rect', polygons=[], layers=lay, npth=False,
                         drill=0.0, drill_w=0.0, drill_h=0.0, net=net, clearance=a.clearance, rr=0))

    # every point to connect should land on a pad of its connection's net
    side = {'top': 'F', 'bottom': 'B'}
    unmatched = 0
    for c in srj['connections']:
        net = net_name[find(c['name'])]
        for p in c['pointsToConnect']:
            x, y, s = p['x'] - x0, y1 - p['y'], side.get(p.get('layer'))
            if not any(q['net'] == net and (s is None or q['layers'] in (s, 'all'))
                       and abs(q['x'] - x) <= q['w'] / 2 + 1e-6 and abs(q['y'] - y) <= q['h'] / 2 + 1e-6 for q in pads):
                unmatched += 1

    w_b, h_b = b['maxX'] - b['minX'], b['maxY'] - b['minY']
    count = collections.Counter(q['net'] for q in pads if q['net'])
    routable = [n for n, k in count.items() if k >= 2]
    rules = [width, a.clearance, vd, vdrill]
    model = dict(
        name=a.srj.stem, size=[w_b, h_b], origin=[0.0, 0.0],
        outline=[[0.0, 0.0], [w_b, 0.0], [w_b, h_b], [0.0, h_b]],
        layers=['F.Cu', 'B.Cu'], rules=rules, pads=pads, ignored_pads=[], footprint_count=0,
        routable_nets=len(routable), required_pad_connections=sum(count[n] - 1 for n in routable),
        differential_pairs=[], copper_graphics=[], mask_graphics=[], edge_graphics=[], outline_holes=[],
        classes={'Default': rules}, net_class={n: 'Default' for n in sorted(count)}, rule_areas=[], rf_nets=[],
        global_rules=dict(min_clearance=a.clearance, min_track_width=width, min_via_diameter=vd,
                          min_through_hole_diameter=vdrill, min_via_annular_width=(vd - vdrill) / 2,
                          min_hole_clearance=a.clearance, min_hole_to_hole=0.25, min_copper_edge_clearance=a.edge),
        cold_routing=dict(planes={}, layer_widths={}, layer_min_widths={}, disallow={}, class_clearance={}),
        source=dict(srj=a.srj.name, converted_by='bench/pcb/tscircuit/srj_to_model.py'),
    )
    meta = dict(layers=['F.Cu', 'B.Cu'], layer_names={'F.Cu': 'top', 'B.Cu': 'bottom'}, plane_layers=[], planes={},
                shift=dict(minX=x0, maxY=y1),
                connections={c: dict(net=net_name[find(c)], kind='signal') for c in conns},
                counts=dict(obstacles=len(obstacles), connections=len(conns), signal_connections=len(conns), plane_drops=0))
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / 'model.json').write_text(json.dumps(model))
    (a.out / 'srj-meta.json').write_text(json.dumps(meta, indent=1))
    if a.srj.resolve() != (a.out / 'srj.json').resolve():
        shutil.copyfile(a.srj, a.out / 'srj.json')
    print(json.dumps(dict(pads=len(pads), nets=len(count), routable_nets=len(routable),
                          pad_connections=model['required_pad_connections'], connections=len(conns),
                          unmatched_points=unmatched, notes=dict(notes),
                          rules=dict(width=width, clearance=a.clearance, via=a.via, edge=a.edge))))


if __name__ == '__main__':
    main()
