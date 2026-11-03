#!/usr/bin/env python3
"""topor: read a board's inner-layer pours as plane regions. Run with KiCad Python on a case input dir.

An inner layer whose pours cover at least half the board is a plane layer: no signal routing there, and
a via of a pour's net joins that net's plane only inside the pour (clear of higher-priority pours of
other nets and of fill keepouts, e.g. the cutout under an antenna). Partial pours (a power island on a
shared layer) are regions like any other. Written into model.json cold_routing: planes (net -> main
layer), plane_regions, fill_keepouts.
"""
import json, pathlib, sys
import pcbnew

case = pathlib.Path(sys.argv[1])
m = json.loads((case / 'model.json').read_text())
b = pcbnew.LoadBoard(str(case / 'unrouted.kicad_pcb'))
ox, oy = m['origin']
MM = pcbnew.ToMM
inner = m['layers'][1:-1]
board_area = m['size'][0] * m['size'][1]


def pts(chain):
    return [[MM(chain.CPoint(i).x) - ox, MM(chain.CPoint(i).y) - oy] for i in range(chain.PointCount())]


regions, keepouts = [], []
cover = {l: 0.0 for l in inner}
for z in b.Zones():
    ls = [b.GetStandardLayerName(l) for l in z.GetLayerSet().CuStack()]
    o = z.Outline()
    for l in ls:
        if l not in inner:
            continue
        if z.GetIsRuleArea():
            if z.GetDoNotAllowZoneFills():
                for i in range(o.OutlineCount()):
                    keepouts.append(dict(layer=l, points=pts(o.Outline(i))))
            continue
        if not z.GetNetname():
            continue
        for i in range(o.OutlineCount()):
            regions.append(dict(net=z.GetNetname(), layer=l, priority=z.GetAssignedPriority(), points=pts(o.Outline(i)),
                                holes=[pts(o.Hole(i, j)) for j in range(o.HoleCount(i))]))
        cover[l] += o.Area() / 1e12
plane_layers = [l for l in inner if cover[l] >= 0.5 * board_area]
regions = [r for r in regions if r['layer'] in plane_layers]
keepouts = [k for k in keepouts if k['layer'] in plane_layers]


def area(p):
    return abs(sum(p[i][0] * p[(i + 1) % len(p)][1] - p[(i + 1) % len(p)][0] * p[i][1] for i in range(len(p)))) / 2


best = {}
for r in regions:
    a = area(r['points'])
    if r['net'] not in best or a > best[r['net']][1]:
        best[r['net']] = (r['layer'], a)
cr = m.setdefault('cold_routing', {})
cr['planes'] = {n: l for n, (l, _) in best.items()}
cr['plane_layers'] = plane_layers
cr['plane_regions'] = regions
cr['fill_keepouts'] = keepouts
(case / 'model.json').write_text(json.dumps(m, indent=2))
print(json.dumps(dict(plane_layers=plane_layers, coverage={l: round(cover[l] / board_area, 2) for l in inner},
                      planes=cr['planes'], regions=len(regions), fill_keepouts=len(keepouts))))
