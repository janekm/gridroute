#!/usr/bin/env python3
"""topor: read a board's inner-layer pours as plane regions. Run with KiCad Python on a case input dir.

An inner layer whose pours cover at least half the board is a plane layer: no signal routing there, and
a via of a pour's net joins that net's plane only inside the pour (clear of higher-priority pours of
other nets and of fill keepouts, e.g. the cutout under an antenna). Partial pours (a power island on a
shared layer) are regions like any other. Written into model.json cold_routing: planes (net -> main
layer), plane_regions, fill_keepouts.

--outer-pours (experimental): a plane net whose inner pours reach under half its pads (a module's
separate ground, poured on the outer layer under it) also gets its outer-layer pours as regions, so an
SMD pad in one joins the plane on its own layer. On EXT this saves GND1 some 450 mm of track, but other
nets' wires on that layer then cut the pour apart: it needs the pour kept whole first.
"""
import json, pathlib, sys
import pcbnew

case = pathlib.Path(sys.argv[1])
use_outer = '--outer-pours' in sys.argv[2:]
m = json.loads((case / 'model.json').read_text())
b = pcbnew.LoadBoard(str(case / 'unrouted.kicad_pcb'))
ox, oy = m['origin']
MM = pcbnew.ToMM
inner = m['layers'][1:-1]
board_area = m['size'][0] * m['size'][1]


def pts(chain):
    return [[MM(chain.CPoint(i).x) - ox, MM(chain.CPoint(i).y) - oy] for i in range(chain.PointCount())]


regions, keepouts, outer = [], [], []
cover = {l: 0.0 for l in inner}
for z in b.Zones():
    ls = [b.GetStandardLayerName(l) for l in z.GetLayerSet().CuStack()]
    o = z.Outline()
    for l in ls:
        if l not in inner:
            if not z.GetIsRuleArea() and z.GetNetname():
                for i in range(o.OutlineCount()):
                    outer.append(dict(net=z.GetNetname(), layer=l, priority=z.GetAssignedPriority(), points=pts(o.Outline(i)),
                                      holes=[pts(o.Hole(i, j)) for j in range(o.HoleCount(i))]))
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


def inside(p, x, y):
    c, j = False, len(p) - 1
    for i in range(len(p)):
        if (p[i][1] > y) != (p[j][1] > y) and x < (p[j][0] - p[i][0]) * (y - p[i][1]) / (p[j][1] - p[i][1]) + p[i][0]:
            c = not c
        j = i
    return c


poor = []
for net in sorted({r['net'] for r in regions}):
    pads = [p for p in m['pads'] if p.get('net') == net]
    covered = sum(1 for p in pads if any(r['net'] == net and inside(r['points'], p['x'], p['y']) for r in regions))
    if use_outer and pads and covered < 0.5 * len(pads):
        poor.append((net, covered, len(pads)))
        regions += [r for r in outer if r['net'] == net]

best = {}
# a net's plane layer is an inner one (an outer pour that serves as its plane is a region, not a layer)
for r in regions:
    if r['layer'] not in plane_layers:
        continue
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
                      planes=cr['planes'], regions=len(regions), fill_keepouts=len(keepouts),
                      outer_pours_for=[dict(net=n, pads_on_inner_pours=c, pads=t) for n, c, t in poor])))
