#!/usr/bin/env python3
"""topor: merge a board's overlapping rule areas into their union, so each keeps its true shape. Run with KiCad
Python on a case input dir.

The engine joins overlapping obstacles into their convex hull (touching copper of one pad stack, say). That is
wrong for keepouts drawn as many overlapping strips, e.g. a band along a round board edge built of quads: their
hull covers the whole board, and the band itself stops keeping anything out. Here the rule areas that forbid the
same things on the same layers are united with KiCad's own polygon operations (the shape KiCad's DRC tests
against), and model.json `rule_areas` is rewritten from the board: one entry per united outline. Only where the
hull is much larger than the union (over 1.5 times its area) is that worth it; overlapping areas whose hull is
about their own shape stay as drawn, and a board with none of the first kind keeps its rule areas untouched. A
union with a hole (a closed ring) cannot be one engine obstacle: its areas stay as drawn, with a warning.
"""
import json, pathlib, sys
import pcbnew

case = pathlib.Path(sys.argv[1])
m = json.loads((case / 'model.json').read_text())
b = pcbnew.LoadBoard(str(case / 'unrouted.kicad_pcb'))
ox, oy = m['origin']
MM = pcbnew.ToMM


def ring(chain):
    return [[MM(chain.CPoint(i).x) - ox, MM(chain.CPoint(i).y) - oy] for i in range(chain.PointCount())]


groups = {}
for z in b.Zones():
    if not z.GetIsRuleArea():
        continue
    layers = tuple(b.GetStandardLayerName(l) for l in z.GetLayerSet().CuStack())
    key = (layers, bool(z.GetDoNotAllowTracks()), bool(z.GetDoNotAllowVias()))
    if not layers or not (key[1] or key[2]):
        continue
    groups.setdefault(key, []).append(z)


def area(pts):
    return abs(sum(pts[i][0] * pts[(i + 1) % len(pts)][1] - pts[(i + 1) % len(pts)][0] * pts[i][1] for i in range(len(pts)))) / 2


def hull(pts):
    pts = sorted(map(tuple, pts))
    if len(pts) < 3:
        return pts
    cross = lambda o, a, b: (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, hi = [], []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(hi) >= 2 and cross(hi[-2], hi[-1], p) <= 0:
            hi.pop()
        hi.append(p)
    return lo[:-1] + hi[:-1]


areas, warnings, united = [], [], 0
for (layers, tracks, vias), zs in groups.items():
    u = pcbnew.SHAPE_POLY_SET()
    for z in zs:
        u.BooleanAdd(z.Outline())
    u.Simplify()
    for i in range(u.OutlineCount()):
        out = ring(u.Outline(i))
        su = pcbnew.SHAPE_POLY_SET(u.Outline(i))
        part = [z for z in zs if su.Contains(z.Outline().BBox().Centre())]
        drawn = [ring(z.Outline().Outline(k)) for z in part for k in range(z.Outline().OutlineCount())] or [out]
        if u.HoleCount(i):
            warnings.append('%d rule areas on %s unite to a ring with a hole: kept as drawn' % (len(part), ','.join(layers)))
            polys = drawn
        elif len(part) > 1 and area(hull(out)) > 1.5 * area(out):
            polys = [out]
            united += 1
        else:
            polys = drawn
        for pg in polys:
            areas.append(dict(layers=list(layers), tracks=tracks, vias=vias, polygons=[pg]))
before = len(m.get('rule_areas', []))
if united:
    m['rule_areas'] = areas
m['rule_areas_merged'] = True
(case / 'model.json').write_text(json.dumps(m, indent=2))
print(json.dumps(dict(rule_areas_before=before, after=len(m['rule_areas']), united=united, points=[len(a['polygons'][0]) for a in areas] if united else [], warnings=warnings)))
