#!/usr/bin/env python3
"""topor: add a board's small copper zones (local islands that carry connectivity: switch nodes, antenna
feeds) to a prepared model as fixed copper. Run with KiCad Python on a case input dir.

A zone connects only what KiCad's filler makes of it, so the zones are filled here exactly as the export
will refill them (on the unrouted board, before any track exists) and the fill itself is the copper: each
island that reaches pads of its net becomes convex pieces tied to those pads. A routed board that keeps
clear of that copper refills to the same copper. A zone whose fill is empty (e.g. a thermal relief gap
wider than the zone) or reaches no pad connects nothing and is left to the router as plain space.

Big pours, inferred planes and islands of plane nets (optional copper that wires may cut) are left
out: the export refills them around the routed copper.
"""
import json, pathlib, sys
import pcbnew

case = pathlib.Path(sys.argv[1])
m = json.loads((case / 'model.json').read_text())
b = pcbnew.LoadBoard(str(case / 'unrouted.kicad_pcb'))
ox, oy = m['origin']
planes = m.get('cold_routing', {}).get('planes', {})
MM = pcbnew.ToMM
# pieces stop this far short of their own pads and of each other: touching copper, not overlapping hulls
PAD_GAP, SHRINK = 0.01, 0.003


def area(pts):
    return abs(sum(pts[i][0] * pts[(i + 1) % len(pts)][1] - pts[(i + 1) % len(pts)][0] * pts[i][1] for i in range(len(pts)))) / 2


def ring(chain):
    return [(MM(chain.CPoint(i).x) - ox, MM(chain.CPoint(i).y) - oy) for i in range(chain.PointCount())]


def convex_pieces(rings):
    """Convex pieces covering a polygon set (outlines, holes and fracture slits alike, even-odd): vertical
    slabs between vertex x positions cut it into trapezoids, which grow rightwards while they stay convex."""
    edges = []
    for r in rings:
        for i in range(len(r)):
            a, c = r[i], r[(i + 1) % len(r)]
            if abs(a[0] - c[0]) > 1e-9:
                edges.append((a, c) if a[0] < c[0] else (c, a))
    xs = sorted({round(p[0], 7) for r in rings for p in r})
    y_at = lambda e, x: e[0][1] + (e[1][1] - e[0][1]) * (x - e[0][0]) / (e[1][0] - e[0][0])
    cross = lambda o, a, c: (a[0] - o[0]) * (c[1] - o[1]) - (a[1] - o[1]) * (c[0] - o[0])
    done, open_ = [], []  # piece: [bottom chain left->right, top chain left->right]
    for x0, x1 in zip(xs, xs[1:]):
        if x1 - x0 < 1e-6:
            continue
        xm = (x0 + x1) / 2
        act = sorted((e for e in edges if e[0][0] < xm < e[1][0]), key=lambda e: y_at(e, xm))
        nxt = []
        for lo, hi in zip(act[0::2], act[1::2]):
            b0, b1, t0, t1 = (x0, y_at(lo, x0)), (x1, y_at(lo, x1)), (x0, y_at(hi, x0)), (x1, y_at(hi, x1))
            if t0[1] - b0[1] < 1e-6 and t1[1] - b1[1] < 1e-6:
                continue
            for p in open_:
                bot, top = p
                if abs(bot[-1][0] - x0) < 1e-7 and abs(bot[-1][1] - b0[1]) < 1e-7 and abs(top[-1][1] - t0[1]) < 1e-7 \
                        and (len(bot) < 2 or cross(bot[-2], bot[-1], b1) >= -1e-9) and (len(top) < 2 or cross(top[-2], top[-1], t1) <= 1e-9):
                    open_.remove(p)
                    bot.append(b1)
                    top.append(t1)
                    nxt.append(p)
                    break
            else:
                nxt.append([[b0, b1], [t0, t1]])
        done += open_
        open_ = nxt
    done += open_
    out = []
    for bot, top in done:
        poly = bot + top[::-1]
        clean = [q for i, q in enumerate(poly) if abs(cross(poly[i - 1], q, poly[(i + 1) % len(poly)])) > 1e-10]
        if len(clean) < 3 or area(clean) < 1e-4:
            continue
        cx, cy = sum(q[0] for q in clean) / len(clean), sum(q[1] for q in clean) / len(clean)
        shrunk = []
        for q in clean:
            d = ((q[0] - cx) ** 2 + (q[1] - cy) ** 2) ** 0.5
            k = max(0.0, 1 - SHRINK / d) if d > 0 else 0.0
            shrunk.append([round(cx + (q[0] - cx) * k, 5), round(cy + (q[1] - cy) * k, 5)])
        out.append(shrunk)
    return out


filler = pcbnew.ZONE_FILLER(b)
filler.Fill(b.Zones())
max_err = b.GetDesignSettings().m_MaxError
zones = []
skipped = []
for z in b.Zones():
    if z.GetIsRuleArea() or not z.GetNetname():
        continue
    ls = list(z.GetLayerSet().CuStack())
    if len(ls) != 1 or b.GetStandardLayerName(ls[0]) not in m['layers']:
        continue
    lid, layer, net = ls[0], b.GetStandardLayerName(ls[0]), z.GetNetname()
    # a plane net connects through its plane: its islands are optional copper that wires may cut
    if net in planes:
        skipped.append((net, layer, 'plane net island'))
        continue
    o = z.Outline()
    a = o.Area() / 1e12
    npts = sum(o.Outline(i).PointCount() for i in range(o.OutlineCount()))
    if a > 25.0 or npts > 12:
        skipped.append((net, layer, 'pour (%.1f mm2, %d points)' % (a, npts)))
        continue
    own = [p for p in b.GetPads() if p.GetNetname() == net and p.IsOnLayer(lid)]
    fill = z.GetFilledPolysList(lid)
    if fill.OutlineCount() == 0:
        skipped.append((net, layer, 'fills empty'))
        continue
    clearance = MM(z.GetLocalClearance() or 0)
    for k in range(fill.OutlineCount()):
        isl = pcbnew.SHAPE_POLY_SET()
        isl.AddOutline(fill.Outline(k))
        for j in range(fill.HoleCount(k)):
            isl.AddHole(fill.Hole(k, j))
        reached = []
        for p in own:
            ps = pcbnew.SHAPE_POLY_SET()
            p.TransformShapeToPolygon(ps, lid, pcbnew.FromMM(PAD_GAP), max_err, pcbnew.ERROR_OUTSIDE)
            hit = pcbnew.SHAPE_POLY_SET(isl)
            hit.BooleanIntersection(ps)
            if hit.OutlineCount() and hit.Area() > 1e3:
                reached.append((p, ps))
        if not reached:
            skipped.append((net, layer, 'island %d reaches no pad' % k))
            continue
        for _, ps in reached:
            isl.BooleanSubtract(ps)
        rings = []
        for i in range(isl.OutlineCount()):
            rings.append(ring(isl.Outline(i)))
            rings += [ring(isl.Hole(i, j)) for j in range(isl.HoleCount(i))]
        pieces = convex_pieces(rings)
        if not pieces:
            continue
        zones.append(dict(net=net, layer=layer, island=len(zones), clearance=clearance, area=round(isl.Area() / 1e12, 4),
                          pads=['%s.%s' % (p.GetParentFootprint().GetReference(), p.GetNumber()) for p, _ in reached], pieces=pieces))
m['zones'] = zones
(case / 'model.json').write_text(json.dumps(m, indent=2))
print(json.dumps(dict(zones=[(z['net'], z['layer'], z['pads'], len(z['pieces'])) for z in zones], skipped=skipped)))
