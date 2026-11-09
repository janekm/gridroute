"""The eBUS_IF inline board widened at its top edge, with a printed inverted-F antenna in the new strip (fieldsolve).

The board and the enclosure are stretched by W mm at KiCad y = CUT (just inside the top edge): everything at
y < CUT moves out by W, and the section at CUT is extruded across the gap. For the board this keeps the
original copper and adds a bare W mm strip above the old top edge (the outline's corner arcs move with it);
for the enclosure it moves the top wall out by W and lengthens what crosses y = CUT (the lid's three ribs at
x 19.4, 29.3 and 44.7 and the base's four PCB standoffs) by W. The ground pour on all four layers ends 0.3 mm
inside the old edge (y 0.3), which becomes the antenna's ground edge.

The inverted-F (F.Cu, net ANT): an arm of width AW at 0.5 mm from the new top edge, from the short at XS to
the open end at XE; the short runs from the arm down into the ground pour; the feed stub sits S mm right of
the short, over J3's centre pin (x 24.2), and the port spans its 0.45 mm gap to the pour. Stitching vias
(0.3 mm, metal posts through the board) tie the four ground layers together at the short and the feed.
Decoupling, pin loads and screws are as in ebus_fpc_specs.inst; the cable and the FPC are left out.

    python3 ebus_wide_ifa.py prep W                              writes gbr_wW/ and encl_wW/ from gbr/ and $ENC
    python3 ebus_wide_ifa.py spec NAME W XE S [AW] [ff] [dry] [cell=MM] [line]   writes NAME.json
    ... line: feed from J3's centre pin through a 50 ohm coplanar line instead of at the antenna's feed leg
"""
import json, math, os, shutil, struct, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ebus_fpc_specs as e

CUT = 1.0137        # KiCad y of the stretch plane (a value no mesh line or STL vertex sits on)
BW, BH, R = 50.0, 28.0, 4.0
Y_GND = 0.3         # ground pour edge, all layers
X_FEED = 24.2       # J3 centre pin


# ---- board outline -------------------------------------------------------------------------------------------
def outline(w, n_arc=24):
    """The stretched outline (KiCad frame, y down) as a closed polyline: the 4 mm corner arcs are cut at y = CUT
    and their outer parts moved up by w, with straight sides across the gap."""
    def arc(cx, cy, a0, a1):
        return [(cx + R * math.cos(a0 + (a1 - a0) * k / n_arc), cy + R * math.sin(a0 + (a1 - a0) * k / n_arc)) for k in range(n_arc + 1)]

    xc = R - math.sqrt(R * R - (R - CUT) ** 2)            # the arcs' x at y = CUT
    a_tl = math.atan2(CUT - R, xc - R) + 2 * math.pi      # in (pi, 1.5 pi)
    a_tr = math.atan2(CUT - R, R - xc) + 2 * math.pi      # in (1.5 pi, 2 pi)
    pts = arc(R, R, math.pi, a_tl)                        # (0, 4) -> (xc, CUT)
    pts += arc(R, R - w, a_tl, 1.5 * math.pi)             # (xc, CUT - w) -> (4, -w)
    pts += arc(BW - R, R - w, 1.5 * math.pi, a_tr)        # (46, -w) -> (50 - xc, CUT - w)
    pts += arc(BW - R, R, a_tr, 2 * math.pi)              # (50 - xc, CUT) -> (50, 4)
    pts += arc(BW - R, BH - R, 0.0, 0.5 * math.pi)        # (50, 24) -> (46, 28)
    pts += arc(R, BH - R, 0.5 * math.pi, math.pi)         # (4, 28) -> (0, 24)
    pts.append(pts[0])
    return pts


def write_outline(path, pts):
    num = lambda v: f'{round(v * 1e6):d}'
    s = '%TF.FileFunction,Profile,NP*%\n%FSLAX46Y46*%\n%MOMM*%\n%LPD*%\nG01*\n%ADD10C,0.100000*%\nD10*\n'
    for i, (x, y) in enumerate(pts):
        s += f'X{num(x)}Y{num(-y)}{"D02" if i == 0 else "D01"}*\n'
    open(path, 'w').write(s + 'M02*\n')


# ---- enclosure -------------------------------------------------------------------------------------------------
def read_stl(p):
    b = open(p, 'rb').read()
    n = struct.unpack('<I', b[80:84])[0]
    out = []
    for i in range(n):
        v = struct.unpack('<12f', b[84 + 50 * i: 84 + 50 * i + 48])
        out.append([tuple(v[3:6]), tuple(v[6:9]), tuple(v[9:12])])
    return out


def write_stl(p, tris):
    with open(p, 'wb') as f:
        f.write(b'stretched'.ljust(80, b' '))
        f.write(struct.pack('<I', len(tris)))
        for t in tris:
            f.write(struct.pack('<12fH', 0, 0, 0, *t[0], *t[1], *t[2], 0))


def stretch(tris, yc, w):
    """Cut a closed mesh at frame y = yc, move the part at y > yc by +w and extrude the section across the gap.
    Edge crossings are computed from the edge's endpoints in a fixed order, so neighbours share them exactly."""
    def cross(p, q):
        a, b = (p, q) if p < q else (q, p)
        s = (yc - a[1]) / (b[1] - a[1])
        return (a[0] + s * (b[0] - a[0]), yc, a[2] + s * (b[2] - a[2]))

    mv = lambda p: (p[0], p[1] + w, p[2])
    out = []
    for t in tris:
        side = [p[1] > yc for p in t]
        if all(side):
            out.append([mv(p) for p in t])
            continue
        if not any(side):
            out.append(list(t))
            continue
        # the lone vertex (on its own side) and the other two, keeping the winding
        k = next(i for i in range(3) if side.count(side[i]) == 1)
        a, b, c = t[k], t[(k + 1) % 3], t[(k + 2) % 3]
        pb, pc = cross(a, b), cross(a, c)
        lone = [a, pb, pc]
        pair = [[pb, b, c], [pb, c, pc]]
        if side[k]:
            out.append([mv(p) for p in lone])
            out += pair
        else:
            out.append(lone)
            out += [[mv(p) for p in q] for q in pair]
        # the gap: the section segment pb-pc extruded from yc to yc + w
        out += [[pb, pc, mv(pc)], [pb, mv(pc), mv(pb)]]
    return out


def prep(w, gbr='gbr'):
    gd = f'gbr_w{w:g}'
    os.makedirs(gd, exist_ok=True)
    for f in os.listdir(gbr):
        if f.endswith('Edge_Cuts.gm1'):
            write_outline(os.path.join(gd, f), outline(w))
        elif f.endswith(('.gtl', '.gbl', '.g1', '.g2', '.gts', '.gbs', '.drl', '.gbrjob')):
            shutil.copy(os.path.join(gbr, f), gd)
    ed = f'encl_w{w:g}'
    os.makedirs(ed, exist_ok=True)
    for f in ('asm-base.stl', 'asm-cover.stl'):
        t = stretch(read_stl(os.path.join(e.ENC, f)), -CUT, w)     # frame y_up: y = -KiCad y
        write_stl(os.path.join(ed, f), t)
    return gd, ed


# ---- antenna ---------------------------------------------------------------------------------------------------
def ifa(w, xe, s, aw=0.5):
    y0 = -w + 0.5
    xf = X_FEED - 0.25
    xs = xf - 0.5 - s
    rects = [
        [xs, y0, xe, y0 + aw],                    # arm
        [xs, y0, xs + 0.5, Y_GND + 0.3],          # short, into the pour
        [xf, y0, xf + 0.5, Y_GND - 0.45],         # feed stub
    ]
    port = {'name': 'ANT', 'layer': 'F.Cu', 'pos': [X_FEED, Y_GND - 0.6], 'neg': [X_FEED, Y_GND + 0.2], 'width': 0.5, 'r': 50}
    vias = [(xs + 0.25, Y_GND + 0.3), (xs - 0.6, Y_GND + 0.3), (xs + 1.1, Y_GND + 0.3),
            (X_FEED - 0.9, Y_GND + 0.3), (X_FEED + 0.9, Y_GND + 0.3)]
    return rects, port, vias


# The feed line from J3's centre pin (pad 23.73-24.73 x 4.71-5.77) up to the feed leg: 50 ohm grounded coplanar
# on F.Cu over In1 (0.21 mm FR4): 0.36 mm wide, 0.2 mm gaps cut in the pour from the ground edge to J3's
# clearance. J3 itself goes (its pads stay as copper); the In1 clearance under it is filled so the line and
# the port have a reference. The port is then at the pin, from F.Cu down to In1.
LINE_W, LINE_G = 0.36, 0.2
J3_PIN = (24.23, 5.24)


def feed_line(w):
    x = X_FEED
    edits = [
        {'layer': 'F.Cu', 'op': 'clear', 'rect': [x - LINE_W / 2 - LINE_G, -0.1, x + LINE_W / 2 + LINE_G, 2.7]},
        {'layer': 'F.Cu', 'op': 'add', 'rect': [x - LINE_W / 2, Y_GND - 0.7, x + LINE_W / 2, 4.9]},
        # J3's In1 clearance (x 23.27-25.13 at y 2.35-4.75, x 22.73-26.79 at y 4.89-6.4), filled up to the pour
        # round it and kept 0.2 mm off the In1 track along its top (y 2.05)
        {'layer': 'In1.Cu', 'op': 'add', 'rect': [23.2, 2.33, 25.2, 4.9]},
        {'layer': 'In1.Cu', 'op': 'add', 'rect': [22.7, 4.8, 26.8, 6.3]},
    ]
    port = {'name': 'ANT', 'at': list(J3_PIN), 'from': 'F.Cu', 'to': 'In1.Cu', 'r': 50}
    return edits, port


def spec(name, w, xe, s, aw=0.5, ff=False, line=False):
    sp = e.inst(name, 0, ff, cable=False)
    sp['gerbers'] = f'gbr_w{w:g}'
    sp['cell'] = 0.05
    sp['layers'] = []
    ed = os.path.abspath(f'encl_w{w:g}')
    sp['bodies'] = [
        {'name': 'base', 'stl': os.path.join(ed, 'asm-base.stl'), 'frame': 'y_up', **e.PCABS},
        {'name': 'cover', 'stl': os.path.join(ed, 'asm-cover.stl'), 'frame': 'y_up', **e.PCABS},
    ] + [b for b in sp['bodies'] if b['name'].startswith('screw')]
    rects, port, vias = ifa(w, xe, s, aw)
    # no net: the short joins the GND pour (copper of two nets that touch is cut apart by the net-aware raster)
    sp['edits'] = [{'layer': 'F.Cu', 'op': 'add', 'rect': r} for r in rects]
    sp['bodies'] += [{'name': f'stitch {i}', 'cylinder': [x, y, 0.15, 0.0, 0.8], 'metal': True} for i, (x, y) in enumerate(vias)]
    sp['ports'] = [port]
    x0 = min(r[0] for r in rects)
    sp['roi'] = [x0 - 2.0, -w - 0.5, xe + 2.0, 3.0]
    if line:
        ed2, sp['ports'] = feed_line(w)
        sp['ports'] = [sp['ports']]
        sp['edits'] += ed2
        sp['roi'][3] = 7.0
    sp['coarse_copper_ok'] = [[-2.0, -w - 2.0, 52.0, 30.0]]
    sp['open_nets_ok'] = sp['open_nets_ok'] + ['/MCU and radio/RF_UFL']
    if not ff:
        sp.pop('farfield', None)
    return sp


if __name__ == '__main__':
    if sys.argv[1] == 'prep':
        print(prep(float(sys.argv[2])))
    else:
        name, w, xe, s = sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
        rest = sys.argv[6:]
        aw = next((float(a) for a in rest if a.replace('.', '', 1).isdigit()), 0.5)
        sp = spec(name, w, xe, s, aw, 'ff' in rest, 'line' in rest)
        cell = next((float(a[5:]) for a in rest if a.startswith('cell=')), None)
        if cell:
            sp['cell'] = cell   # the two-mesh check: about 0.7 x 0.05
        if 'dry' in rest:
            sp['max_steps'] = 20
            sp['scene'] = True
            sp.pop('farfield', None)
        json.dump(sp, open(name + '.json', 'w'), indent=1)
        print('wrote', name + '.json')
