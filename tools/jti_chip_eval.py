"""Johanson 2450AT42E0100001E evaluation board (2450AT42E0100001CE1, datasheet p. 3) as Gerbers for fieldsolve,
and the chip footprint for other boards.

The eval board is 50 x 30 mm, 2-layer, assumed 0.8 mm FR4 (the part is the <= 1.0 mm substrate variant). Top
copper is ground except a 6.26 x 3 mm zone at the centre of the long edge, where the chip sits; the bottom is
solid ground, with ground vias round the edge and the zone. In the zone (frame: mm from the zone's top-left
corner, y into the board): GND pads 4 and 3 at y 0.3-0.9, tied to the top ground by 0.63 mm strips at the zone's
sides; the feed pad 1 and the floating pad 2 at y 1.9-2.5. The feed line and matching network are left out:
the stand-in is fed at pad 1.

Ports for the calibrated stand-in: FEED from pad 1 to the ground below the zone, CHIP from pad 4 to pad 1
(the chip's body, terminated in co-simulation by a series R-L-C).

    python3 jti_chip_eval.py OUTDIR      writes OUTDIR/ (Gerbers, drill) for the eval board
"""
import os, sys

ZW, ZH = 6.26, 3.0
PADS = {  # name: (x0, y0, x1, y1, net) in the zone frame
    '4': (0.63, 0.3, 2.13, 0.9, 'GND'),
    '3': (4.13, 0.3, 5.63, 0.9, 'GND'),
    '1': (0.63, 1.9, 2.13, 2.5, 'CHIP_FEED'),
    '2': (4.13, 1.9, 5.63, 2.5, 'CHIP_NC'),
}
STRIPS = [(0.0, 0.3, 0.63, 0.9), (5.63, 0.3, ZW, 0.9)]  # pads 4/3 to the top ground


def zone_shapes(zx, zy):
    """Pads and tie strips of a chip zone whose top-left corner is at (zx, zy), KiCad frame."""
    pads = {n: (zx + a, zy + b, zx + c, zy + d, net) for n, (a, b, c, d, net) in PADS.items()}
    strips = [(zx + a, zy + b, zx + c, zy + d) for a, b, c, d in STRIPS]
    return pads, strips


def ports(zx, zy):
    """FEED (pad 1 to the ground below the zone) and CHIP (pad 4 to pad 1) in KiCad coordinates."""
    xc = zx + 1.38
    return [
        {'name': 'FEED', 'layer': 'F.Cu', 'pos': [xc, zy + 2.2], 'neg': [xc, zy + ZH + 0.3], 'width': 0.6, 'r': 50},
        {'name': 'CHIP', 'layer': 'F.Cu', 'pos': [xc, zy + 0.6], 'neg': [xc, zy + 2.2], 'width': 0.6, 'r': 50},
    ]


def _num(v):
    return f'{round(v * 1e6):d}'


def rect_region(x0, y0, x1, y1):
    pts = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    s = 'G36*\n' + f'X{_num(pts[0][0])}Y{_num(-pts[0][1])}D02*\n'
    s += ''.join(f'X{_num(x)}Y{_num(-y)}D01*\n' for x, y in pts[1:])
    return s + 'G37*\n'


def gerber(function, objects):
    s = f'%TF.FileFunction,{function}*%\n%FSLAX46Y46*%\n%MOMM*%\n%LPD*%\nG01*\n'
    for net, x0, y0, x1, y1, dark in objects:
        s += '%LPD*%\n' if dark else '%LPC*%\n'
        if net:
            s += f'%TO.N,{net}*%\n'
        s += rect_region(x0, y0, x1, y1)
        if net:
            s += '%TD*%\n'
    return s + 'M02*\n'


def write_eval(out, w=50.0, h=30.0, t=0.8):
    os.makedirs(out, exist_ok=True)
    zx, zy = (w - ZW) / 2, 0.0
    pads, strips = zone_shapes(zx, zy)
    top = [('GND', 0, 0, w, h, True), (None, zx, zy, zx + ZW, zy + ZH, False)]
    top += [('GND', *s, True) for s in strips]
    top += [(net, x0, y0, x1, y1, True) for (x0, y0, x1, y1, net) in pads.values()]
    open(os.path.join(out, 'eval-F_Cu.gtl'), 'w').write(gerber('Copper,L1,Top', top))
    open(os.path.join(out, 'eval-B_Cu.gbl'), 'w').write(gerber('Copper,L2,Bot', [('GND', 0, 0, w, h, True)]))
    edge = '%TF.FileFunction,Profile,NP*%\n%FSLAX46Y46*%\n%MOMM*%\n%ADD10C,0.100000*%\nD10*\n'
    for i, (x, y) in enumerate([(0, 0), (w, 0), (w, h), (0, h), (0, 0)]):
        edge += f'X{_num(x)}Y{_num(-y)}{"D02" if i == 0 else "D01"}*\n'
    open(os.path.join(out, 'eval-Edge_Cuts.gm1'), 'w').write(edge + 'M02*\n')
    # ground vias: round the edge every 2.5 mm, 0.8 mm in, and round the zone 0.7 mm out
    vias = []
    for k in range(int(w / 2.5) + 1):
        x = min(max(0.8 + k * 2.5, 0.8), w - 0.8)
        for y in (0.8, h - 0.8):
            if not (zx - 1.2 < x < zx + ZW + 1.2 and y < 2):
                vias.append((x, y))
    for k in range(int(h / 2.5) + 1):
        y = min(max(0.8 + k * 2.5, 0.8), h - 0.8)
        vias += [(0.8, y), (w - 0.8, y)]
    for x in (zx - 0.7, zx + ZW + 0.7):
        vias += [(x, 0.8), (x, 2.2), (x, ZH + 0.7)]
    vias += [(zx + 1.0 + k * 1.5, ZH + 0.7) for k in range(4)]
    drl = 'M48\n; #@! TF.FileFunction,Plated,1,2,PTH\nFMAT,2\nMETRIC\nT1C0.300\n%\nG90\nG05\nT1\n' + ''.join(f'X{x:.3f}Y{-y:.3f}\n' for x, y in sorted(set(vias))) + 'M30\n'
    open(os.path.join(out, 'eval-PTH.drl'), 'w').write(drl)
    return zx, zy


def chip_standin(zx, zy, z_top, eps=10.0, tand=0.002, h=1.4):
    """A physical stand-in for the 2450AT42E's LTCC element: a shorted patch inside the 5.0 x 2.0 x 1.5 mm body.
    A plate h above the board's top copper (z_top) spans the body; a post from the feed pad 1 and a post from the
    GND pad 4 (both at the same end) feed and short it; a LOAD port runs from pad 3 up to the plate's far end,
    where co-simulation puts a series C (the resonance) and R (the chip's internal loss). Returns (bodies, LOAD
    port); the FEED port is ports(zx, zy)[0]."""
    x1, x3 = zx + 1.38, zx + 4.88
    y_feed, y_gnd = zy + 2.2, zy + 0.6
    zp = z_top + h
    bodies = [
        {'name': 'chip body', 'box': [zx + 0.63, zy + 0.35, z_top, zx + 5.63, zy + 2.45, z_top + 1.5], 'eps': eps, 'tand': tand},
        {'name': 'chip plate', 'box': [x1 - 0.3, y_gnd - 0.3, zp, x3 + 0.3, y_feed + 0.3, zp], 'metal': True},
        {'name': 'chip feed post', 'cylinder': [x1, y_feed, 0.2, z_top, zp], 'metal': True},
        {'name': 'chip short post', 'cylinder': [x1, y_gnd, 0.2, z_top, zp], 'metal': True},
    ]
    load = {'name': 'LOAD', 'p0': [x3, y_gnd, z_top], 'p1': [x3, y_gnd, zp], 'r': 50}
    return bodies, load


if __name__ == '__main__':
    print(write_eval(sys.argv[1]))
