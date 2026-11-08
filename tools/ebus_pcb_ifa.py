"""A printed inverted-F antenna in the eBUS_IF inline board's free top-right area (fieldsolve specs).

The free area, found from the A.4 layout (no parts, tracks or signal vias): a strip along the top edge (KiCad
y < 2.0) from x 29.5 to 46.5, and the top-right corner (x > 46.5, y < 7.0), clear of the H3 screw boss at
(45, 4.3). The keep-out is cleared on all four copper layers. The arm (0.5 mm) runs along the top edge at
y 0.5-1.0 from the short at x 30.0 to x 47, steps down round the rounded corner and runs down the right edge
(x 48.5-49.0) for `d` mm. The feed stub sits `s` mm right of the short; the port spans its gap to the ground
edge. The board, decoupling, pin loads, enclosure and screws are as in ebus_fpc_specs.inst.

    python3 ebus_pcb_ifa.py NAME D S [ff] [dry] [w=ARM_WIDTH] [g=KEEPOUT_DEPTH]
"""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ebus_fpc_specs as e

LAYERS = ('F.Cu', 'In1.Cu', 'In2.Cu', 'B.Cu')
KEEPOUT = [[29.5, -0.5, 46.5, 2.0], [46.5, -0.5, 50.5, 7.0]]
X_SHORT, Y_GND = 30.0, 2.0


def ifa(d, s, w=0.5, yg=None):
    Y_GND = yg if yg is not None else globals()['Y_GND']
    y0, y1 = 0.5, 0.5 + w
    arm = [
        [X_SHORT, y0, 47.0, y1],                 # along the top edge
        [46.5, y0, 47.5, 2.3],                   # step down inside the corner arc
        [46.5, 1.8, 49.0, 2.3],                  # across to the right edge
        [48.5, 1.8, 49.0, 1.8 + max(d, 0.5)],    # down the right edge
    ]
    short = [X_SHORT, y0, X_SHORT + 0.5, Y_GND + 0.1]
    xf = X_SHORT + 0.5 + s
    feed = [xf, y0, xf + 0.5, Y_GND - 0.45]
    port = {'name': 'ANT', 'layer': 'F.Cu', 'pos': [xf + 0.25, Y_GND - 0.6], 'neg': [xf + 0.25, Y_GND + 0.2], 'width': 0.5, 'r': 50}
    return arm + [short, feed], port


def spec(name, d, s, ff=False, w=0.5, depth=None):
    """depth: the keep-out's depth along the top edge (default 2.0 mm, the free area); deeper needs D8, R10-R14,
    U6 and C2 moved out of the strip (a layout change; their copper is simply cleared here)."""
    yg = depth if depth is not None else Y_GND
    sp = e.inst(name, 0, ff, cable=False)
    sp['gerbers'] = 'gbr'
    sp['cell'] = 0.05
    sp['layers'] = []
    sp['bodies'] = [b for b in sp['bodies'] if b['name'] in ('base', 'cover') or b['name'].startswith('screw')]
    keep = [[29.5, -0.5, 46.5, yg], [46.5, -0.5, 50.5, max(7.0, yg)]]
    ed = [{'layer': L, 'op': 'clear', 'rect': r} for L in LAYERS for r in keep]
    rects, port = ifa(d, s, w, yg)
    ed += [{'layer': 'F.Cu', 'op': 'add', 'net': 'ANT', 'rect': r} for r in rects]
    sp['edits'] = ed
    sp['ports'] = [port]
    sp['roi'] = [29.0, -0.5, 50.5, 7.5]
    sp['open_nets_ok'] = sp['open_nets_ok'] + ['/MCU and radio/RF_UFL']
    if not ff:
        sp.pop('farfield', None)
    return sp


if __name__ == '__main__':
    name, d, s = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    w = next((float(a[2:]) for a in sys.argv[4:] if a.startswith('w=')), 0.5)
    g = next((float(a[2:]) for a in sys.argv[4:] if a.startswith('g=')), None)
    sp = spec(name, d, s, 'ff' in sys.argv[4:], w, g)
    if 'dry' in sys.argv[4:]:
        sp['max_steps'] = 20
        sp['scene'] = True
        sp.pop('farfield', None)
    json.dump(sp, open(name + '.json', 'w'), indent=1)
    print('wrote', name + '.json')
