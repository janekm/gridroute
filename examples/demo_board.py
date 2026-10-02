"""A small two-layer board routed with gridroute.board: three parts, footprints defined inline, a GND plane net.

    python3 examples/demo_board.py [out.json]

Prints what was routed and writes the layout JSON (tracks, vias, placement) that a board writer would turn into a
KiCad .kicad_pcb.
"""
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
from gridroute import board  # noqa: E402

# 2-layer rules: one net class, GND carried by a plane on the bottom (its pads get via drops)
board.configure(layers=['F.Cu', 'B.Cu'], planes={'GND': 'B.Cu'}, pitch=0.05,
                classes={'Default': (0.2, 0.2, 0.6, 0.3)}, cache='off')


def chip(n_side, pitch, pad=(0.6, 1.6), span=5.0):
    """A dual-row SMD footprint: n_side pads per side, pin 1 bottom left, numbered counter-clockwise."""
    pads = []
    for k in range(n_side):
        x = (k - (n_side - 1) / 2) * pitch
        pads.append(dict(num=str(k + 1), x=x, y=span / 2, w=pad[0], h=pad[1], rot=0, shape='rect', layers='F',
                         npth=False, drill=0.0, rr=0.0))
        pads.append(dict(num=str(2 * n_side - k), x=x, y=-span / 2, w=pad[0], h=pad[1], rot=0, shape='rect',
                         layers='F', npth=False, drill=0.0, rr=0.0))
    w = n_side * pitch / 2 + 0.5
    return dict(pads=pads, courtyard=[-w, -span / 2 - 1.2, w, span / 2 + 1.2])


footprints = {
    'demo:SOIC8': chip(4, 1.27),
    'demo:R0603': dict(pads=[dict(num='1', x=-0.8, y=0, w=0.9, h=0.95, rot=0, shape='roundrect', layers='F',
                                  npth=False, drill=0.0, rr=0.2),
                             dict(num='2', x=0.8, y=0, w=0.9, h=0.95, rot=0, shape='roundrect', layers='F',
                                  npth=False, drill=0.0, rr=0.2)],
                       courtyard=[-1.5, -0.75, 1.5, 0.75]),
}
comps = {'U1': dict(footprint='demo:SOIC8'), 'U2': dict(footprint='demo:SOIC8'), 'R1': dict(footprint='demo:R0603'),
         'R2': dict(footprint='demo:R0603')}
pin_net = {('U1', '1'): 'SDA', ('U2', '8'): 'SDA', ('R1', '1'): 'SDA',
           ('U1', '2'): 'SCL', ('U2', '7'): 'SCL', ('R2', '1'): 'SCL',
           ('R1', '2'): 'VCC', ('R2', '2'): 'VCC', ('U1', '8'): 'VCC', ('U2', '1'): 'VCC',
           ('U1', '4'): 'GND', ('U2', '4'): 'GND', ('U1', '3'): 'EN', ('U2', '3'): 'EN'}

t = time.perf_counter()
bd = board.Board(40.0, 30.0, comps, pin_net, footprints)
bd.place('U1', 10.0, 15.0)
bd.place('U2', 30.0, 15.0, 180)
bd.place('R1', 20.0, 6.0)
bd.place('R2', 20.0, 24.0, 90)
bd.fanout_all(['U1', 'U2'])                       # GND pads: a short stub and a via to the plane
for net in ('SDA', 'SCL', 'EN', 'VCC'):
    ok = bd.connect(net, margin=4.0)
    print('%-4s %s' % (net, 'routed' if ok else 'FAILED'))
print('clearance check:', bd.check() or 'clean')
out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'demo_layout.json')
bd.save(out, [40.0, 30.0])
print('%d tracks, %d vias in %.2f s (native kernels: %s, GPU: %s) -> %s' % (
    len(bd.tracks), len(bd.vias), time.perf_counter() - t, board._GR is not None,
    board._GR.device_name() if board._GR is not None and board._GR.metal_available() else 'none', out))
