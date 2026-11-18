import unittest
from types import SimpleNamespace

from gridroute import necks

CU = ['F.Cu', 'In1.Cu', 'B.Cu']


def pad(ref, num, x, y, net, layers='F', w=1.0, h=1.0, drill=0.0):
    return dict(ref=ref, num=num, x=x, y=y, w=w, h=h, shape='rect', layers=layers, net=net, npth=False, drill=drill)


def board(pads, zones=(), tracks=(), vias=()):
    parts = {}
    for p in pads:
        parts.setdefault(p['ref'], []).append(p)
    return SimpleNamespace(parts={r: SimpleNamespace(placed=True, pads=lambda ps=ps: list(ps)) for r, ps in parts.items()},
                           zones=list(zones), tracks=list(tracks), vias=list(vias))


def via(net, x, y):
    return dict(net=net, x=x, y=y, d=0.5, drill=0.3)


class NeckTests(unittest.TestCase):
    def width(self, bd, reach=0.0, **kw):
        r = necks.analyse(bd, 'V', ['C1'], CU, pitch=0.05, reach=reach, **kw)
        return {t['ref']: t for t in r['pads']}

    def test_pour_with_a_waist(self):
        # two 6 mm squares on In1 joined by a 0.6 mm waist; source and load drop in through vias
        waist = [(0, 0), (6, 0), (6, 2.7), (9, 2.7), (9, 0), (15, 0), (15, 6), (9, 6), (9, 3.3), (6, 3.3), (6, 6), (0, 6)]
        bd = board([pad('C1', '1', 2, 3, 'V'), pad('U1', '1', 13, 3, 'V')],
                   zones=[dict(net='V', layer='In1.Cu', pts=waist, priority=1)], vias=[via('V', 2, 3), via('V', 13, 3)])
        t = self.width(bd, via_width=9)['U1']
        self.assertEqual(t['width'], 0.6)
        self.assertEqual(t['neck'][2], 'In1.Cu')
        self.assertTrue(5.5 < t['neck'][0] < 9.5 and 2.5 < t['neck'][1] < 3.5)

    def test_vias_of_another_net_pinch_a_strip(self):
        # a 3 mm strip with two foreign vias side by side across it: 0.5 mm of copper left between their antipads
        strip = [(0, 0), (20, 0), (20, 3), (0, 3)]
        z = [dict(net='V', layer='In1.Cu', pts=strip, priority=1)]
        ends = [pad('C1', '1', 1.5, 1.5, 'V'), pad('U1', '1', 18.5, 1.5, 'V')]
        own = [via('V', 1.5, 1.5), via('V', 18.5, 1.5)]
        clear = self.width(board(ends, z, vias=own), via_width=9, reach=1.0)['U1']      # reach: past the pad itself
        self.assertEqual(clear['width'], 3.0)
        pinched = self.width(board(ends, z, vias=own + [via('G', 10, 0.8), via('G', 10, 2.2)]), via_width=9, reach=1.0)['U1']
        self.assertEqual(pinched['width'], 0.5)
        self.assertTrue(9 < pinched['neck'][0] < 11)

    def test_a_via_counts_as_its_barrel(self):
        strip = [(0, 0), (20, 0), (20, 3), (0, 3)]
        bd = board([pad('C1', '1', 1.5, 1.5, 'V'), pad('U1', '1', 18.5, 1.5, 'V')],
                   [dict(net='V', layer='In1.Cu', pts=strip, priority=1)], vias=[via('V', 1.5, 1.5), via('V', 18.5, 1.5)])
        t = self.width(bd)['U1']
        self.assertEqual(t['width'], 0.8)              # pi x 0.3 mm = 0.94: the 0.8 step
        self.assertTrue(t['via'])

    def test_track_width_and_open_net(self):
        bd = board([pad('C1', '1', 0, 0, 'V'), pad('U1', '1', 10, 0, 'V'), pad('U2', '1', 10, 5, 'V')],
                   tracks=[dict(net='V', layer='F.Cu', width=0.4, pts=[(0, 0), (10, 0)])])
        t = self.width(bd)
        self.assertEqual(t['U1']['width'], 0.4)
        self.assertEqual(t['U2']['width'], 0.0)        # not connected
        self.assertIn('open', necks.report(necks.analyse(bd, 'V', ['C1'], CU, reach=0.0)))


if __name__ == '__main__':
    unittest.main()
