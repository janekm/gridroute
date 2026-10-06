"""Small physical fixtures for the routing policy, independent of PCBench."""
import unittest
from gridroute import board
from gridroute.router import RoutingController,Options


def channel(pitch):
    board.configure(layers=['F.Cu'],pitch=pitch,margin=pitch/2,pad_grow=pitch/2,
                    classes={'Default':(.15,.15,.45,.2)},cache='off',spec=0,search='exact')
    pad=dict(num='1',x=0,y=0,w=.3,h=.15,rot=0,shape='rect',layers='F',npth=False,drill=0)
    bd=board.Board(5,5,{r:{'footprint':'f'} for r in ['A','B']},
                   {('A','1'):'N',('B','1'):'N'},{'f':{'pads':[pad],'courtyard':[0,0,0,0]}})
    bd.place('A',.8,2.525);bd.place('B',4.2,2.525)
    bd.keepout_polygon([(1.5,0),(3.5,0),(3.5,2.275),(1.5,2.275)],layers=['F.Cu'])
    bd.keepout_polygon([(1.5,2.775),(3.5,2.775),(3.5,5),(1.5,5)],layers=['F.Cu'])
    return bd


class ControllerTests(unittest.TestCase):
    def test_static_cache_isolated_from_routed_copper_and_bounded(self):
        calls=[]
        def factory(pitch):calls.append(pitch);return channel(pitch)
        ctl=RoutingController(factory)
        first=ctl._fresh(.05);first.add_track('N','F.Cu',[(.8,2.525),(4.2,2.525)],.15)
        second=ctl._fresh(.05)
        self.assertEqual(calls,[.05]);self.assertEqual(second.tracks,[])
        self.assertTrue(second.unrouted());self.assertFalse(first.unrouted())
        ctl=RoutingController(factory,Options(static_cache_bytes=0))
        ctl._fresh(.05);ctl._fresh(.05);self.assertEqual(len(calls),3)

    def test_initial_scene_can_skip_unused_template_allocation(self):
        ctl=RoutingController(channel)
        ctl._fresh(.05,cache_template=False)
        self.assertEqual(ctl._templates,{})
        self.assertEqual(ctl._template_bytes,0)
        ctl._fresh(.025)
        self.assertIn(.025,ctl._templates)

    def test_external_repair_preserves_source_and_rebuilds_selected_net(self):
        source=channel(.05)
        source.add_track('N','F.Cu',[(.8,2.525),(1.1,2.525)],.15)
        old=RoutingController._state_key(source)
        ctl=RoutingController(channel,Options(pitch=.025,fallback_pitches=(),fanout=False,cleanup=False))
        candidate,events=ctl.repair(source,['N'])
        self.assertEqual(candidate.unrouted(),[])
        self.assertEqual(RoutingController._state_key(source),old)
        self.assertTrue(any(e['stage']=='external_feedback' for e in events))
        with self.assertRaises(ValueError):ctl.repair(source,['UNKNOWN'])

    def test_finer_fallback_opens_channel_without_changing_physical_rules(self):
        coarse=RoutingController(channel,Options(fanout=False,cleanup=False,fallback_pitches=(),recovery_passes=0))
        bd,_=coarse.run();self.assertTrue(bd.unrouted())
        ctl=RoutingController(channel,Options(fanout=False,cleanup=False,recovery_passes=0))
        fine,events=ctl.run();self.assertEqual(fine.unrouted(),[])
        self.assertTrue(any(e['stage']=='refine' and e['accepted'] for e in events))
        self.assertTrue(all(t['width']==.15 for t in fine.tracks))
        self.assertEqual(fine.clearance('N'),.15)
        first=RoutingController._state_key(fine)
        repeated,repeated_events=ctl.run()
        self.assertEqual(first,RoutingController._state_key(repeated))
        self.assertEqual(len(events),len(repeated_events))

    def test_invalid_budgets_and_pitches_fail_before_building(self):
        for kw in [dict(pitch=float('nan')),dict(fallback_pitches=(.05,)),
                   dict(candidate_paths=0),dict(max_expansions=0),dict(deadline_seconds=-1),
                   dict(blocker_slacks=()),dict(blocker_slacks=(float('nan'),)),
                   dict(recovery_passes=-1),dict(static_cache_bytes=-1),dict(max_fanout_radius=-1)]:
            with self.subTest(kw=kw),self.assertRaises(ValueError):
                RoutingController(lambda pitch:self.fail('Factory must not run'),Options(**kw))

if __name__=='__main__':unittest.main()
