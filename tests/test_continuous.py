import unittest
from gridroute.continuous import Shape,distance,check
from test_router import channel

def capsule(a,b,r=.1):return Shape('N',1,(a,b),r,False,(),'track',('track',0))
class ContinuousTests(unittest.TestCase):
 def test_diagonal_intersection_and_gap(self):
  self.assertEqual(distance(capsule((0,0),(1,1)),capsule((0,1),(1,0))),0)
  self.assertAlmostEqual(distance(capsule((0,0),(1,0)),capsule((0,.25),(1,.25))),.05)
 def test_polygon_hole_has_no_copper(self):
  p=Shape('N',1,((0,0),(4,0),(4,4),(0,4)),0,True,(((1,1),(3,1),(3,3),(1,3)),),'pad',('pad','P','1',0))
  self.assertAlmostEqual(distance(p,capsule((2,2),(2,2))),.9)
  self.assertEqual(distance(p,capsule((.5,.5),(.5,.5))),0)
 def test_subgrid_gap_is_not_connected(self):
  bd=channel(.05);bd.add_track('N','F.Cu',[(.8,2.525),(3.96,2.525)],.15)
  self.assertEqual(check(bd)['missing_nets'],['N'])
  bd.add_track('N','F.Cu',[(3.96,2.525),(4.2,2.525)],.15)
  self.assertEqual(check(bd)['missing_nets'],[])
 def test_physical_clearance_between_tracks(self):
  bd=channel(.05);bd.tracks += [dict(net='A',layer='F.Cu',pts=[(1,1),(4,1)],width=.15),dict(net='B',layer='F.Cu',pts=[(1,1.2),(4,1.2)],width=.15)]
  findings=check(bd)['violations']
  self.assertTrue(any(v['kind']=='clearance' and set(v['nets'])=={'A','B'} for v in findings))
 def test_orphan_copper_is_a_missing_connection(self):
  bd=channel(.05);bd.add_track('N','F.Cu',[(.8,2.525),(4.2,2.525)],.15)
  bd.add_track('N','F.Cu',[(1,1),(1.5,1)],.15)
  result=check(bd)
  self.assertEqual(len(result['pad_components']['N']),1)
  self.assertEqual(result['copper_components']['N'],2)
  self.assertEqual(result['missing_nets'],['N'])
 def test_repair_uses_continuous_feedback_and_preserves_source(self):
  from gridroute.continuous import repair_connectivity
  from gridroute.router import Options,RoutingController
  bd=channel(.05);bd.add_track('N','F.Cu',[(.8,2.525),(3.96,2.525)],.15)
  original=RoutingController._state_key(bd)
  result,report,attempts=repair_connectivity(bd,channel,Options(fanout=False),pitches=(.025,))
  self.assertEqual(report['missing_nets'],[])
  self.assertEqual(report['violations'],[])
  self.assertTrue(attempts[0]['accepted'])
  self.assertEqual(RoutingController._state_key(bd),original)
 def test_failed_repair_is_rolled_back(self):
  from gridroute.continuous import repair_connectivity
  from gridroute.router import Options
  bd=channel(.05)
  result,report,attempts=repair_connectivity(bd,channel,Options(fanout=False,recovery_passes=0),pitches=(.05,))
  self.assertIs(result,bd);self.assertEqual(report['missing_nets'],['N'])
  self.assertFalse(attempts[0]['accepted'])
 def test_rejected_fine_candidate_restores_returned_board_configuration(self):
  from gridroute.continuous import repair_connectivity
  from gridroute.router import Options
  from gridroute import board
  bd=channel(.05)
  result,report,attempts=repair_connectivity(bd,channel,Options(fanout=False,recovery_passes=0,max_expansions=1),pitches=(.025,))
  self.assertIs(result,bd)
  self.assertEqual(board.G,bd.config['pitch'])
  self.assertTrue(bd.unrouted())
 def test_rotated_slot_distance_is_rigid_motion_invariant(self):
  from gridroute.board import rot
  for angle in [0,30,90,137]:
   a=Shape(None,3,(rot(-1,0,angle),rot(1,0,angle)),.2,False,(),'hole',('pad','H','1',0))
   b=capsule(rot(-1,.6,angle),rot(1,.6,angle),.1)
   self.assertAlmostEqual(distance(a,b),.3)
 def test_hole_clearance_is_independent_of_copper_net_clearance(self):
  bd=channel(.05)
  bd.vias=[dict(net='N',x=1.,y=1.,d=.4,drill=.3),dict(net='N',x=1.45,y=1.,d=.4,drill=.3)]
  bd.config['hole_to_hole']=.2
  result=check(bd)
  self.assertTrue(any(v['kind']=='hole_to_hole' and abs(v['actual']-.15)<1e-8 for v in result['violations']))
if __name__=='__main__':unittest.main()
