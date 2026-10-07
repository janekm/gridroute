"""Board-level geometry invariants; run with PYTHONPATH=python python -m unittest discover -s tests."""
import unittest
import numpy as np
from gridroute import board


class PadSeedTests(unittest.TestCase):
    def setUp(self):
        board.configure(layers=['F.Cu','In1.Cu','B.Cu'],cache='off',spec=0)

    def make_board(self, shape, layers, x, y, w, h, rr=0):
        pad=dict(num='1',x=0,y=0,w=w,h=h,rot=0,shape=shape,layers=layers,npth=False,drill=.3 if layers=='all' else 0,rr=rr)
        bd=board.Board(5,5,{'P1':{'footprint':'f'}},{('P1','1'):'N'},{'f':{'pads':[pad],'courtyard':[0,0,0,0]}})
        bd.place('P1',x,y)
        return bd,bd.parts['P1'].pad('1')

    def test_compact_seeds_preserve_full_copper_island(self):
        # Off-grid edges and minimum-width pads exercise raster phase, not just
        # round-number centres. Compare complete masks, not only pad counts.
        for shape in ['rect','roundrect','circle','oval']:
            for layers in ['F','B','all']:
                for phase in [0,.013,.024,.026,.049]:
                    with self.subTest(shape=shape,layers=layers,phase=phase):
                        bd,p=self.make_board(shape,layers,2+phase,2-phase,.25,.9,.08 if shape=='roundrect' else 0)
                        win=(0,0,bd.ny-1,bd.nx-1)
                        pw,mask=bd._pad_mask(p,0)
                        ii,jj=np.nonzero(mask)
                        exhaustive=[(L,int(i+pw[0]),int(j+pw[1])) for L in board.pad_layers(p) for i,j in zip(ii,jj)]
                        compact=bd._pad_seeds(p,win)
                        self.assertEqual(len(compact),len(board.pad_layers(p)))
                        np.testing.assert_array_equal(bd._island(bd.net_id['N'],compact,win),bd._island(bd.net_id['N'],exhaustive,win))

    def test_clipped_window_retains_visible_pad_cells(self):
        bd,p=self.make_board('rect','F',2,2,1,1)
        win=bd._win(1.5,1.5,1.8,2.5)
        seeds=bd._pad_seeds(p,win)
        self.assertGreater(len(seeds),1)
        self.assertTrue(bd._island(bd.net_id['N'],seeds,win).any())


class RecoveryInvariantTests(unittest.TestCase):
    def setUp(self):
        board.configure(layers=['F.Cu','B.Cu'],cache='off',spec=0,pitch=.05)

    def make(self):
        pad=dict(num='1',x=0,y=0,w=.4,h=.4,rot=0,shape='circle',layers='F',npth=False,drill=0)
        bd=board.Board(5,5,{r:{'footprint':'f'} for r in ['A','B']},{('A','1'):'N',('B','1'):'N'},
                       {'f':{'pads':[pad],'courtyard':[0,0,0,0]}})
        bd.place('A',1,2);bd.place('B',4,2)
        return bd

    def test_replacement_repaints_copper_and_logs(self):
        bd=self.make();bd.add_track('N','F.Cu',[(1,2),(4,2)],.15)
        tracks=list(bd.tracks);self.assertEqual(bd.unrouted(),[])
        bd.replace_copper([],[]);self.assertEqual(len(bd.unrouted()),1)
        bd.replace_copper(tracks,[]);self.assertEqual(bd.unrouted(),[])
        self.assertEqual(len([op for op in bd.ops if op[3] is not None]),1)

    def test_native_unit_precision_never_rounds_width_below_minimum(self):
        bd=self.make();bd.add_track('N','F.Cu',[(1.000026,2),(3.999974,2)],.37592)
        self.assertEqual(bd.tracks[0]['width'],.37592)
        self.assertEqual(bd.tracks[0]['pts'][0],[1.000026,2])
        tracks=list(bd.tracks);before=bd.core.copy();bd.replace_copper(tracks,[])
        np.testing.assert_array_equal(before,bd.core)

    def test_safe_via_star_stays_inside_existing_copper(self):
        board.configure(layers=['F.Cu','B.Cu'],cache='off',spec=0,pitch=.025,
                        classes={'Default':(.25,.2,.6,.3)})
        bd=self.make();bd.add_track('N','F.Cu',[(1,2),(2.24,2)],.25)
        bd.add_track('N','F.Cu',[(2.76,2),(4,2)],.25);bd.add_via('N',2.5,2,d=.6,drill=.3)
        stars=bd.via_replacement_tracks(bd.vias[0]);self.assertEqual(len(stars),2)
        # Sample copper added by the proposed stars against the exact union of
        # the old tracks and via, independently of the routing raster.
        for x in np.linspace(2.1,2.9,161):
            for y in np.linspace(1.85,2.15,61):
                added=any(board._seg_dist(x,y,*s['pts'])<=s['width']/2 for s in stars)
                old=(np.hypot(x-2.5,y-2)<=.3 or any(board._seg_dist(x,y,*t['pts'])<=t['width']/2+1e-12 for t in bd.tracks))
                self.assertFalse(added and not old)
        bd.replace_copper(bd.tracks+stars,[]);self.assertEqual(bd.unrouted(),[])

    def test_single_layer_via_can_be_required_bridge(self):
        bd=self.make();bd.add_track('N','F.Cu',[(1,2),(2.3,2)],.15)
        bd.add_track('N','F.Cu',[(2.7,2),(4,2)],.15);bd.add_via('N',2.5,2,d=.6,drill=.3)
        self.assertEqual(bd.unrouted(),[])
        bd.prune_vias()
        self.assertEqual(bd.unrouted(),[])
        self.assertEqual(len(bd.vias),1)

    def test_tangent_tracks_cannot_justify_removing_via(self):
        bd=self.make();bd.add_track('N','F.Cu',[(1,2),(2.375,2)],.25)
        bd.add_track('N','F.Cu',[(2.625,2),(4,2)],.25);bd.add_via('N',2.5,2,d=.6,drill=.3)
        v=bd.vias[0];self.assertTrue(bd.via_bridges_copper(v))
        bd.prune_vias();self.assertEqual(len(bd.vias),1)
        # A real overlapping joint needs no copper from the redundant via.
        bd.replace_copper([],[]);bd.add_track('N','F.Cu',[(1,2),(2.5,2)],.25)
        bd.add_track('N','F.Cu',[(2.5,2),(4,2)],.25);bd.add_via('N',2.5,2,d=.6,drill=.3)
        self.assertFalse(bd.via_bridges_copper(bd.vias[0]));bd.prune_vias();self.assertEqual(len(bd.vias),0)

    def test_trim_stops_inside_pad_not_at_tangent_capsule(self):
        bd=self.make();bd.add_track('N','F.Cu',[(.5,2),(4.5,2)],.25)
        bd.trim_dangling();t=bd.tracks[0]
        self.assertAlmostEqual(board._pad_distance(*t['pts'][0],bd.parts['A'].pad('1')),0.)
        self.assertAlmostEqual(board._pad_distance(*t['pts'][-1],bd.parts['B'].pad('1')),0.)

    def test_hole_mask_includes_same_net_pth_and_slot_ends(self):
        pad=dict(num='1',x=0,y=0,w=1.2,h=.6,rot=0,shape='oval',layers='all',npth=False,
                 drill=.3,drill_w=.9,drill_h=.3)
        bd=board.Board(5,5,{'P':{'footprint':'f'}},{('P','1'):'N'},
                       {'f':{'pads':[pad],'courtyard':[0,0,0,0]}})
        bd.place('P',2.5,2.5)
        win=(0,0,bd.ny-1,bd.nx-1);mask=bd._near_holes(win,.15+.25)
        self.assertTrue(mask[50,65])
        self.assertFalse(mask[50,70])
        self.assertFalse(bd._near_holes(win,.4,exclude_net='N').any())
        bd.add_via('N',4,4,d=.5,drill=.3)
        self.assertTrue(bd._near_holes(win,.4)[80,80])
        bd._rip_ops([len(bd.ops)-1]);self.assertFalse(bd._near_holes(win,.4)[80,80])

    def test_component_labels_count_islands_not_pads(self):
        bd=self.make();self.assertEqual(bd.pad_components('N'),[0,1])
        bd.add_track('N','F.Cu',[(1,2),(4,2)],.15)
        self.assertEqual(bd.pad_components('N'),[0,0])

class RotatedGeometryTests(unittest.TestCase):
    def test_rotated_rectangle_and_slot_are_not_bounding_boxes(self):
        board.configure(layers=['F.Cu','B.Cu'],pitch=.025,cache='off',spec=0)
        p=dict(num='1',x=0,y=0,w=1.2,h=.3,rot=45,shape='rect',layers='all',npth=False,
               drill=.3,drill_w=.9,drill_h=.3)
        bd=board.Board(5,5,{'P':{'footprint':'f'}},{('P','1'):'N'},
                       {'f':{'pads':[p],'courtyard':[0,0,0,0]}});bd.place('P',2.5,2.5)
        pad=bd.parts['P'].pad('1');nid=bd.net_id['N']
        self.assertEqual(bd.core[0,88,112],nid)  # along the rotated pad
        self.assertEqual(bd.core[0,112,112],0)   # in its bounding box, outside copper
        full=(0,0,bd.ny-1,bd.nx-1);pw,mask=bd._pad_mask(pad,0);ii,jj=np.nonzero(mask)
        seeds=[(L,int(i+pw[0]),int(j+pw[1])) for L in board.pad_layers(pad) for i,j in zip(ii,jj)]
        np.testing.assert_array_equal(bd._island(nid,seeds,full),bd._island(nid,bd._pad_seeds(pad,full),full))
        holes=bd._near_holes(full,.4)
        self.assertTrue(holes[80,120]);self.assertFalse(holes[120,120])
        self.assertAlmostEqual(board._pad_distance(2.8,2.2,pad),0)
        self.assertGreater(board._pad_distance(2.8,2.8,pad),.25)

class ObstacleTests(unittest.TestCase):
    def setUp(self):
        board.configure(layers=['F.Cu','B.Cu'],pitch=.025,classes={'Default':(.127,.127,.3556,.254)},cache='off',spec=0)

    def test_small_annulus_preserves_hole_to_foreign_track_clearance(self):
        bd=board.Board(5,5,{}, {('unused','1'):'A',('unused','2'):'B'}, {})
        bd.add_via('A',2.5,2.5)
        win=(0,0,bd.ny-1,bd.nx-1)
        blocked=bd._blocked(0,bd.net_id['B'],win,bd._reach('B',.127/2))
        self.assertTrue(blocked[100,116])  # 0.4 mm centre distance is too close to this drill
        self.assertFalse(blocked[100,122])
        bd._rip_ops([len(bd.ops)-1]);self.assertFalse(bd._blocked(0,bd.net_id['B'],win,bd._reach('B',.127/2))[100,116])

    def test_polygon_holes_and_via_only_exclusions(self):
        bd=board.Board(5,5,{}, {('unused','1'):'N'}, {})
        outer=[(1,1),(4,1),(4,4),(1,4)];hole=[(2,2),(3,2),(3,3),(2,3)]
        bd.keepout_polygon(outer,layers=['F.Cu'],holes=[hole],tracks=False,vias=True)
        self.assertFalse((bd.occ[:,50,50]==board.KEEP).any())
        self.assertEqual(bd.no_via[50,50],1);self.assertEqual(bd.no_via[100,100],0)
        bd.add_track('N','B.Cu',[(1.25,1.25),(1.5,1.25)],.127)
        bd.replace_copper([],[])
        self.assertEqual(bd.no_via[50,50],1)

    def test_circular_stroke_keeps_interior_open_after_repaint(self):
        bd=board.Board(8,8,{}, {('unused','1'):'N'}, {})
        bd.keepout_ring(4,4,2,.2,clearance=.5)
        self.assertEqual(bd.core[0,160,240],board.KEEP)
        self.assertEqual(bd.core[0,160,160],0)
        self.assertEqual(bd.no_via[160,160],0)
        before=bd.occ.copy();bd.add_track('N','F.Cu',[(3,4),(4,4)],.127)
        bd.replace_copper([],[]);np.testing.assert_array_equal(before,bd.occ)

class LocalHoleMaskTests(unittest.TestCase):
    def test_local_windows_equal_full_distance_field(self):
        import math
        rng=np.random.default_rng(709)
        for pitch in [.05,.025]:
            board.configure(pitch=pitch,cache='off',spec=0)
            bd=board.Board(7,6,{}, {}, {})
            holes=[(float(x),float(y),float(w),float(h),float(a),n)
                   for x,y,w,h,a,n in zip(rng.uniform(-1,8,40),rng.uniform(-1,7,40),
                                        rng.uniform(.2,1,40),rng.uniform(.2,1,40),
                                        rng.choice([0,30,45,90,135],40),['A','B']*20)]
            bd._hole_pads=holes
            for win in [(0,0,bd.ny-1,bd.nx-1),(13,17,81,93)]:
                for radius in [.1,.35,.6]:
                    for exclude in [None,'A']:
                        ys,xs=bd._grid(win);expected=np.zeros((win[2]-win[0]+1,win[3]-win[1]+1),bool)
                        for x,y,w,h,angle,net in holes:
                            if exclude is not None and net==exclude:continue
                            c,s=math.cos(math.radians(angle)),math.sin(math.radians(angle))
                            u,v=(xs-x)*c-(ys-y)*s,(xs-x)*s+(ys-y)*c
                            dx=np.maximum(np.abs(u)-max(0,(w-h)/2),0)
                            dy=np.maximum(np.abs(v)-max(0,(h-w)/2),0)
                            expected|=dx*dx+dy*dy<(min(w,h)/2+radius)**2
                        np.testing.assert_array_equal(bd._near_holes(win,radius,exclude),expected)

class ClassClearanceTests(unittest.TestCase):
    def test_drill_clearance_does_not_subtract_netclass_surplus(self):
        board.configure(pitch=.025,margin=0,pad_grow=0,cache='off',spec=0,
                        classes={'Default':(.2,.1,.4,.2),'Wide':(.2,.2,.4,.2)},net_class=lambda n:'Wide' if n=='A' else 'Default')
        bd=board.Board(6,6,{}, {('A','1'):'A',('B','1'):'B'}, {})
        bd.add_track('B','F.Cu',[(2,3),(4,3)],.2)
        win=(0,0,bd.ny-1,bd.nx-1)
        blocked=bd._near_foreign_copper(bd.net_id['A'],win,.35)
        self.assertTrue(blocked[137,120]);self.assertFalse(blocked[139,120])
        self.assertFalse(bd._near_foreign_copper(bd.net_id['B'],win,.35).any())

    def test_classes_use_maximum_not_sum_of_surpluses(self):
        board.configure(pitch=.025,margin=0,pad_grow=0,hole_clear=.1,cache='off',spec=0,
                        classes={'Default':(.2,.1,.4,.2),'Wide':(.2,.2,.4,.2)},net_class=lambda n:'Wide')
        bd=board.Board(5,5,{}, {('A','1'):'A',('B','1'):'B'}, {})
        bd.add_via('B',2.5,2.5,d=.4,drill=.2)
        win=(0,0,bd.ny-1,bd.nx-1);r=bd._reach('A',.1)
        blocked=bd._blocked(0,bd.net_id['A'],win,r)
        self.assertTrue(blocked[100,119])   # 0.475 < .2 radius + .1 halfwidth + .2 clearance
        self.assertFalse(blocked[100,122])  # 0.550 > .500; old double-surplus threshold was .600
        np.testing.assert_array_equal(blocked,bd._blocked_np(0,bd.net_id['A'],win,r))

    def test_wider_class_adds_pad_grow_to_pads_only(self):
        board.configure(pitch=.025,margin=0,pad_grow=.05,cache='off',spec=0,
                        classes={'Default':(.2,.1,.4,.2),'Wide':(.2,.2,.4,.2)},net_class=lambda n:'Wide' if n=='A' else 'Default')
        pad=dict(num='1',x=0,y=0,w=.4,h=.4,rot=0,shape='rect',layers='F',npth=False,drill=0)
        bd=board.Board(6,6,{'P':{'footprint':'f'}},{('P','1'):'C',('A','1'):'A',('B','1'):'B'},
                       {'f':{'pads':[pad],'courtyard':[0,0,0,0]}});bd.place('P',3,1)
        bd.add_track('B','F.Cu',[(2,3),(4,3)],.2)
        win=(0,0,bd.ny-1,bd.nx-1);r=bd._reach('A',.1)
        blocked=bd._blocked(0,bd.net_id['A'],win,r)
        # a Default track: .1 half width + .1 half width + .2 clearance, no pad grow (a hand-drawn neck stays reachable)
        self.assertTrue(blocked[135,120]);self.assertFalse(blocked[137,120])
        # a Default pad keeps its pad grow under the wider class: .2 edge + .1 + .2 + .05
        self.assertTrue(blocked[61,120]);self.assertFalse(blocked[63,120])
        np.testing.assert_array_equal(blocked,bd._blocked_np(0,bd.net_id['A'],win,r))

    def test_layer_clearance_floor_and_foreign_class_both_apply(self):
        board.configure(layers=['F.Cu','In1.Cu','B.Cu'],pitch=.025,margin=0,pad_grow=0,hole_clear=.1,cache='off',spec=0,
                        classes={'Default':(.2,.1,.4,.2),'Wide':(.2,.2,.4,.2)},
                        net_class=lambda n:'Wide' if n=='B' else 'Default',layer_clearances={'In1.Cu':.3})
        bd=board.Board(5,5,{}, {('A','1'):'A',('B','1'):'B'}, {});bd.add_via('B',2.5,2.5,d=.4,drill=.2)
        win=(0,0,bd.ny-1,bd.nx-1)
        surface=bd._blocked(0,bd.net_id['A'],win,bd._reach('A',.1,'F.Cu'))
        inner=bd._blocked(1,bd.net_id['A'],win,bd._reach('A',.1,'In1.Cu'))
        self.assertFalse(surface[100,122]);self.assertTrue(inner[100,122])
        self.assertFalse(inner[100,126])
        np.testing.assert_array_equal(inner,bd._blocked_np(1,bd.net_id['A'],win,bd._reach('A',.1,'In1.Cu')))

class PadConstraintTests(unittest.TestCase):
    def setUp(self):
        board.configure(pitch=.025,margin=0,pad_grow=0,cache='off',spec=0,
                        classes={'Default':(.2,.2,.5,.25)})

    def make(self,pad):
        p=dict(num='1',x=0,y=0,w=.5,h=.5,rot=0,shape='rect',layers='F',npth=False,drill=0,**pad)
        bd=board.Board(6,6,{'P':{'footprint':'f'}},{('P','1'):'A',('other','1'):'B'},
                       {'f':{'pads':[p],'courtyard':[0,0,0,0]}});bd.place('P',3,3)
        return bd

    def test_local_pad_clearance_does_not_create_false_copper(self):
        bd=self.make({'clearance':1.})
        w=(0,0,bd.ny-1,bd.nx-1);blocked=bd._blocked(0,bd.net_id['B'],w,bd._reach('B',.1))
        self.assertTrue(blocked[120,168]);self.assertFalse(blocked[120,178])
        self.assertEqual(bd.core[0,120,160],0)

    def test_mask_aperture_only_obstacle_on_its_surface(self):
        bd=self.make({'mask_layers':['B.Cu'],'mask_margins':{'B.Cu':.4}})
        w=(0,0,bd.ny-1,bd.nx-1)
        self.assertTrue(bd._blocked(1,bd.net_id['B'],w,bd._reach('B',.1))[120,146])
        self.assertFalse(bd._blocked(0,bd.net_id['B'],w,bd._reach('B',.1))[120,146])
        self.assertEqual(bd.core[1,120,120],0)

    def test_polygon_pad_with_hole_and_offset(self):
        p=dict(num='1',x=0,y=0,w=2,h=2,rot=0,shape='polygon',layers='F',npth=False,drill=0,
               polygons=[dict(points=[(-1,-1),(1,-1),(1,1),(-1,1)],holes=[[(-.3,-.3),(.3,-.3),(.3,.3),(-.3,.3)]])])
        bd=board.Board(6,6,{'P':{'footprint':'f'}},{('P','1'):'A'}, {'f':{'pads':[p],'courtyard':[0,0,0,0]}});bd.place('P',3,3)
        pad=bd.parts['P'].pad('1');self.assertEqual(bd.core[0,120,120],0);self.assertEqual(bd.core[0,120,152],bd.net_id['A'])
        self.assertAlmostEqual(board._pad_distance(3,3,pad),.3)
        self.assertAlmostEqual(board._pad_distance(4.2,3,pad),.2)
        self.assertGreater(len(bd._pad_seeds(pad,(0,0,bd.ny-1,bd.nx-1))),1)

    def test_polygon_pad_local_offset_is_translated(self):
        p=dict(num='1',x=.7,y=.4,w=.4,h=.4,rot=0,shape='polygon',layers='F',npth=False,drill=0,
               polygons=[dict(points=[(-.2,-.2),(.2,-.2),(.2,.2),(-.2,.2)],holes=[])])
        bd=board.Board(6,6,{'P':{'footprint':'f'}},{('P','1'):'A'}, {'f':{'pads':[p],'courtyard':[0,0,0,0]}});bd.place('P',3,3)
        pad=bd.parts['P'].pad('1')
        self.assertAlmostEqual(board._pad_distance(3.7,3.4,pad),0.)
        self.assertGreater(board._pad_distance(3,3,pad),.5)

if __name__=='__main__':unittest.main()
