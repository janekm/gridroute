"""Bounded, deterministic whole-board routing built on the Board primitives.

The board factory is called with a pitch in millimetres. It must return a fresh,
placed Board with identical physical rules and geometry at every pitch. Numerical
raster allowances may scale with pitch; physical track widths/clearances may not.
The current Board toolkit has process-global configuration: use one controller
per process, with speculation disabled, and do not access other Boards mid-run.
"""
from dataclasses import dataclass
import copy
import hashlib
import math
import time
from . import board as geometry


@dataclass(frozen=True)
class Options:
    pitch: float = .05
    fallback_pitches: tuple = (.025,)
    retain_fine_grid: bool = False
    verify_at_finest: bool = False
    late_fanout: bool = False
    priority_nets: tuple = ()
    pre_route_nets: tuple = ()
    direct_before_fanout_mm: float = 0.
    plateau_budget: int = 0
    candidate_paths: int = 1
    cache_static: bool = True
    static_cache_bytes: int = 128 * 1024 * 1024
    fanout: bool = True
    max_fanout_radius: float = 2.5
    recovery_passes: int = 2
    max_blocker_nets: int = 6
    max_blocker_objects: int = 64
    blocker_slacks: tuple = (.025, .1, .3)
    deadline_seconds: float = 30.
    max_expansions: int = 4_000_000
    cleanup: bool = True
    power_nets: tuple = ('GND', '+3V3', '+5V', 'VCC')


def _segment_distance(a, b, c, d):
    def orient(a, b, c):
        return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    if orient(a,b,c)*orient(a,b,d)<0 and orient(c,d,a)*orient(c,d,b)<0:
        return 0.
    return min(geometry._seg_dist(*a,c,d),geometry._seg_dist(*b,c,d),
               geometry._seg_dist(*c,a,b),geometry._seg_dist(*d,a,b))


def _segments(item):
    if 'pts' in item:
        return list(zip(item['pts'],item['pts'][1:]))
    p=(item['x'],item['y'])
    return [(p,p)]


class RoutingController:
    """Escalate only failed nets: refinement, then actual-obstacle local repair.

    No footprints move. Default recovery is transactional and must reduce the
    number of disconnected copper components at pads. An explicit plateau
    budget permits bounded, non-repeating equal-score moves. Results remain
    candidates until a native CAD DRC accepts them.
    """
    def __init__(self, make_board, options=None):
        self.make_board=make_board
        self.options=options or Options()
        if self.options.pitch <= 0 or any(p <= 0 or p >= self.options.pitch for p in self.options.fallback_pitches):
            raise ValueError('Fallback pitches must be positive and finer than the primary pitch')
        if not all(math.isfinite(p) for p in (self.options.pitch,*self.options.fallback_pitches)):
            raise ValueError('Grid pitches must be finite')
        if self.options.candidate_paths<1 or self.options.max_expansions<1 or self.options.plateau_budget<0:
            raise ValueError('Search limits must be positive and plateau budget nonnegative')
        if (self.options.recovery_passes<0 or self.options.max_blocker_nets<0 or
                self.options.max_blocker_objects<0 or self.options.static_cache_bytes<0):
            raise ValueError('Recovery and cache budgets must be nonnegative')
        if not self.options.blocker_slacks or any(not math.isfinite(x) or x<0 for x in self.options.blocker_slacks):
            raise ValueError('At least one finite, nonnegative blocker slack is required')
        if (not math.isfinite(self.options.max_fanout_radius) or self.options.max_fanout_radius<0 or
                not math.isfinite(self.options.direct_before_fanout_mm) or self.options.direct_before_fanout_mm<0):
            raise ValueError('Fanout and direct-routing distances must be finite and nonnegative')
        if not math.isfinite(self.options.deadline_seconds) or self.options.deadline_seconds<=0:
            raise ValueError('Routing deadline must be finite and positive')
        self.events=[]
        self._static_failures=set()
        self._seen_states=set()
        self._plateaus=0
        self._refinement_attempted=False
        self._templates={}
        self._template_bytes=0
        self.started=0.

    @staticmethod
    def _activate(bd):
        geometry.configure(**bd.config)

    def _fresh(self, pitch, source=None, only_net=None, *, cache_template=True):
        template=self._templates.get(pitch)
        if template is not None:
            bd=copy.deepcopy(template)
            self._activate(bd)
        else:
            bd=self.make_board(pitch)
            size=sum(getattr(v,'nbytes',0) for v in vars(bd).values())
            if cache_template and self.options.cache_static and self._template_bytes+size<=self.options.static_cache_bytes:
                self._templates[pitch]=copy.deepcopy(bd)
                self._template_bytes+=size
        if bd.config['spec']:
            raise ValueError('RoutingController requires spec=0 because Board configuration is process-global')
        if source is not None:
            bd.replace_copper([t for t in source.tracks if only_net is None or t['net']==only_net],
                              [v for v in source.vias if only_net is None or v['net']==only_net])
        return bd

    def _event(self, stage, **kw):
        self.events.append(dict(stage=stage,elapsed_seconds=time.perf_counter()-self.started,**kw))

    def _time_left(self):
        return time.perf_counter()-self.started < self.options.deadline_seconds

    def _connect(self, bd, net):
        return bd.connect(net,layers=bd.config['layers'],strict=True,max_exp=self.options.max_expansions)

    def _fanout_net(self,bd,net):
        for p in sorted(bd.pads_of(net),key=lambda p:(min(p['w'],p['h']),p['ref'],p['num'])):
            if not self._time_left():break
            if p['layers']!='all':
                bd.fanout(p['ref'],p['num'],width=p.get('escape_width') or bd.width(net),max_r=self.options.max_fanout_radius)

    @staticmethod
    def _score(bd):
        # Count distinct components, not missing pads: many missing pads can be
        # one disconnected island. Fall back to pad deficits for unusual shapes.
        deficit=0
        for net in bd.nets:
            ps=bd.pads_of(net)
            if len(ps)<2 or net.startswith('unconnected-'):
                continue
            labels=bd.pad_components(net)
            deficit+=max(0,len(set(labels))-1)
        return deficit

    @staticmethod
    def _order(bd):
        def key(n):
            ps=bd.pads_of(n)
            span=math.hypot(max(p['x'] for p in ps)-min(p['x'] for p in ps),
                            max(p['y'] for p in ps)-min(p['y'] for p in ps))
            return (len(ps)>8,span,n)
        return sorted([n for n in bd.nets if len(bd.pads_of(n))>=2 and not n.startswith('unconnected-')],key=key)

    def _refine_net(self, bd, net):
        for pitch in self.options.fallback_pitches:
            if pitch>=bd.pitch or not self._time_left():
                continue
            self._refinement_attempted=True
            fine=self._fresh(pitch,bd)
            nt,nv=len(fine.tracks),len(fine.vias)
            if self.options.late_fanout:self._fanout_net(fine,net)
            ok=self._connect(fine,net)
            if ok and self.options.retain_fine_grid:
                self._event('refine',net=net,pitch=pitch,accepted=True,retained_grid=True)
                return fine
            self._activate(bd)
            if ok:
                for t in fine.tracks[nt:]:bd.add_track(t['net'],t['layer'],t['pts'],t['width'])
                for v in fine.vias[nv:]:bd.add_via(v['net'],v['x'],v['y'],d=v['d'],drill=v['drill'])
            self._event('refine',net=net,pitch=pitch,accepted=ok)
            if ok:return bd
        return bd

    @staticmethod
    def _via_layers(bd,v):
        layers=set()
        for t in bd.tracks:
            if t['net']==v['net'] and any(geometry._seg_dist(v['x'],v['y'],a,b)<(v['d']+t['width'])/2-1e-3 for a,b in _segments(t)):
                layers.add(t['layer'])
        for p in bd.pads_of(v['net']):
            if geometry._pad_distance(v['x'],v['y'],p)<v['d']/2-1e-3:
                layers.update(geometry.LAYERS[L] for L in geometry.pad_layers(p))
        return layers

    def _cleanup(self,bd):
        removed=0
        for v in sorted(list(bd.vias),key=lambda v:(v['x'],v['y'],v['net'])):
            if v['net'] in bd.config['planes'] or len(self._via_layers(bd,v))>=2:
                continue
            if bd.via_bridges_copper(v):
                stars=bd.via_replacement_tracks(v)
                if stars is None:continue
                for t in stars:bd.add_track(t['net'],t['layer'],t['pts'],t['width'])
            idx=next((i for i,op in enumerate(bd.ops) if op[3] is v),None)
            if idx is None:continue
            before=bd.pad_components(v['net'])
            bd._rip_ops([idx]);k=len(bd.ops);nv=len(bd.vias)
            if bd.pad_components(v['net'])!=before:self._connect(bd,v['net'])
            if bd.pad_components(v['net'])==before and len(bd.vias)==nv:
                removed+=1
            else:
                bd._rip_ops(list(range(k,len(bd.ops))))
                bd.add_via(v['net'],v['x'],v['y'],d=v['d'],drill=v['drill'])
        old_tracks,old_vias=list(bd.tracks),list(bd.vias);before=self._score(bd)
        length=bd.trim_dangling()
        if self._score(bd)>before:
            bd.replace_copper(old_tracks,old_vias);length=0.
        self._event('cleanup',removed_vias=removed,trimmed_mm=length)

    def _blockers(self,bd,net,candidate,slack):
        def box(segments):
            points=[p for ab in segments for p in ab]
            return min(p[0] for p in points),min(p[1] for p in points),max(p[0] for p in points),max(p[1] for p in points)
        def separated(a,b,d):
            return a[2]+d<b[0] or b[2]+d<a[0] or a[3]+d<b[1] or b[3]+d<a[1]
        def segments_with_boxes(segments):
            return [(a,b,(min(a[0],b[0]),min(a[1],b[1]),max(a[0],b[0]),max(a[1],b[1]))) for a,b in segments]
        prepared=[]
        for c in candidate:
            segments=_segments(c);prepared.append((c,segments_with_boxes(segments),box(segments)))
        gone=[]
        for i,op in enumerate(bd.ops):
            item=op[3]
            if item is None or item['net']==net:continue
            item_segments=_segments(item);item_box=box(item_segments);item_segments=segments_with_boxes(item_segments)
            for c,segments,candidate_box in prepared:
                if 'pts' in c and 'pts' in item and item['layer']!=c['layer']:continue
                radius=(c.get('width',c.get('d'))+item.get('width',item.get('d')))/2
                common_layers=([c['layer']] if 'pts' in c else [item['layer']] if 'pts' in item else bd.config['layers'])
                clearance=max([bd.clearance(net),bd.clearance(item['net'])]+
                              [bd.config['layer_clearances'].get(l,0.) for l in common_layers])
                limit=radius+clearance
                if 'drill' in c:limit=max(limit,c['drill']/2+item.get('width',item.get('d'))/2+bd.config['hole_clear'])
                if 'drill' in item:limit=max(limit,item['drill']/2+c.get('width',c.get('d'))/2+bd.config['hole_clear'])
                if 'drill' in c and 'drill' in item:limit=max(limit,(c['drill']+item['drill'])/2+bd.config['hole_to_hole'])
                reach=limit+slack
                if separated(item_box,candidate_box,reach):continue
                if any(not separated(abox,ubox,reach) and _segment_distance(a,b,u,v)<reach
                       for a,b,abox in segments for u,v,ubox in item_segments):
                    gone.append(i);break
        return gone

    @staticmethod
    def _state_key(bd):
        tracks=sorted((t['net'],t['layer'],t['width'],tuple(map(tuple,t['pts']))) for t in bd.tracks)
        vias=sorted((v['net'],v['x'],v['y'],v['d'],v['drill']) for v in bd.vias)
        return hashlib.blake2b(repr((tracks,vias)).encode(),digest_size=16).digest()

    def _candidates(self,bd,net):
        own = repr(([t for t in bd.tracks if t['net']==net],[v for v in bd.vias if v['net']==net]))
        failure_key = (net,bd.pitch,own)
        if failure_key in self._static_failures:
            return []
        # Remove dynamic copper only in a temporary scene to find a legal path
        # through static geometry. Include own copper to preserve its islands.
        static=self._fresh(bd.pitch,bd,only_net=net)
        own_tracks,own_vias=list(static.tracks),list(static.vias)
        protected=set();candidates=[]
        for attempt in range(self.options.candidate_paths):
            if not self._time_left():break
            self._activate(static)
            static.replace_copper(own_tracks+[t for t in bd.tracks if t['net'] in protected],
                                  own_vias+[v for v in bd.vias if v['net'] in protected])
            nt,nv=len(static.tracks),len(static.vias)
            if self.options.late_fanout:self._fanout_net(static,net)
            ok=self._connect(static,net)
            candidate=static.tracks[nt:]+static.vias[nv:]
            self._activate(bd)
            if not ok or not candidate:
                if not ok and attempt==0:self._static_failures.add(failure_key)
                self._event('static_candidate',net=net,success=ok,attempt=attempt)
                break
            gone=self._blockers(bd,net,candidate,self.options.blocker_slacks[0])
            counts={}
            for i in gone:
                n=bd.ops[i][3]['net'];counts[n]=counts.get(n,0)+1
            candidates.append((len(counts),len(gone),candidate))
            if not counts:break
            # Preserve the most disruptive blocker in the next temporary scene.
            # This creates an alternative corridor without softening clearance.
            protected.add(max(counts,key=lambda n:(counts[n],n)))
        self._activate(bd)
        return [c for _,_,c in sorted(candidates,key=lambda c:c[:2])]

    def _repair(self,bd,net):
        for candidate in self._candidates(bd,net):
            for slack in self.options.blocker_slacks:
                if not self._time_left():break
                gone=self._blockers(bd,net,candidate,slack)
                affected=sorted({bd.ops[i][3]['net'] for i in gone})
                if not gone or len(gone)>self.options.max_blocker_objects or len(affected)>self.options.max_blocker_nets:
                    continue
                old_tracks,old_vias=list(bd.tracks),list(bd.vias);before=self._score(bd)
                bd._rip_ops(gone)
                if self.options.late_fanout:self._fanout_net(bd,net)
                first=self._connect(bd,net)
                for n in affected:self._connect(bd,n)
                after=self._score(bd);plateau=False
                accept=first and after<before
                if first and after==before and self._plateaus<self.options.plateau_budget:
                    key=self._state_key(bd)
                    if key not in self._seen_states:
                        accept=plateau=True;self._plateaus+=1
                if accept:self._seen_states.add(self._state_key(bd))
                self._event('blocker_repair',net=net,removed_objects=len(gone),affected=affected,slack=slack,
                            before=before,after=after,accepted=accept,target_connected=first,plateau=plateau)
                if accept:return True
                bd.replace_copper(old_tracks,old_vias)
        return False

    def _begin(self):
        self.events=[];self._static_failures.clear();self._seen_states.clear();self._plateaus=0
        self._refinement_attempted=False
        self._templates.clear();self._template_bytes=0
        self.started=time.perf_counter()

    def repair(self, source, nets, pitch=None):
        """Create a candidate for nets rejected by an external connectivity oracle.

        The source is not modified. Its selected nets are removed and rebuilt,
        even when the raster considers them complete. Callers must validate the
        returned candidate and retain the source when native acceptance worsens.
        This deliberately does not equate raster connectivity with CAD acceptance.
        """
        nets=list(dict.fromkeys(nets));pitch=self.options.pitch if pitch is None else pitch
        if not math.isfinite(pitch) or pitch<=0:
            raise ValueError('Repair pitch must be finite and positive')
        if set(nets)-set(source.nets):
            raise ValueError('Repair targets must be existing net names')
        self._begin();bd=self._fresh(pitch,source,cache_template=False)
        bd.replace_copper([t for t in bd.tracks if t['net'] not in nets],
                          [v for v in bd.vias if v['net'] not in nets])
        self._event('external_feedback',nets=nets,pitch=pitch)
        for net in nets:
            if not self._time_left():break
            if self.options.fanout:self._fanout_net(bd,net)
            if not self._connect(bd,net):self._repair(bd,net)
        for _ in range(self.options.recovery_passes):
            missing=sorted({p[0] for p in bd.unrouted()})
            if not missing or not self._time_left():break
            for net in missing:
                if not self._time_left():break
                if not self._connect(bd,net):self._repair(bd,net)
        if self.options.cleanup:self._cleanup(bd)
        self._event('finished',deficit=self._score(bd),missing=bd.unrouted())
        return bd,self.events

    def run(self):
        self._begin()
        # The ordinary coarse pass may finish without rebuilding any scene.
        # Allocate a reusable template only when a retry actually needs one.
        bd=self._fresh(self.options.pitch,cache_template=False)
        for net in self.options.pre_route_nets:
            if net not in bd.nets or not self._time_left():continue
            if self.options.fanout:self._fanout_net(bd,net)
            if not self._connect(bd,net):bd=self._refine_net(bd,net)
        if self.options.pre_route_nets:self._event('preroute',deficit=self._score(bd))
        if self.options.direct_before_fanout_mm>0:
            for net in self._order(bd):
                if not self._time_left():break
                ps=bd.pads_of(net)
                span=math.hypot(max(p['x'] for p in ps)-min(p['x'] for p in ps),max(p['y'] for p in ps)-min(p['y'] for p in ps))
                if len(ps)<=4 and span<=self.options.direct_before_fanout_mm:self._connect(bd,net)
            self._event('direct',deficit=self._score(bd))
        if self.options.fanout:
            nets=bd.nets
            if self.options.direct_before_fanout_mm>0 or self.options.pre_route_nets:
                nets=[n for n in nets if len(bd.pads_of(n))>=2 and not n.startswith('unconnected-') and not bd._net_complete(n)]
            ps=[p for n in nets for p in bd.pads_of(n) if p['layers']!='all']
            ps.sort(key=lambda p:(p['net'] not in self.options.power_nets,min(p['w'],p['h']),p['ref'],p['num']))
            for p in ps:
                if not self._time_left():break
                bd.fanout(p['ref'],p['num'],width=p.get('escape_width') or bd.width(p['net']),max_r=self.options.max_fanout_radius)
        order=self._order(bd)
        preferred=[n for n in self.options.priority_nets if n in order]
        order=list(dict.fromkeys(preferred+order))
        for n in order:
            if not self._time_left():break
            if not self._connect(bd,n):bd=self._refine_net(bd,n)
        self._event('initial',deficit=self._score(bd),missing=bd.unrouted())
        if self.options.cleanup or bd.unrouted():
            if (self.options.fallback_pitches and bd.pitch>min(self.options.fallback_pitches)
                    and (self._refinement_attempted or self.options.verify_at_finest or bd.unrouted())):
                bd=self._fresh(min(self.options.fallback_pitches),bd)
            if self.options.cleanup:self._cleanup(bd)
        for iteration in range(self.options.recovery_passes):
            missing=sorted({p[0] for p in bd.unrouted()})
            if not missing or not self._time_left():break
            self._seen_states.add(self._state_key(bd))
            for net in missing:
                if not self._time_left():break
                if self.options.late_fanout:self._fanout_net(bd,net)
                if not self._connect(bd,net):self._repair(bd,net)
        if self.options.cleanup:self._cleanup(bd)
        self._event('finished',deficit=self._score(bd),missing=bd.unrouted())
        return bd,self.events
