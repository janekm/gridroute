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
import numpy as np
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
    seed: int = 0
    negotiate: bool = False
    negotiate_rounds: int = 30
    soft_cost_mm: float = .5
    history_cost_mm: float = .25
    max_rips: int = 8
    outward_drops: bool = False
    escape_reserve_mm: float = 0.
    outer_layer_cost: float = 1.   # per-mm cost multiplier on the outer (component) layers: >1 sends long runs inside
    pad_entry: str = 'any'          # 'axial': leave/enter SMD pads straight through their middle (Board.pad_entry)
    pad_entry_halo: float = .1      # axial: straight run beyond the pad edge (mm, at least one grid cell)
    escape_plan: bool = False
    escape_hold: bool = False        # reserve dense pins' escapes (corridor, planned via site) until their net connects
    escape_local_mm: float = 6.
    escape_buses: bool = True
    escape_dogbones: bool = True
    escape_refs: tuple = ()          # restrict the plan to these packages (original refs); () = every dense one
    neck_escapes: bool = True
    neck_nets_first: bool = True
    max_neck_escape_mm: float = 1.5
    fine_pitch_mm: float = .65
    fix_clearances: bool = True
    join_components: bool = True
    polish_via_cost_mm: float = 0.
    polish_seconds: float = 20.
    plane_drop_radii: tuple = (2., 2.9)  # within unrouted()'s 3 mm drop window


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
        bd.pad_entry=self.options.pad_entry;bd.pad_entry_halo=self.options.pad_entry_halo
        if source is not None:
            bd.replace_copper([t for t in source.tracks if only_net is None or t['net']==only_net],
                              [v for v in source.vias if only_net is None or v['net']==only_net])
            for kind,args,_,_ in source.reservations():bd._op(kind,args,None)
        return bd

    def _event(self, stage, **kw):
        self.events.append(dict(stage=stage,elapsed_seconds=time.perf_counter()-self.started,**kw))

    def _time_left(self):
        return time.perf_counter()-self.started < self.options.deadline_seconds

    def _connect(self, bd, net, layers=None):
        if self.options.neck_escapes:self._neck_net(bd,net)
        kw={}
        if self.options.outer_layer_cost!=1.:
            ls=bd.config['layers'];kw['layer_cost']={ls[0]:self.options.outer_layer_cost,ls[-1]:self.options.outer_layer_cost}
        ok=bd.connect(net,layers=layers or bd.config['layers'],strict=True,max_exp=self.options.max_expansions,**kw)
        # Two components: connect already tried that pair. Refinement boards are too costly to search twice.
        # not inside a refinement retry (too costly there); a retained fine grid is the normal board, though
        if (not ok and self.options.join_components and not getattr(self,'_refining',False)
                and len(set(bd.pad_components(net)))>=3):
            ok=self._join_components(bd,net)
        if ok and self.options.escape_hold and bd.reservations([net]):bd.release_reservations([net])
        return ok

    def _join_components(self,bd,net):
        """Board.connect grows one tree from the first pad, so a boxed-in first pad fails the whole net. Join the
        other copper components to the largest one instead; each pair is tried once, within a bounded budget."""
        pads=bd.pads_of(net)
        if len(pads)<2:return True
        groups={}
        for p,l in zip(pads,bd.pad_components(net,pads)):groups.setdefault(l,[]).append(p)
        rest=sorted(groups.values(),key=lambda g:(-len(g),g[0]['ref'],g[0]['num']))
        attempts=4*len(rest)
        while len(rest)>1 and attempts>0 and self._time_left():
            main=rest.pop(0);left=[]
            for g in rest:
                if attempts<=0 or not self._time_left():left.append(g);continue
                attempts-=1
                if bd._route(net,g,main,layers=bd.config['layers'],max_exp=self.options.max_expansions):main=main+g
                else:left.append(g)
            if not left:return True
            # the leftovers may still join each other (e.g. when the main group itself was boxed in)
            rest=sorted(left,key=lambda g:(-len(g),g[0]['ref'],g[0]['num']))
            if len(main)<sum(len(g) for g in rest):rest.append(main)
        return len(set(bd.pad_components(net,pads)))<2

    def _neck_net(self,bd,net):
        """Escape stubs for the net's neck-down pads that have none yet (placed against the current copper)."""
        if not bd.__dict__.get('necks'):return   # only boards with neck-down rules (Board.add_neck)
        pads=[p for p in bd.pads_of(net) if p.get('escape_width') and p['layers']!='all']
        if not pads or net in bd.config['planes'] or bd._net_complete(net):return
        starts={tuple(t['pts'][0]) for t in bd.tracks if t['net']==net}
        for p in pads:
            if (round(p['x'],6),round(p['y'],6)) not in starts:self._neck_escape(bd,p)

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
        planes=bd.config['planes']
        for net in bd.nets:
            ps=bd.pads_of(net)
            if net in planes:
                deficit+=sum(1 for p in ps if p['layers']!='all' and not bd.plane_reached(p))
                continue
            if len(ps)<2 or net.startswith('unconnected-'):
                continue
            labels=bd.pad_components(net)
            deficit+=max(0,len(set(labels))-1)
        return deficit

    def _order(self,bd):
        seed=self.options.seed
        def jitter(n):
            # deterministic per (seed, net): a portfolio explores orders, not run-to-run noise
            return 1.+int.from_bytes(hashlib.blake2b(f'{seed}:{n}'.encode(),digest_size=4).digest(),'big')/2**32 if seed else 1.
        def key(n):
            ps=bd.pads_of(n)
            span=math.hypot(max(p['x'] for p in ps)-min(p['x'] for p in ps),
                            max(p['y'] for p in ps)-min(p['y'] for p in ps))
            return (len(ps)>8,span*jitter(n),n)
        planes=bd.config['planes']
        return sorted([n for n in bd.nets if len(bd.pads_of(n))>=2 and not n.startswith('unconnected-') and n not in planes],key=key)

    def _refine_net(self, bd, net):
        self._refining=True
        try:
            return self._refine_net_at(bd, net)
        finally:
            self._refining=False

    def _refine_net_at(self, bd, net):
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

    def _candidates(self,bd,net,action=None):
        kind=getattr(action,'key',None);action=action or (lambda b:self._connect(b,net))
        own = repr(([t for t in bd.tracks if t['net']==net],[v for v in bd.vias if v['net']==net]))
        failure_key = (net,bd.pitch,own,kind)
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
            if self.options.late_fanout and net not in bd.config['planes']:self._fanout_net(static,net)
            ok=action(static)
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

    def _repair(self,bd,net,action=None):
        plane=net in bd.config['planes']
        for candidate in self._candidates(bd,net,action):
            for slack in self.options.blocker_slacks:
                if not self._time_left():break
                gone=self._blockers(bd,net,candidate,slack)
                affected=sorted({bd.ops[i][3]['net'] for i in gone})
                if not gone or len(gone)>self.options.max_blocker_objects or len(affected)>self.options.max_blocker_nets:
                    continue
                old_tracks,old_vias=list(bd.tracks),list(bd.vias);before=self._score(bd)
                bd._rip_ops(gone)
                if plane:first=action(bd)
                else:
                    if self.options.late_fanout:self._fanout_net(bd,net)
                    first=self._connect(bd,net)
                for n in affected:
                    if n in bd.config['planes']:self._drop_net(bd,n)
                    else:self._connect(bd,n)
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

    def _drop(self,bd,p):
        """Connect an SMD pad of a plane net to its plane with a short stub and a via.

        Larger radii are tried only when the short drop fails. A via in the pad is
        the last resort, and only for pads that can hold the whole via (exposed pads)."""
        if p['layers']=='all' or bd.plane_reached(p):return True
        width=p.get('escape_width') or bd.width(p['net'])
        # Like a dog-bone: leave the package outwards so the via does not
        # sit in the escape channels of the neighbouring pins.
        outward=self._outward_dirs(bd,p) if self.options.outward_drops else None
        axial=[(1.,0.),(-1.,0.),(0.,1.),(0.,-1.)] if self.options.pad_entry=='axial' else None
        for dirs in ([outward] if outward else [])+([axial] if axial else [])+[None]:
            for r in self.options.plane_drop_radii:
                if bd.fanout(p['ref'],p['num'],width=width,max_r=r,dirs=dirs):return True
        if min(p['w'],p['h'])>=bd.via_size(p['net'])[0]+.1:
            if bd.fanout(p['ref'],p['num'],width=width,max_r=min(p['w'],p['h'])/2,in_pad=True):return True
        return self._tie_to_dropped(bd,p)

    def _tie_to_dropped(self,bd,p,reach=6.):
        """No room for a via at the pad: route it to the nearest pads of its net that already reach the plane
        (their islands include the via), as a designer ties a crowded pin to a neighbour's via."""
        targets=sorted((q for q in bd.pads_of(p['net']) if (q['ref'],q['num'])!=(p['ref'],p['num'])
                        and math.hypot(q['x']-p['x'],q['y']-p['y'])<=reach
                        and (q['layers']=='all' or bd.plane_reached(q))),
                       key=lambda q:math.hypot(q['x']-p['x'],q['y']-p['y']))[:4]
        return bool(targets) and bd._route(p['net'],[p],targets,layers=bd.config['layers'],
                                           max_exp=self.options.max_expansions)

    @staticmethod
    def _outward_dirs(bd,p,min_cos=.35):
        """Unit directions pointing away from the pad's package centre (None for a central pad)."""
        centres=bd.__dict__.get('_package_centres')
        if centres is None:
            groups={}
            for part in bd.parts.values():
                if part.placed:
                    for q in part.pads():groups.setdefault(q['original_ref'],[]).append((q['x'],q['y']))
            centres=bd._package_centres={k:(sum(x for x,_ in v)/len(v),sum(y for _,y in v)/len(v),len(v))
                                         for k,v in groups.items()}
        cx,cy,n=centres.get(p['original_ref'],(p['x'],p['y'],1))
        vx,vy=p['x']-cx,p['y']-cy;d=math.hypot(vx,vy)
        if n<2 or d<max(p['w'],p['h'])/2:return None
        dirs=[(math.cos(a*math.pi/4),math.sin(a*math.pi/4)) for a in range(8)]
        return [u for u in dirs if (u[0]*vx+u[1]*vy)/d>=min_cos] or None

    def _reserve_escapes(self,bd):
        """Reserve a short outward corridor at every connected pin of fine-pitch packages for the initial pass,
        so earlier nets' vias and tracks do not seal later pins in (released before recovery)."""
        from .escape import dense_packages,package_centre
        count=0
        # pins whose escape is already planned (held site, neck stub or dog-bone copper starting at the pad)
        planned={tuple(op[1][2]) for op in bd.reservations()}|{tuple(t['pts'][0]) for t in bd.tracks}
        G=bd.pitch
        for ref,(ps,pitch) in dense_packages(bd,self.options.fine_pitch_mm).items():
            cx,cy=package_centre(bd,ref)
            for q in ps:
                net=q['net']
                if not net or net.startswith('unconnected-') or len(bd.pads_of(net))<2 and net not in bd.config['planes']:
                    continue
                if (round(q['x'],6),round(q['y'],6)) in planned or (q['x'],q['y']) in planned:
                    continue
                vx,vy=q['x']-cx,q['y']-cy
                if q['w']>q['h']*1.2 or (q['h']<=q['w']*1.2 and abs(vx)>=abs(vy)):
                    d=(math.copysign(1,vx),0.);half=q['w']/2;short=q['h']
                    if abs(vx)<half:continue
                else:
                    d=(0.,math.copysign(1,vy));half=q['h']/2;short=q['w']
                    if abs(vy)<half:continue
                reach=half+self.options.escape_reserve_mm
                r=min(q.get('escape_width') or bd.width(net),short)/2
                L=geometry.pad_layers(q)[0];layer=geometry.LAYERS[L]
                # clip where it would come within clearance of anything already there (other reservations too)
                win=bd._win(q['x']-reach-G,q['y']-reach-G,q['x']+reach+G,q['y']+reach+G)
                ok=bd.track_ok(net,L,2*r,win);k=0
                while k*G/2<reach:
                    x,y=q['x']+d[0]*(k+1)*G/2,q['y']+d[1]*(k+1)*G/2
                    if not ok[int(round(y/G))-win[0],int(round(x/G))-win[1]]:break
                    k+=1
                if k*G/2<=half:continue
                bd.reserve(net,layer,(q['x'],q['y']),(q['x']+d[0]*k*G/2,q['y']+d[1]*k*G/2),r);count+=1
        self._event('reserve_escapes',corridors=count)

    def _neck_escape(self,bd,p):
        """A straight escape stub at the pad's permitted neck-down width (p['escape_width']) from the pad centre
        to the nearest grid point where the net's class width fits, outwards first. It starts on the pad, so it
        touches the package courtyard as neck-down rules require."""
        net=p['net'];ew=p.get('escape_width')
        if not ew or not net or net in bd.config['planes'] or p['layers']=='all':return False
        L=geometry.pad_layers(p)[0];layer=geometry.LAYERS[L];full=bd.width(net,layer)
        if ew>=full:return False
        G=bd.pitch;R=self.options.max_neck_escape_mm+max(p['w'],p['h'])/2
        win=bd._win(p['x']-R-G,p['y']-R-G,p['x']+R+G,p['y']+R+G);i0,j0=win[0],win[1];nid=bd.net_id[net]
        thin=bd._blocked(L,nid,win,bd._reach(net,ew/2,layer))
        fat=bd._blocked(L,nid,win,bd._reach(net,full/2,layer))
        outward=self._outward_dirs(bd,p,min_cos=.7) or []
        dirs=outward+[d for d in self._outward_dirs(bd,p,min_cos=.3) or [] if d not in outward]
        def cell(x,y):return int(round(y/G))-i0,int(round(x/G))-j0
        for dx,dy in dirs:
            step=G/2;n=int(R/step)
            for k in range(1,n+1):
                x,y=p['x']+dx*k*step,p['y']+dy*k*step
                i,j=cell(x,y)
                if not(0<=i<thin.shape[0] and 0<=j<thin.shape[1]) or thin[i,j]:break
                if fat[i,j]:continue
                ex,ey=(j+j0)*G,(i+i0)*G
                m=max(int(math.hypot(ex-p['x'],ey-p['y'])/(G/2)),1)
                if any(thin[cell(p['x']+(ex-p['x'])*t/m,p['y']+(ey-p['y'])*t/m)] for t in range(m+1)):break
                bd.add_track(net,layer,[(p['x'],p['y']),(ex,ey)],ew)
                return True
        return False

    def _neck_escapes(self,bd):
        count=0
        for part in bd.parts.values():
            if not part.placed:continue
            for q in part.pads():
                if q.get('escape_width') and q['net'] and len(bd.pads_of(q['net']))>=2 and self._neck_escape(bd,q):count+=1
        self._event('neck_escapes',stubs=count)

    def _drop_net(self,bd,net):
        ok=True
        for p in sorted(bd.pads_of(net),key=lambda p:(min(p['w'],p['h']),p['ref'],p['num'])):
            if not self._time_left():break
            if not self._drop(bd,p):ok=False
        return ok

    def _plane_action(self,net,p):
        def act(b):return self._drop(b,b.parts[p['ref']].pad(p['num']))
        act.key=('drop',p['ref'],p['num'])
        return act

    def _recover_plane(self,bd,net):
        for p in bd.pads_of(net):
            if not self._time_left():break
            if p['layers']=='all' or bd.plane_reached(p):continue
            if not self._drop(bd,p):self._repair(bd,net,self._plane_action(net,p))

    def _soft_hook(self,bd,history):
        """Cost of states that only other nets' routed copper blocks (passable in a negotiated search)."""
        unit=self.options.soft_cost_mm/bd.pitch;hunit=self.options.history_cost_mm/bd.pitch
        def soft(static,net,nid,win,Ls,half,blk):
            i0,j0,i1,j1=win
            cost=np.zeros(blk.shape,dtype=np.float32)
            hist=history[:,i0:i1+1,j0:j1+1]
            self._activate(bd)
            try:
                for L in Ls:
                    wh=half if half is not None else bd.width(net,geometry.LAYERS[L])/2
                    hard=bd._blocked(L,nid,win,bd._reach(net,wh,geometry.LAYERS[L]))
                    contested=hard&~blk[L]
                    cost[L][contested]=unit*(1+hist[L][contested])
            finally:
                self._activate(static)
            cost+=hunit*hist
            return cost
        return soft

    def _negotiated_candidate(self,bd,net,history):
        static=self._fresh(bd.pitch,bd,only_net=net)
        nt,nv=len(static.tracks),len(static.vias)
        static.soft_cost=self._soft_hook(bd,history);static.contested=[]
        try:
            ok=self._connect(static,net)
        finally:
            self._activate(bd)
        if not ok:return None,[]
        return static.tracks[nt:]+static.vias[nv:],static.contested

    def _negotiate(self,bd):
        """Negotiated rip-up and reroute (PathFinder-style) for nets that still fail.

        A failed net searches the static scene, where other nets' routed copper is
        passable at a cost that grows with a per-cell history of contention. The
        objects its path crosses are ripped, the net is routed, and the victims are
        rerouted or queued. Temporary regressions are allowed; the best state is kept."""
        shape=(geometry.NL,bd.ny,bd.nx)
        history=np.zeros(shape,dtype=np.float32)
        best=(self._score(bd),list(bd.tracks),list(bd.vias))
        rips={};fails={};dead=set();planes=bd.config['planes']
        for round_ in range(self.options.negotiate_rounds):
            missing=sorted({p[0] for p in bd.unrouted()}-dead,key=lambda n:(-fails.get(n,0),n))
            if not missing or not self._time_left():break
            for net in missing:
                if not self._time_left():break
                if net in planes:
                    self._recover_plane(bd,net);continue
                if self._connect(bd,net):continue
                fails[net]=fails.get(net,0)+1
                candidate,contested=self._negotiated_candidate(bd,net,history)
                if candidate is None:
                    dead.add(net);self._event('negotiate_dead',net=net);continue
                for L,i,j in contested:history[L,i,j]+=1
                gone=self._blockers(bd,net,candidate,self.options.blocker_slacks[0])
                affected=sorted({bd.ops[i][3]['net'] for i in gone})
                if any(rips.get(n,0)>=self.options.max_rips for n in affected):
                    self._event('negotiate_skip',net=net,affected=affected);continue
                old_tracks,old_vias=list(bd.tracks),list(bd.vias)
                bd._rip_ops(gone)
                ok=self._connect(bd,net)
                if not ok:
                    bd.replace_copper(old_tracks,old_vias)
                    self._event('negotiate',net=net,removed_objects=len(gone),affected=affected,connected=False)
                    continue
                for n in affected:
                    rips[n]=rips.get(n,0)+1
                    if n in planes:self._drop_net(bd,n)
                    else:self._connect(bd,n)
                score=self._score(bd)
                self._event('negotiate',net=net,removed_objects=len(gone),affected=affected,connected=True,
                            deficit=score,round=round_)
                if score<best[0]:best=(score,list(bd.tracks),list(bd.vias))
        if self._score(bd)>best[0]:bd.replace_copper(best[1],best[2])
        self._event('negotiated',deficit=best[0],dead=sorted(dead))
        return bd

    def _fix_clearances(self,bd,rounds=3):
        """Remove raster near-misses found by the continuous checker: rip one side of each finding (the track
        first, then the other object, then both), reroute with a wider raster margin, and keep the result only if
        findings drop without losing connectivity."""
        from .continuous import check
        deficit=lambda r:sum(max(0,c-1) for c in r['copper_components'].values())
        report=check(bd)
        for _ in range(rounds):
            if not report['violations']:break   # bounded by rounds, not the routing deadline
            pairs=[[x for x in (v['a'],v['b']) if x[0] in ('track','via')] for v in report['violations']]
            pairs=[sorted(x,key=lambda y:y[0]!='track') for x in pairs if x]
            if not pairs:break
            lookup=lambda x:(bd.tracks if x[0]=='track' else bd.vias)[x[1]]
            choices=[[x[0] for x in pairs],[x[-1] for x in pairs],[y for x in pairs for y in x]]
            old_tracks,old_vias=list(bd.tracks),list(bd.vias);accept=False
            for choice,repair in [(c,False) for c in choices]+[(choices[0],True)]:
                for grow in ((.5,.25) if not repair else (.5,)):
                    items={id(lookup(x)):lookup(x) for x in choice}
                    nets=sorted({x['net'] for x in items.values()})
                    bd._rip_ops([i for i,op in enumerate(bd.ops) if op[3] is not None and id(op[3]) in items])
                    margin=geometry.MARGIN;geometry.MARGIN=margin+bd.pitch*grow
                    try:
                        for n in nets:
                            if n in bd.config['planes']:self._drop_net(bd,n)
                            elif not self._connect(bd,n) and repair:self._repair(bd,n)
                    finally:
                        geometry.MARGIN=margin
                    checked=check(bd)
                    accept=len(checked['violations'])<len(report['violations']) and deficit(checked)<=deficit(report)
                    self._event('fix_clearances',nets=nets,margin_grow=grow,repair=repair,before=len(report['violations']),
                                after=len(checked['violations']),accepted=accept)
                    if accept:break
                    bd.replace_copper(old_tracks,old_vias)
                if accept:break
            if not accept:
                # Last resort: an illegal connection is worse than a missing one. Remove the offending tracks
                # and let the net stay open (cleanup trims the remains).
                items={id(lookup(x)):lookup(x) for x in choices[0] if x[0]=='track'}
                if items:
                    bd._rip_ops([i for i,op in enumerate(bd.ops) if op[3] is not None and id(op[3]) in items])
                    checked=check(bd)
                    accept=not checked['violations'] and deficit(checked)<=deficit(report)+len(items)
                    self._event('fix_clearances',nets=sorted({x['net'] for x in items.values()}),dropped=True,
                                before=len(report['violations']),after=len(checked['violations']),accepted=accept)
                    if not accept:bd.replace_copper(old_tracks,old_vias)
                if not accept:break
            report=checked
        return bd

    def _polish(self,bd):
        """Human-like cleanup of complete nets: reroute each against the finished board with a via cost and keep
        it only if it is complete and cheaper (Board.relax), within a time budget of its own."""
        started=time.perf_counter();nets=[n for n in self._order(bd) if bd._net_complete(n)]
        vias,length=len(bd.vias),sum(math.dist(a,b) for t in bd.tracks for a,b in zip(t['pts'],t['pts'][1:]))
        saved=[]
        for n in sorted(nets,key=lambda n:-sum(1 for v in bd.vias if v['net']==n)):
            if time.perf_counter()-started>self.options.polish_seconds:break
            if not any(v['net']==n for v in bd.vias):continue
            saved+=bd.relax(nets=[n],passes=1,via_cost=self.options.polish_via_cost_mm,layers=bd.config['layers'])
        self._event('polish',nets=len(saved),vias_before=vias,vias_after=len(bd.vias),
                    length_before=round(length,2),
                    length_after=round(sum(math.dist(a,b) for t in bd.tracks for a,b in zip(t['pts'],t['pts'][1:])),2))
        return bd

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
                if net in bd.config['planes']:
                    self._recover_plane(bd,net);continue
                if not self._connect(bd,net):self._repair(bd,net)
        if self.options.cleanup:self._cleanup(bd)
        self._event('finished',deficit=self._score(bd),missing=bd.unrouted())
        return bd,self.events

    def run(self):
        self._begin()
        # The ordinary coarse pass may finish without rebuilding any scene.
        # Allocate a reusable template only when a retry actually needs one.
        bd=self._fresh(self.options.pitch,cache_template=False)
        packages={}
        if self.options.escape_hold:
            # Every dense pin's escape, planned together before anything else is routed: neck-down power
            # stubs and plane dog-bones as copper, signal via sites reserved until their net connects.
            # Power nets, plane drops and other signals then route around the reservations.
            from .escape import dense_packages,plan_dogbones
            held=dense_packages(bd,self.options.fine_pitch_mm)
            if self.options.escape_refs:held={k:v for k,v in held.items() if k in self.options.escape_refs}
            placed,unplaced=plan_dogbones(bd,held,self._time_left,commit=False)
            self._event('escape_hold',sites=placed,unplaced=unplaced)
        if self.options.escape_plan:
            # buses between facing fine-pitch rows need the lateral room the reservations below would take
            from .escape import dense_packages,route_local_buses
            packages=dense_packages(bd,self.options.fine_pitch_mm)
            if self.options.escape_refs:
                packages={k:v for k,v in packages.items() if k in self.options.escape_refs}
            local=route_local_buses(bd,lambda n,ls:self._connect(bd,n,ls),packages,self.options.escape_local_mm,
                                    self._time_left) if self.options.escape_buses else []
            self._event('local_buses',packages=sorted(packages),routed=local,deficit=self._score(bd))
        if self.options.escape_reserve_mm>0:self._reserve_escapes(bd)
        if self.options.neck_escapes and self.options.neck_nets_first and bd.__dict__.get('necks'):
            # Class-width power escapes between fine pins have the fewest options: route them first.
            necked=[n for n in self._order(bd) if any(p.get('escape_width') for p in bd.pads_of(n))]
            for net in necked:
                if not self._time_left():break
                if not self._connect(bd,net):bd=self._refine_net(bd,net)
            self._event('neck_nets',nets=len(necked),deficit=self._score(bd))
        if packages and self.options.escape_dogbones:
            from .escape import plan_dogbones
            placed,unplaced=plan_dogbones(bd,packages,self._time_left)
            self._event('dogbones',placed=placed,unplaced=unplaced)
        for net in sorted(bd.config['planes']):
            self._drop_net(bd,net)
        if bd.config['planes']:self._event('plane_drops',deficit=self._score(bd))
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
        if self.options.escape_reserve_mm>0 and not self.options.escape_hold:bd.release_reservations()
        self._event('initial',deficit=self._score(bd),missing=bd.unrouted())
        if self.options.cleanup or bd.unrouted():
            if (self.options.fallback_pitches and bd.pitch>min(self.options.fallback_pitches)
                    and (self._refinement_attempted or self.options.verify_at_finest or bd.unrouted())):
                bd=self._fresh(min(self.options.fallback_pitches),bd)
            if self.options.cleanup:self._cleanup(bd)
        if self.options.negotiate and bd.unrouted():bd=self._negotiate(bd)
        # Fewest attempts first (then name), so a deadline does not starve nets late in the alphabet.
        attempts={}
        for iteration in range(self.options.recovery_passes):
            missing=sorted({p[0] for p in bd.unrouted()})
            if not missing or not self._time_left():break
            self._seen_states.add(self._state_key(bd))
            for net in sorted(missing,key=lambda n:(attempts.get(n,0),n)):
                if not self._time_left():break
                attempts[net]=attempts.get(net,0)+1
                if net in bd.config['planes']:
                    self._recover_plane(bd,net);continue
                if self.options.late_fanout:self._fanout_net(bd,net)
                if not self._connect(bd,net):self._repair(bd,net)
        if bd.reservations():self._event('released',nets=sorted({op[1][0] for op in bd.reservations()}),count=bd.release_reservations())
        if self.options.cleanup:self._cleanup(bd)
        if self.options.polish_via_cost_mm>0:bd=self._polish(bd)
        if self.options.fix_clearances:bd=self._fix_clearances(bd)
        self._event('finished',deficit=self._score(bd),missing=bd.unrouted(),relaxed_entries=getattr(bd,'relaxed_entries',0))
        return bd,self.events
