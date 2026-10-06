"""In-process continuous copper checks, independent of the routing raster and CAD.

This is a bounded subset of DRC: physical pad connectivity, copper clearance,
hole-to-copper clearance and hole-to-hole clearance. It does not evaluate arbitrary
CAD rule expressions, solder-mask webs, board edges, zones or differential pairs.
Plane nets (configure(planes=...)) count every via and plated pad as joined by the plane.
A final native CAD validation is still required. Geometry uses millimetres.
"""
from dataclasses import dataclass
from collections import defaultdict
import math,time
from .board import _seg_dist,rot

EPS=1e-6  # one native KiCad nanometre in millimetres

@dataclass
class Shape:
    net: object
    layers: int
    points: tuple
    radius: float
    closed: bool
    holes: tuple
    role: str
    parent: tuple
    clearance: float=0.

    def __post_init__(self):
        if not self.points or not math.isfinite(self.radius) or self.radius<0 or not all(math.isfinite(x) for p in self.points for x in p):
            raise ValueError('Invalid continuous geometry')
        xs,ys=zip(*self.points);r=self.radius
        self.box=(min(xs)-r,min(ys)-r,max(xs)+r,max(ys)+r)
        self.edges=tuple(zip(self.points,self.points[1:]+self.points[:1])) if self.closed else ((self.points[0],self.points[-1]),)
        for h in self.holes:self.edges+=tuple(zip(h,h[1:]+h[:1]))


def _inside(p,s):
    if not s.closed:return False
    def ring(points):
        x,y=p;value=False
        for (ax,ay),(bx,by) in zip(points,points[1:]+points[:1]):
            if (ay>y)!=(by>y) and x<(bx-ax)*(y-ay)/(by-ay)+ax:value=not value
        return value
    return ring(s.points) and not any(ring(h) for h in s.holes)


def _segment_distance(a,b,c,d):
    def orient(p,q,r):return (q[0]-p[0])*(r[1]-p[1])-(q[1]-p[1])*(r[0]-p[0])
    if orient(a,b,c)*orient(a,b,d)<0 and orient(c,d,a)*orient(c,d,b)<0:return 0.
    return min(_seg_dist(*a,c,d),_seg_dist(*b,c,d),_seg_dist(*c,a,b),_seg_dist(*d,a,b))


def distance(a,b):
    """Exact distance between capsules, rounded rectangles and polygon copper."""
    if _inside(a.points[0],b) or _inside(b.points[0],a):return 0.
    raw=min(_segment_distance(p,q,u,v) for p,q in a.edges for u,v in b.edges)
    return max(0.,raw-a.radius-b.radius)


def _pairs(shapes,reach):
    """Sweep-line broad phase; exact distances are evaluated only for nearby boxes."""
    active=[]
    for i in sorted(range(len(shapes)),key=lambda i:shapes[i].box[0]):
        a=shapes[i];x0,y0,x1,y1=a.box
        active=[j for j in active if shapes[j].box[2]+reach>=x0]
        for j in active:
            b=shapes[j]
            if a.layers & b.layers and b.box[3]+reach>=y0 and y1+reach>=b.box[1]:yield j,i
        active.append(i)


def shapes_from_board(bd):
    layers=bd.config['layers'];all_layers=(1<<len(layers))-1;shapes=[]
    def add(net,mask,points,radius,closed,holes,role,parent,clearance=0.):
        shapes.append(Shape(net,mask,tuple(map(tuple,points)),radius,closed,
                            tuple(tuple(map(tuple,h)) for h in holes),role,parent,clearance))
    def capsule(net,mask,x,y,w,h,angle,role,parent,clearance=0.):
        dx=max(0.,(w-h)/2);dy=max(0.,(h-w)/2)
        a=rot(-dx,-dy,angle);b=rot(dx,dy,angle)
        add(net,mask,[(x+a[0],y+a[1]),(x+b[0],y+b[1])],min(w,h)/2,False,(),role,parent,clearance)
    for part in bd.parts.values():
        if not part.placed:continue
        for pi,p in enumerate(part.pads()):
            parent=('pad',p['ref'],p['num'],pi);net=p['net'];mask=all_layers if p['layers']=='all' else 1<<(0 if p['layers']=='F' else len(layers)-1)
            if not p['npth']:
                if p['shape']=='polygon':
                    for poly in p['polygons']:add(net,mask,poly['points'],0.,True,poly.get('holes',()),'pad',parent,p.get('clearance',0.))
                else:
                    rr=min(p['w'],p['h'])/2 if p['shape'] in ('circle','oval') else min(p.get('rr',0.),p['w']/2,p['h']/2)
                    w,h=p['w']/2-rr,p['h']/2-rr;pts=[]
                    for x,y in [(-w,-h),(w,-h),(w,h),(-w,h)]:
                        dx,dy=rot(x,y,p.get('angle',0.));pts.append((p['x']+dx,p['y']+dy))
                    add(net,mask,pts,rr,True,(),'pad',parent,p.get('clearance',0.))
            if p['drill']:
                capsule(net,all_layers,p['x'],p['y'],p.get('drill_w',p['drill']),p.get('drill_h',p['drill']),p.get('angle',0.),'hole',parent)
    for ti,t in enumerate(bd.tracks):
        for a,b in zip(t['pts'],t['pts'][1:]):add(t['net'],1<<layers.index(t['layer']),[a,b],t['width']/2,False,(),'track',('track',ti))
    for vi,v in enumerate(bd.vias):
        parent=('via',vi);xy=[(v['x'],v['y'])]
        add(v['net'],all_layers,xy,v['d']/2,False,(),'via',parent)
        add(v['net'],all_layers,xy,v['drill']/2,False,(),'hole',parent)
    return shapes


def check(bd,clearances=True):
    """Return pad connectivity and new-copper standard clearance findings.

    Static pad/pad defects are not emitted. The findings are a subset of DRC,
    never a claim that arbitrary CAD constraints passed. Counters expose cost.
    """
    started=time.perf_counter();shapes=shapes_from_board(bd);copper=[s for s in shapes if s.role!='hole'];parents=list(range(len(copper)))
    def root(i):
        while parents[i]!=i:parents[i]=parents[parents[i]];i=parents[i]
        return i
    def union(i,j):
        i,j=root(i),root(j)
        if i!=j:parents[j]=i
    grouped=defaultdict(list)
    for i,s in enumerate(copper):grouped[s.net].append(i)
    connectivity_pairs=0
    for net,ids in grouped.items():
        if net is None:continue
        owner={}
        for i in ids:
            if copper[i].parent in owner:union(i,owner[copper[i].parent])
            else:owner[copper[i].parent]=i
        scene=[copper[i] for i in ids]
        for a,b in _pairs(scene,EPS):
            ia,ib=ids[a],ids[b]
            if root(ia)==root(ib):continue
            connectivity_pairs+=1
            if distance(scene[a],scene[b])<=EPS:union(ia,ib)
    # A plane net's vias and plated holes all reach its plane layer (not modelled as copper here).
    for net in bd.config.get('planes',{}):
        anchors=[i for i in grouped.get(net,()) if copper[i].role=='via' or
                 (copper[i].role=='pad' and copper[i].layers==(1<<len(bd.config['layers']))-1)]
        for i in anchors[1:]:union(anchors[0],i)
    components={};copper_components={};missing=[]
    for net,ids in grouped.items():
        if net is None or str(net).startswith('unconnected-'):continue
        pads=defaultdict(set)
        for i in ids:
            if copper[i].role=='pad':pads[root(i)].add(copper[i].parent[1:3])
        if not pads:continue
        components[net]=[sorted(x) for x in pads.values()]
        copper_components[net]=len({root(i) for i in ids})
        if copper_components[net]>1:missing.append(net)
    violations=[];clearance_pairs=0;cfg=bd.config
    if clearances:
        def net_clear(s):return cfg['classes'][bd.cls[s.net]][1] if s.net in bd.cls else min(r[1] for r in cfg['classes'].values())
        reach=max([cfg['hole_clear'],cfg['hole_to_hole']]+[r[1] for r in cfg['classes'].values()]+list(cfg['layer_clearances'].values())+[s.clearance for s in shapes])
        for i,j in _pairs(shapes,reach):
            a,b=shapes[i],shapes[j]
            if a.parent==b.parent:continue
            if a.parent[0]=='pad' and b.parent[0]=='pad':continue
            holes=(a.role=='hole')+(b.role=='hole')
            if holes==2:required=cfg['hole_to_hole'];kind='hole_to_hole'
            elif holes:
                if a.net==b.net and a.net is not None:continue
                required=cfg['hole_clear'];kind='hole_clearance'
            else:
                if a.net==b.net and a.net is not None:continue
                common=a.layers & b.layers
                floor=max([0.]+[cfg['layer_clearances'].get(l,0.) for k,l in enumerate(cfg['layers']) if common&(1<<k)])
                required=max(net_clear(a),net_clear(b),a.clearance,b.clearance,floor);kind='clearance'
            clearance_pairs+=1;actual=distance(a,b)
            if actual+EPS<required:violations.append(dict(kind=kind,a=a.parent,b=b.parent,required=required,actual=actual,nets=[a.net,b.net]))
    return dict(missing_nets=sorted(missing),pad_components=components,copper_components=copper_components,violations=violations,
                seconds=time.perf_counter()-started,primitives=len(shapes),connectivity_pairs=connectivity_pairs,
                clearance_pairs=clearance_pairs,coverage=['continuous_pad_connectivity','copper_clearance','hole_clearance','hole_to_hole'],
                requires_final_native_validation=True)


def repair_connectivity(source,factory,options=None,pitches=(.025,.0125),deadline_seconds=30.,max_repair_nets=4):
    """Bounded transactional repair driven only by continuous in-process checks.

    Returns (board, report, attempts). The source is never mutated. Repairs are
    limited to max_repair_nets (default four) and the shared time budget. A candidate
    is retained only if copper-component deficit improves (or a clearance fault is fixed),
    connectivity never worsens, and the candidate has no standard
    clearance findings. Unsupported constraints still require final native DRC.
    No CAD process is called. Global Board configuration makes this serial-only.
    """
    from .router import RoutingController,Options
    from dataclasses import replace
    if not math.isfinite(deadline_seconds) or deadline_seconds<=0:raise ValueError('Invalid deadline')
    if any(not math.isfinite(p) or p<=0 for p in pitches):raise ValueError('Invalid pitch')
    if not isinstance(max_repair_nets,int) or max_repair_nets<1:raise ValueError('Invalid repair net bound')
    started=time.perf_counter();best=source;report=check(best);attempts=[]
    deficit=lambda r:sum(max(0,c-1) for c in r['copper_components'].values())
    for pitch in pitches:
        if not report['missing_nets'] and not report['violations']:break
        remaining=deadline_seconds-(time.perf_counter()-started)
        if remaining<=0:break
        opts=options or Options(candidate_paths=3,plateau_budget=12,late_fanout=True,
                                recovery_passes=3,max_expansions=500000,cleanup=False)
        opts=replace(opts,pitch=pitch,fallback_pitches=(),deadline_seconds=min(opts.deadline_seconds,remaining),cleanup=False)
        nets=set(report['missing_nets'])
        nets.update(n for v in report['violations'] for n in v['nets'] if n in best.net_id)
        if len(nets)>max_repair_nets:
            attempts.append(dict(accepted=False,reason='repair_net_bound',nets=len(nets),limit=max_repair_nets));break
        candidate,events=RoutingController(factory,opts).repair(best,sorted(nets))
        checked=check(candidate)
        accepted=(deficit(checked)<=deficit(report) and not checked['violations'] and
                  (deficit(checked)<deficit(report) or bool(report['violations'])))
        attempts.append(dict(pitch=pitch,accepted=accepted,missing_nets=checked['missing_nets'],
                             pad_deficit=deficit(checked),violations=len(checked['violations']),
                             check_seconds=checked['seconds'],events=events))
        if accepted:best,report=candidate,checked
    RoutingController._activate(best)
    return best,report,attempts
