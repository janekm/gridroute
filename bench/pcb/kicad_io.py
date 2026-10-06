"""Extract exact placed pads with pcbnew, or write gridroute copper back to KiCad."""
import collections
import hashlib
import json
import math
import os
import pathlib
import re
import sys
import uuid

import pcbnew

ROOT = pathlib.Path(__file__).resolve().parent
MM = pcbnew.ToMM


def normalize_nil_uuids(path, source_hash):
    """Give invalid legacy object IDs deterministic identities, without moving copper.

    KiCad DRC resolves a nil footprint ID to an arbitrary footprint, producing
    unstable item descriptions and misleading baseline differences.
    """
    count=0
    def replace(match):
        nonlocal count
        count+=1
        return '(uuid "'+str(uuid.uuid5(uuid.NAMESPACE_URL,source_hash+':nil:'+str(count)))+'")'
    text=path.read_text()
    fixed=re.sub(r'\(uuid\s+"00000000-0000-0000-0000-000000000000"\)',replace,text)
    if count:path.write_text(fixed)
    return count


def project(path, rules, classes=None, net_class=None):
    width, clearance, diameter, drill = rules
    classes = classes or {'Default':rules}
    data = {
        'meta': {'filename': path.with_suffix('.kicad_pro').name, 'version': 1},
        'board': {'design_settings': {'rules': {
            'min_clearance': min(r[1] for r in classes.values()), 'min_track_width': min(r[0] for r in classes.values()),
            'min_via_diameter': min(r[2] for r in classes.values()), 'min_through_hole_diameter': min(r[3] for r in classes.values()),
            'min_via_annular_width': min((r[2]-r[3])/2 for r in classes.values()),
            'min_hole_clearance': 0.25, 'min_hole_to_hole': 0.25,
            'min_copper_edge_clearance': 0.5}}},
        'net_settings': {'classes': [{'name': 'Default', 'clearance': clearance,
            'track_width': width, 'via_diameter': diameter, 'via_drill': drill}],
            'meta': {'version': 3}, 'netclass_assignments': {},
            'netclass_patterns': [{'netclass': 'Default', 'pattern': '*'}]},
    }
    data['net_settings']['classes'] = [dict(name=n,clearance=r[1],track_width=r[0],via_diameter=r[2],via_drill=r[3]) for n,r in classes.items()]
    data['net_settings']['netclass_patterns'] = [dict(netclass=c,pattern=n) for n,c in sorted((net_class or {}).items()) if c!='Default']
    path.with_suffix('.kicad_pro').write_text(json.dumps(data, indent=2))
    if len(classes)>1:
        rules_text=['(version 1)']
        for name,r in classes.items():
            rules_text.append('(rule '+json.dumps(name+' benchmark minimum width')+' (constraint track_width (min '+str(r[0])+'mm)) (condition '+json.dumps("A.NetClass == '"+name+"'")+'))')
        path.with_suffix('.kicad_dru').write_text('\n'.join(rules_text)+'\n')


def prepare(name, source_dir):
    source = pathlib.Path(source_dir) / name / 'processed.kicad_pcb'
    out = ROOT / 'cases' / name
    out.mkdir(parents=True, exist_ok=True)
    b = pcbnew.LoadBoard(str(source))
    zones = list(b.Zones())
    rdl = json.loads((source.parent / 'final.json').read_text())
    routing_layers=[b.GetStandardLayerName(b.GetLayerID(n)) for n in rdl['layers']]
    nc = b.GetDesignSettings().m_NetSettings.GetDefaultNetclass()
    rules = [MM(nc.GetTrackWidth()), MM(nc.GetClearance()), MM(nc.GetViaDiameter()), MM(nc.GetViaDrill())]
    active_classes = [c for c in rdl['rules']['net_classes'] if c['indices']]
    available={str(n):[MM(c.GetTrackWidth()),MM(c.GetClearance()),MM(c.GetViaDiameter()),MM(c.GetViaDrill())] for n,c in b.GetAllNetClasses().items()}
    classes={'Default':rules};net_class={}
    code_name={n.GetNetCode():n.GetNetname() for n in b.GetNetsByNetcode().values()}
    for i,c in enumerate(active_classes):
        target=[c['width'],c['clearance'],c['via_diameter']]
        matches=[n for n,r in available.items() if all(abs(a-z)<1e-6 for a,z in zip(r[:3],target))]
        assert matches,('RDL class not present in board',target)
        cls_name=sorted(matches,key=lambda n:(n!='Default',n))[0];classes[cls_name]=available[cls_name]
        for code in c['indices']:net_class[code_name[code]]=cls_name
    poly = pcbnew.SHAPE_POLY_SET()
    assert b.GetBoardPolygonOutlines(poly, False), 'invalid board outline'
    assert poly.OutlineCount() == 1, 'multiple disjoint board outlines need explicit support'
    chain = poly.Outline(0)
    outline = [[MM(chain.CPoint(i).x), MM(chain.CPoint(i).y)] for i in range(chain.PointCount())]
    x0, y0 = min(p[0] for p in outline), min(p[1] for p in outline)
    size = [max(p[0] for p in outline) - x0, max(p[1] for p in outline) - y0]
    outline_holes=[]
    for j in range(poly.HoleCount(0)):
        c=poly.Hole(0,j)
        outline_holes.append([[MM(c.CPoint(i).x)-x0,MM(c.CPoint(i).y)-y0] for i in range(c.PointCount())])
    pads = []
    ignored = []
    shape_names = {pcbnew.PAD_SHAPE_RECT: 'rect', pcbnew.PAD_SHAPE_ROUNDRECT: 'roundrect',
                   pcbnew.PAD_SHAPE_CIRCLE: 'circle', pcbnew.PAD_SHAPE_OVAL: 'oval'}
    for fp_index, fp in enumerate(b.GetFootprints()):
        for pad_index, p in enumerate(fp.Pads()):
            ls = p.GetLayerSet()
            layers = [name for name in routing_layers if ls.Contains(b.GetLayerID(name))]
            if not layers:
                ignored.append([fp.GetReference(), p.GetNumber(), 'no copper layers'])
                continue
            tht = p.GetAttribute() in (pcbnew.PAD_ATTRIB_PTH, pcbnew.PAD_ATTRIB_NPTH)
            assert tht or layers in [['F.Cu'], ['B.Cu']], layers
            shape = p.GetShape(pcbnew.F_Cu)
            angle = p.GetOrientation().AsDegrees() % 180
            s = p.GetSize(pcbnew.F_Cu)
            w, h = MM(s.x), MM(s.y)
            if abs(angle - 90) < 0.01:
                w, h = h, w
            drill = p.GetDrillSize()
            dw, dh = MM(drill.x), MM(drill.y)
            if abs(angle - 90) < .01: dw, dh = dh, dw
            polygons=[]
            offset=p.GetOffset(pcbnew.F_Cu)
            if shape not in shape_names or offset.x or offset.y:
                ps=pcbnew.SHAPE_POLY_SET()
                p.TransformShapeToPolygon(ps,p.GetLayer(),0,pcbnew.FromMM(.002),pcbnew.ERROR_OUTSIDE)
                def relative_points(c):
                    return [[MM(c.CPoint(i).x-p.GetPosition().x),MM(c.CPoint(i).y-p.GetPosition().y)] for i in range(c.PointCount())]
                for i in range(ps.OutlineCount()):
                    polygons.append(dict(points=relative_points(ps.Outline(i)),holes=[relative_points(ps.Hole(i,j)) for j in range(ps.HoleCount(i))]))
                assert polygons,('unsupported empty pad geometry',shape)
            pads.append(dict(ref='pad_%d_%d' % (fp_index, pad_index), original_ref=fp.GetReference(),
                original_pin=p.GetNumber(), num='1', x=MM(p.GetPosition().x)-x0, y=MM(p.GetPosition().y)-y0,
                w=w, h=h, rot=0, angle=0 if min(abs(angle),abs(angle-90),abs(angle-180))<.01 else angle,
                shape='polygon' if polygons else shape_names[shape], polygons=polygons,layers='all' if tht else ('F' if layers==['F.Cu'] else 'B'),
                npth=p.GetAttribute()==pcbnew.PAD_ATTRIB_NPTH,
                drill=min(dw,dh), drill_w=dw, drill_h=dh, net=p.GetNetname() if p.GetNetCode() else None,
                mask_layers=[c for c,mask in [('F.Cu',pcbnew.F_Mask),('B.Cu',pcbnew.B_Mask)] if ls.Contains(mask)],
                mask_margins={c:MM(p.GetSolderMaskExpansion(b.GetLayerID(c))) for c in ['F.Cu','B.Cu']},
                clearance=MM(p.GetOwnClearance(p.GetLayer())),
                rr=MM(p.GetRoundRectCornerRadius(pcbnew.F_Cu)) if shape==pcbnew.PAD_SHAPE_ROUNDRECT else 0))
    net_counts = collections.Counter(p['net'] for p in pads if p['net'])
    copper_graphics = []
    mask_graphics = []
    edge_graphics = []
    items = list(b.GetDrawings()) + [g for f in b.GetFootprints() for g in f.GraphicalItems()]
    for g in items:
        mask_layer = {pcbnew.F_Mask:'F.Cu',pcbnew.B_Mask:'B.Cu'}.get(g.GetLayer())
        margin_layer = g.GetLayer()==pcbnew.Margin
        if not g.IsOnCopperLayer() and mask_layer is None and not margin_layer:
            continue
        target=edge_graphics if margin_layer else mask_graphics if mask_layer else copper_graphics
        if isinstance(g,pcbnew.PCB_SHAPE) and g.GetShape()==pcbnew.SHAPE_T_CIRCLE and not g.IsSolidFill():
            center=g.GetCenter();cx,cy=MM(center.x)-x0,MM(center.y)-y0
            radius,width=MM(g.GetRadius()),MM(g.GetWidth())
            def circle(r):return [[cx+r*math.cos(i*math.tau/192),cy+r*math.sin(i*math.tau/192)] for i in range(192)]
            target.append(dict(layer=mask_layer or b.GetStandardLayerName(g.GetLayer()),
                               circle=dict(x=cx,y=cy,radius=radius,width=width),
                               points=circle(radius+width/2),holes=[circle(max(0,radius-width/2))]))
            continue
        polygons = pcbnew.SHAPE_POLY_SET()
        g.TransformShapeToPolygon(polygons, g.GetLayer(), 0, pcbnew.FromMM(.002), pcbnew.ERROR_OUTSIDE)
        def points(chain):
            return [[MM(chain.CPoint(i).x)-x0, MM(chain.CPoint(i).y)-y0] for i in range(chain.PointCount())]
        for i in range(polygons.OutlineCount()):
            target.append(dict(layer=mask_layer or b.GetStandardLayerName(g.GetLayer()), points=points(polygons.Outline(i)),
                holes=[points(polygons.Hole(i,j)) for j in range(polygons.HoleCount(i))]))
    model = dict(name=name, size=size, origin=[x0,y0], outline=[[x-x0,y-y0] for x,y in outline],
                 layers=routing_layers, rules=rules, pads=pads, ignored_pads=ignored,
                 footprint_count=len(list(b.GetFootprints())), routable_nets=sum(n>=2 for n in net_counts.values()),
                 required_pad_connections=sum(n-1 for n in net_counts.values() if n>=2),
                 differential_pairs=rdl['differential_pairs'],
                 source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(), kicad_version=pcbnew.GetBuildVersion())
    model['copper_graphics'] = copper_graphics
    model['mask_graphics'] = mask_graphics
    model['edge_graphics'] = edge_graphics
    model['outline_holes'] = outline_holes
    model['classes'],model['net_class'] = classes,net_class
    model['source_min_track_width'] = MM(b.GetDesignSettings().m_TrackMinWidth)
    (out/'model.json').write_text(json.dumps(model,indent=2))
    reference = out/'reference.kicad_pcb'
    pcbnew.SaveBoard(str(reference),b)
    normalize_nil_uuids(reference,model['source_sha256'])
    project(reference,rules,classes,net_class)
    for t in list(b.GetTracks()): b.Remove(t)
    for z in zones:
        if z.GetIsRuleArea():
            raise ValueError('PCBench importer needs explicit rule-area extraction for this case')
        b.Remove(z)
    stripped = out/'unrouted.kicad_pcb'
    pcbnew.SaveBoard(str(stripped),b)
    model['normalized_nil_uuids']=normalize_nil_uuids(stripped,model['source_sha256'])
    (out/'model.json').write_text(json.dumps(model,indent=2))
    project(stripped,rules,classes,net_class)
    print(name, len(pads), model['routable_nets'], model['required_pad_connections'], flush=True)


def export(name,out):
    out=pathlib.Path(out);case=out/'input' if (out/'input/model.json').exists() else ROOT/'cases'/name
    model=json.loads((case/'model.json').read_text());layout=json.loads((out/'layout.json').read_text())
    b=pcbnew.LoadBoard(str(case/'unrouted.kicad_pcb'));zone_refs=list(b.Zones());ox,oy=model['origin']
    point=lambda p:pcbnew.VECTOR2I(pcbnew.FromMM(p[0]+ox),pcbnew.FromMM(p[1]+oy))
    nets={n.GetNetname():n.GetNetCode() for n in b.GetNetsByNetcode().values()}
    for t in layout['tracks']:
        for a,z in zip(t['pts'],t['pts'][1:]):
            if a==z:continue
            tr=pcbnew.PCB_TRACK(b);tr.SetStart(point(a));tr.SetEnd(point(z));tr.SetWidth(pcbnew.FromMM(t['width']));tr.SetLayer(b.GetLayerID(t['layer']));tr.SetNetCode(nets[t['net']]);b.Add(tr)
    for v in layout['vias']:
        vi=pcbnew.PCB_VIA(b);vi.SetPosition(point([v['x'],v['y']]));vi.SetWidth(pcbnew.FromMM(v['d']));vi.SetDrill(pcbnew.FromMM(v['drill']));vi.SetViaType(pcbnew.VIATYPE_THROUGH);vi.SetLayerPair(pcbnew.F_Cu,pcbnew.B_Cu);vi.SetNetCode(nets[v['net']]);b.Add(vi)
    if model.get('plane_export'):pcbnew.ZONE_FILLER(b).Fill(b.Zones())
    dest=out/'routed.kicad_pcb';pcbnew.SaveBoard(str(dest),b)
    import shutil
    for suffix in ['.kicad_pro','.kicad_dru']:
        src=(case/'unrouted.kicad_pcb').with_suffix(suffix)
        if src.exists():shutil.copyfile(src,dest.with_suffix(suffix))
    verify(name,out)
    print('exported',dest,flush=True)


def geometry(b):
    def pad(p):
        s=p.GetSize(pcbnew.F_Cu);d=p.GetDrillSize();o=p.GetOffset(pcbnew.F_Cu)
        return (p.GetNumber(),p.GetNetname(),p.GetPosition().x,p.GetPosition().y,s.x,s.y,d.x,d.y,
                p.GetShape(pcbnew.F_Cu),p.GetOrientation().AsDegrees(),p.GetAttribute(),str(p.GetLayerSet().FmtHex()),o.x,o.y,
                p.GetRoundRectCornerRadius(pcbnew.F_Cu))
    def part(f):
        return (f.GetReference(),f.GetValue(),str(f.GetFPID().GetLibItemName()),f.GetPosition().x,f.GetPosition().y,
                f.GetOrientation().AsDegrees(),f.GetLayer(),tuple(sorted(pad(p) for p in f.Pads())))
    return {'footprints':sorted(part(f) for f in b.GetFootprints()),'layers':[b.GetLayerName(i) for i in b.GetEnabledLayers().CuStack()]}


def verify(name,out):
    out=pathlib.Path(out);case=out/'input' if (out/'input/model.json').exists() else ROOT/'cases'/name
    before=pcbnew.LoadBoard(str(case/'unrouted.kicad_pcb'));after=pcbnew.LoadBoard(str(out/'routed.kicad_pcb'))
    a,z=geometry(before),geometry(after)
    static_equal=immutable_board_tokens(case/'unrouted.kicad_pcb')==immutable_board_tokens(out/'routed.kicad_pcb')
    rules_equal=all(not (case/('unrouted'+suffix)).exists() or (case/('unrouted'+suffix)).read_bytes()==(out/('routed'+suffix)).read_bytes() for suffix in ['.kicad_pro','.kicad_dru'])
    data={'passed':a==z and rules_equal and static_equal,'geometry_equal':a==z,'rule_files_identical':rules_equal,
          'static_board_content_equal':static_equal,
          'footprints':len(a['footprints']),'pads':sum(len(f[-1]) for f in a['footprints']),'layers':a['layers']}
    (out/'geometry-check.json').write_text(json.dumps(data,indent=2));print('geometry',data,flush=True)


def immutable_board_tokens(path):
    """Compare all persisted content except newly generated copper.

    KiCad can add display aliases to the layer table on reload; physical copper
    layers are compared independently. Every footprint attribute, drawing,
    outline, mask aperture and rule-area polygon otherwise remains exact.
    """
    tokens=re.findall(r'"(?:\\.|[^"\\])*"|[()]|[^\s()]+',path.read_text())
    out=[];depth=0;skip_depth=None
    for i,t in enumerate(tokens):
        if t=='(':
            depth+=1
            if skip_depth is None and ((depth==2 and tokens[i+1] in ('segment','via','arc','layers'))
                                       or tokens[i+1] in ('filled_polygon','fill_segments')):skip_depth=depth
        if skip_depth is None:out.append(t)
        if t==')':
            if depth==skip_depth:skip_depth=None
            depth-=1
    # KiCad may reorder whole footprints when UUIDs share a sort key. Object
    # order is not board geometry; preserve each object's full contents while
    # comparing the top-level objects as a canonical multiset.
    chunks=[];depth=0;start=0
    for i,t in enumerate(out):
        if t=='(':
            depth+=1
            if depth==2:start=i
        elif t==')':
            if depth==2:chunks.append(tuple(out[start:i+1]))
            depth-=1
    return tuple(sorted(chunks))


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','export','verify']);p.add_argument('case');p.add_argument('--out');p.add_argument('--source-dir')
    a=p.parse_args()
    if a.command=='prepare':prepare(a.case,a.source_dir)
    elif a.command=='export':export(a.case,a.out)
    else:verify(a.case,a.out)
    sys.stdout.flush();sys.stderr.flush();os._exit(0)
