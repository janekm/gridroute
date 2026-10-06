#!/usr/bin/env python3
"""Route an immutable fixture, saving the exact inputs beside the candidate."""
import argparse,dataclasses,hashlib,json,os,pathlib,shutil,sys,resource,time
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]/'python'))
os.environ.setdefault('GRIDROUTE_CACHE','off')
from model import build
from supplied_model import build as build_supplied
import gridroute as gr
from gridroute.router import Options,RoutingController


def capture_fixture(case,out):
    source=ROOT/'cases'/case;dest=out/'input';dest.mkdir(parents=True,exist_ok=True)
    for name in ['model.json','unrouted.kicad_pcb','unrouted.kicad_pro','unrouted.kicad_dru','baseline-drc.json','variant-policy.json']:
        if (source/name).exists():shutil.copyfile(source/name,dest/name)
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.iterdir() if p.is_file() and p.name!='sha256.json'}
    (dest/'sha256.json').write_text(json.dumps(hashes,indent=2))
    return json.loads((dest/'model.json').read_text())


def run(case,out,options,margin_factor=.5,copper_grow_factor=0.,continuous_check=False):
    if any((out/name).exists() for name in ['routing.json','layout.json','acceptance.json','routed.kicad_pcb']):
        raise FileExistsError('Use a fresh output directory; benchmark evidence must not be overwritten')
    out.mkdir(parents=True,exist_ok=True);model=capture_fixture(case,out)
    builder=build_supplied if 'design_rules' in model and 'rule_areas' in model else build
    ctl=RoutingController(lambda pitch:builder(model,pitch=pitch,margin=pitch*margin_factor,copper_grow=pitch*copper_grow_factor),options)
    started=time.perf_counter();bd,events=ctl.run();continuous=None;attempts=[]
    if continuous_check:
        from gridroute.continuous import repair_connectivity
        factory=lambda pitch:builder(model,pitch=pitch,margin=pitch*.75,copper_grow=pitch*.5)
        bd,continuous,attempts=repair_connectivity(bd,factory,deadline_seconds=min(30.,options.deadline_seconds))
    elapsed=time.perf_counter()-started
    from gridroute.continuous import check
    final=check(bd);final_summary=dict(deficit=sum(max(0,c-1) for c in final['copper_components'].values()),
                                       violations=len(final['violations']),missing_nets=final['missing_nets'])
    bd.save(str(out/'layout.json'),model['size'])
    sources=[ROOT/'model.py',ROOT/'supplied_model.py',ROOT.parents[1]/'python/gridroute/board.py',ROOT.parents[1]/'python/gridroute/router.py',ROOT.parents[1]/'python/gridroute/continuous.py',ROOT/'run.py']
    data={'case':case,'model_sha256':hashlib.sha256((out/'input/model.json').read_bytes()).hexdigest(),
          'source_sha256':{str(p.relative_to(ROOT.parents[1])):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
          'native':gr.available(),'metal':gr.metal_available(),'device':gr.device_name(),'options':dataclasses.asdict(options),
          'margin_factor':margin_factor,'copper_grow_factor':copper_grow_factor,
          'events':events,'final_check':final_summary,'continuous_check':continuous,'continuous_attempts':attempts,'routing_seconds':elapsed,'missing':bd.unrouted(),
          'tracks':len(bd.tracks),'vias':len(bd.vias),'routable_nets':model['routable_nets'],
          'complete_nets':model['routable_nets']-len({p[0] for p in bd.unrouted()}),
          'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
          'validation':'candidate; native CAD validation required'}
    (out/'routing.json').write_text(json.dumps(data,indent=2))
    print(json.dumps({k:v for k,v in data.items() if k not in ['events','missing','source_sha256','options']}),flush=True)
    return data


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('--out',type=pathlib.Path,required=True)
    p.add_argument('--deadline',type=float,default=60);p.add_argument('--no-fanout',action='store_true')
    p.add_argument('--no-refinement',action='store_true');p.add_argument('--pitch',type=float,default=.05)
    p.add_argument('--options-json',default='{}',help='Options fields, plus margin_factor and copper_grow_factor')
    a=p.parse_args();kw=dict(pitch=a.pitch,fallback_pitches=() if a.no_refinement else (a.pitch/2,),fanout=not a.no_fanout,deadline_seconds=a.deadline)
    kw.update(json.loads(a.options_json));margin=kw.pop('margin_factor',.5);grow=kw.pop('copper_grow_factor',0.);continuous=kw.pop('continuous_check',False)
    run(a.case,a.out,Options(**kw),margin,grow,continuous)
