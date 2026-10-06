#!/usr/bin/env python3
"""Turn native connectivity feedback into a bounded routing candidate."""
import argparse,dataclasses,json,pathlib,shutil,sys,hashlib
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]/'python'))
from model import build
from supplied_model import build as supplied
from gridroute.router import Options,RoutingController
import gridroute as gr


def repair(source,out,pitch=.025,deadline=30):
    if out.exists():raise FileExistsError('Use a fresh repair output directory')
    out.mkdir(parents=True);shutil.copytree(source/'input',out/'input')
    model=json.loads((out/'input/model.json').read_text());layout=json.loads((source/'layout.json').read_text())
    prior=json.loads((source/'routing.json').read_text());native=json.loads((source/'routed-drc.json').read_text())
    nets=[]
    known={p['net'] for p in model['pads'] if p['net']}
    for gap in native.get('unconnected_items',[]):
        for item in gap.get('items',[]):
            matches=[n for n in known if '['+n+']' in item.get('description','')]
            if len(matches)>1:raise ValueError('Ambiguous native item net name')
            if matches and matches[0] not in nets:nets.append(matches[0])
    if not nets:raise ValueError('Native feedback did not identify a known net; do not guess')
    builder=supplied if 'design_rules' in model and 'rule_areas' in model else build
    factory=lambda p:builder(model,pitch=p,margin=p*.75,copper_grow=p*.5)
    source_board=factory(pitch);source_board.replace_copper(layout['tracks'],layout['vias'])
    opts=Options(pitch=pitch,fallback_pitches=(),deadline_seconds=deadline,max_expansions=500000,
                 candidate_paths=3,plateau_budget=12,late_fanout=True,recovery_passes=3,cleanup=False)
    controller=RoutingController(factory,opts);board,events=controller.repair(source_board,nets)
    board.save(str(out/'layout.json'),model['size']);repair_seconds=events[-1]['elapsed_seconds']
    files=[ROOT/'repair.py',ROOT/'model.py',ROOT/'supplied_model.py',ROOT.parents[1]/'python/gridroute/board.py',ROOT.parents[1]/'python/gridroute/router.py']
    result=dict(case=model['name'],mode='native_feedback',source_candidate=str(source),
                native_feedback_nets=nets,options=dataclasses.asdict(opts),events=events,missing=board.unrouted(),
                routing_seconds=prior['routing_seconds']+repair_seconds,repair_seconds=repair_seconds,
                metal=gr.metal_available(),device=gr.device_name(),
                source_sha256={str(p.relative_to(ROOT.parents[1])):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                validation='candidate; fresh native DRC and geometry acceptance required')
    (out/'routing.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k not in ['events','options','source_sha256']}),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path,required=True)
    p.add_argument('--pitch',type=float,default=.025);p.add_argument('--deadline',type=float,default=30);a=p.parse_args()
    repair(a.source,a.out,a.pitch,a.deadline)
