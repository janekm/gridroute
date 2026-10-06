#!/usr/bin/env python3
"""Serial native acceptance harness. A process exit is never routing acceptance."""
import argparse,collections,json,pathlib,subprocess,sys,time,os
ROOT=pathlib.Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=pathlib.Path,required=True);p.add_argument('--cases',nargs='*')
    p.add_argument('--kicad-python',default='/Applications/KiCad.app/Contents/Frameworks/Python.framework/Versions/Current/bin/python3')
    p.add_argument('--kicad-cli',default='/Applications/KiCad.app/Contents/MacOS/kicad-cli');p.add_argument('--deadline',type=float,default=60)
    p.add_argument('--options-json',default='{}')
    p.add_argument('--native-feedback',action='store_true',help='Retry at most four native gaps, at two bounded pitches')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    cases=json.loads((ROOT/'manifest.json').read_text())['cases'];rows=[]
    chosen=[c for c in cases if c['id'] in a.cases] if a.cases else [c for c in cases if c.get('enabled',True)]
    if a.cases and set(a.cases)-{c['id'] for c in chosen}:p.error('Unknown case name')
    if any((a.out/c['id']).exists() for c in chosen):p.error('Use a fresh output directory for every suite run')
    env=dict(os.environ,PYTHONHASHSEED='0',GRIDROUTE_CACHE='off',GRIDROUTE_NO_BUILD='1')
    for case in chosen:
        name=case['id'];out=a.out/name;out.mkdir();options=dict(case.get('routing_options',{}));options.update(json.loads(a.options_json));options.setdefault('deadline_seconds',a.deadline)
        commands=[([sys.executable,str(ROOT/'run.py'),name,'--out',str(out),'--options-json',json.dumps(options)],options['deadline_seconds']+120),
                  ([a.kicad_python,str(ROOT/'kicad_io.py'),'export',name,'--out',str(out)],60),
                  ([a.kicad_cli,'pcb','drc','--format','json','--severity-all','--all-track-errors','-o',str(out/'routed-drc.json'),str(out/'routed.kicad_pcb')],120),
                  ([sys.executable,str(ROOT/'accept.py'),name,'--out',str(out)],30)]
        row={'case':name,'steps':[]};start=time.monotonic()
        with (out/'run.log').open('w') as log:
            for cmd,limit in commands:
                try:code=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=limit).returncode
                except subprocess.TimeoutExpired:code=124
                row['steps'].append(dict(command=cmd,returncode=code))
                if code:break
        selected=out
        if a.native_feedback and (out/'acceptance.json').exists():
            current=json.loads((out/'acceptance.json').read_text());row['repair_attempts']=[]
            for pitch in [.025,.0125]:
                if not (0<current['native_unconnected_items']<=4 and not current['new_errors'] and current['geometry']['passed']):break
                trial=out/('repair-'+str(pitch))
                repair_commands=[([sys.executable,str(ROOT/'repair.py'),'--source',str(selected),'--out',str(trial),'--pitch',str(pitch),'--deadline','30'],90),
                                 ([a.kicad_python,str(ROOT/'kicad_io.py'),'export',name,'--out',str(trial)],60),
                                 ([a.kicad_cli,'pcb','drc','--format','json','--severity-all','--all-track-errors','-o',str(trial/'routed-drc.json'),str(trial/'routed.kicad_pcb')],120),
                                 ([sys.executable,str(ROOT/'accept.py'),name,'--out',str(trial)],30)]
                with (out/('repair-'+str(pitch)+'.log')).open('w') as log:
                    for cmd,limit in repair_commands:
                        try:code=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=limit).returncode
                        except subprocess.TimeoutExpired:code=124
                        if code:break
                improved=False
                if (trial/'acceptance.json').exists():
                    candidate=json.loads((trial/'acceptance.json').read_text())
                    improved=(candidate['geometry']['passed'] and not candidate['new_errors'] and
                              candidate['native_unconnected_items']<current['native_unconnected_items'])
                    if improved:selected=trial;current=candidate
                row['repair_attempts'].append(dict(pitch=pitch,improved=improved,returncode=code,out=str(trial)))
        row['wall_seconds']=time.monotonic()-start;row['selected_output']=str(selected)
        if (selected/'acceptance.json').exists():
            result=json.loads((selected/'acceptance.json').read_text());row.update(passed=result['passed'],unconnected=result['native_unconnected_items'],new_errors=len(result['new_errors']),new_warnings=len(result['new_warnings']))
        else:row['passed']=False;row['failure']='No native acceptance result'
        if (selected/'routing.json').exists():row['routing_seconds']=json.loads((selected/'routing.json').read_text())['routing_seconds']
        rows.append(row);(a.out/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps({k:v for k,v in row.items() if k!='steps'}),flush=True)
    return 0 if rows and all(r['passed'] for r in rows) else 1

if __name__=='__main__':raise SystemExit(main())
