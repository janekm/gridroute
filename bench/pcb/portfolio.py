#!/usr/bin/env python3
"""Route one fixture with several option variants in parallel processes and keep the best candidate.

Each variant runs bench/pcb/run.py into OUT/variants/<name>. The best is chosen by the in-process continuous
check recorded in routing.json (fewest clearance findings, then the lowest copper-component deficit, then fewest
vias) and its result files are copied to OUT, so OUT looks like a single run.py result. Native CAD validation is
still required. Board configuration is process-global, so variants are separate processes.
"""
import argparse,json,os,pathlib,shutil,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parent

DEFAULT=[('base',{}),('s1',{'seed':1}),('s2',{'seed':2}),('s3',{'seed':3}),
         ('late',{'fanout':False,'late_fanout':True}),('late_s1',{'fanout':False,'late_fanout':True,'seed':1})]

def rank(out):
    data=json.loads((out/'routing.json').read_text());f=data.get('final_check') or {}
    return (f.get('violations',10**6),f.get('deficit',10**6),data['vias'])

def main():
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('--out',type=pathlib.Path,required=True)
    p.add_argument('--options-json',default='{}',help='options shared by every variant')
    p.add_argument('--variants-json',help='[[name, {options}], ...]; default: a seed/fanout portfolio')
    p.add_argument('--jobs',type=int,default=min(6,os.cpu_count() or 1))
    a=p.parse_args()
    if any((a.out/n).exists() for n in ['routing.json','layout.json']):p.error('Use a fresh output directory')
    shared=json.loads(a.options_json);variants=json.loads(a.variants_json) if a.variants_json else DEFAULT
    (a.out/'variants').mkdir(parents=True,exist_ok=True);env=dict(os.environ,PYTHONHASHSEED='0')
    pending=list(variants);running=[];done={};started=time.perf_counter()
    while pending or running:
        while pending and len(running)<a.jobs:
            name,opts=pending.pop(0);o=dict(shared);o.update(opts);dest=a.out/'variants'/name
            log=open(a.out/'variants'/(name+'.log'),'w')
            running.append((name,dest,log,subprocess.Popen([sys.executable,str(ROOT/'run.py'),a.case,'--out',str(dest),
                                                            '--options-json',json.dumps(o)],stdout=log,stderr=subprocess.STDOUT,env=env)))
        time.sleep(.5)
        for r in list(running):
            if r[3].poll() is not None:
                running.remove(r);r[2].close();done[r[0]]=r[3].returncode
    ok=[n for n,c in done.items() if c==0 and (a.out/'variants'/n/'routing.json').exists()]
    if not ok:print(json.dumps(dict(case=a.case,error='no variant finished',codes=done)));return 1
    ranking=sorted(ok,key=lambda n:rank(a.out/'variants'/n));best=a.out/'variants'/ranking[0]
    for item in ['input','layout.json']:
        src=best/item
        (shutil.copytree if src.is_dir() else shutil.copyfile)(src,a.out/item)
    data=json.loads((best/'routing.json').read_text())
    data['portfolio']=dict(selected=ranking[0],ranking=[[n,*rank(a.out/'variants'/n)] for n in ranking],
                           wall_seconds=time.perf_counter()-started,jobs=a.jobs)
    (a.out/'routing.json').write_text(json.dumps(data,indent=2))
    print(json.dumps(dict(case=a.case,selected=ranking[0],ranking=data['portfolio']['ranking'],
                          wall_seconds=round(data['portfolio']['wall_seconds'],1))),flush=True)
    return 0

if __name__=='__main__':raise SystemExit(main())
