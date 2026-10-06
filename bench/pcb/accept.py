"""Compare native CAD findings with the immutable stripped-input baseline."""
import collections,json,pathlib,argparse

def key(v):
 return (v['type'],v.get('severity'),tuple(sorted((i.get('description',''),round(i.get('pos',{}).get('x',0),4),round(i.get('pos',{}).get('y',0),4)) for i in v.get('items',[]))))

def accept(case,out):
 root=out/'input' if (out/'input/model.json').exists() else pathlib.Path(__file__).resolve().parent/'cases'/case
 baseline=json.loads((root/'baseline-drc.json').read_text());drc=json.loads((out/'routed-drc.json').read_text())
 old=collections.Counter(key(v) for v in baseline['violations']);added=[]
 for v in drc['violations']:
  k=key(v)
  if old[k]:old[k]-=1
  else:added.append(v)
 geometry=json.loads((out/'geometry-check.json').read_text())
 missing=len(drc.get('unconnected_items',[]));errors=[v for v in added if v.get('severity')=='error']
 result={'case':case,'passed':missing==0 and not errors and geometry['passed'],
         'native_unconnected_items':missing,'new_errors':errors,'new_warnings':[v for v in added if v.get('severity')!='error'],
         'baseline_violation_counts':dict(collections.Counter(v['type'] for v in baseline['violations'])),
         'final_violation_counts':dict(collections.Counter(v['type'] for v in drc['violations'])),
         'geometry':geometry,'differential_pairs':'not qualified; no coupled-routing or skew gate configured'}
 (out/'acceptance.json').write_text(json.dumps(result,indent=2))
 print(json.dumps(dict(case=case,passed=result['passed'],native_unconnected_items=missing,
                       new_errors=dict(collections.Counter(v['type'] for v in errors)),
                       new_warnings=len(result['new_warnings']),geometry_passed=geometry['passed'])),flush=True)
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args()
 raise SystemExit(0 if accept(a.case,a.out)['passed'] else 1)
