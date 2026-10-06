"""Make a separately named fixture with more copper layers (run with KiCad's Python).

usage: make_layer_variant.py SRC_CASE DST_CASE LAYERS
The board, project and rules are copied; the copper layer count is raised (new inner layers are
signal layers, existing zones keep their layers) and rule areas that covered every copper layer are
extended to the new ones, so keep-outs stay keep-outs. model.json gets the new layer list. Run
`kicad-cli pcb drc` on the new unrouted board afterwards to record its baseline. This is a stack-up
change of the design, not an unchanged-input result.
"""
import json,pathlib,shutil,sys,os
import pcbnew
ROOT=pathlib.Path(__file__).resolve().parent
src,dst,count=sys.argv[1],sys.argv[2],int(sys.argv[3])
s,d=ROOT/'cases'/src,ROOT/'cases'/dst
d.mkdir(parents=True,exist_ok=True)   # an existing variant is regenerated from the source case
for f in s.iterdir():
    if f.name.startswith('unrouted.') or f.name=='model.json':shutil.copyfile(f,d/f.name)
b=pcbnew.LoadBoard(str(d/'unrouted.kicad_pcb'))
old=[b.GetLayerName(i) for i in b.GetEnabledLayers().CuStack()]
b.SetCopperLayerCount(count)
new=[b.GetLayerName(i) for i in b.GetEnabledLayers().CuStack()]
extended=0
for z in b.Zones():
    if z.GetIsRuleArea() and set(z.GetLayerSet().CuStack())>={b.GetLayerID(n) for n in old}:
        ls=pcbnew.LSET(z.GetLayerSet());ls.AddLayerSet(pcbnew.LSET.AllCuMask(count));z.SetLayerSet(ls);extended+=1
pcbnew.SaveBoard(str(d/'unrouted.kicad_pcb'),b)
m=json.loads((d/'model.json').read_text())
added=[n for n in new if n not in old]
m['layers']=new
for a in m['rule_areas']:
    if set(a['layers'])>=set(old):a['layers']=list(new)
m['variant']=(m.get('variant','')+f' Stack-up raised from {len(old)} to {count} copper layers; added signal layers {added}.').strip()
m['source_case']=src
(d/'model.json').write_text(json.dumps(m,indent=2))
print(json.dumps(dict(old=old,new=new,rule_areas_extended=extended)),flush=True)
sys.stdout.flush();os._exit(0)
