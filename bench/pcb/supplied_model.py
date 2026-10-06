"""Conservative routing model for a supplied KiCad project with explicit rules.

The current adapter supports this fixture's layer minimums and permitted outer
pad escape widths. RF's local clearance relaxation is not used: retaining the
larger class clearance is conservative and may leave a route incomplete.
Copper planes are stripped for the cold-routing variant; all six layers remain
available. A plane-export variant keeps its inner-layer zones: their nets are
carried by vias to the plane and the plane layers carry no signal tracks. Native DRC against the unchanged supplied project is the acceptance
oracle. This is not differential-pair or impedance qualification.
"""
import copy,types
import gridroute as gr
from gridroute import board
from model import outline_mask

def plane_layers(model):
 """{net: inner layer} for inner-layer zones the exported board keeps and refills."""
 if not model.get('plane_export'):return {}
 return {z['net']:z['layers'][0] for z in model.get('plane_intent',[])
         if len(z['layers'])==1 and z['layers'][0] in model['layers'][1:-1]}

def build(model,pitch=.05,margin=None,copper_grow=0.,planes=None,necks=True):
 """planes: {net: layer}; default plane_layers(model). Plane layers carry no signal tracks."""
 m=copy.deepcopy(model);rules=m['design_rules'];classes=m['classes'];net_classes=m['net_class']
 planes=plane_layers(m) if planes is None else dict(planes)
 signal=[l for l in m['layers'] if l not in planes.values()]
 inner=signal[1:-1]
 layer_widths={n:{l:max(r[0],.152) for l in inner} for n,r in classes.items()}
 board.configure(layers=signal,planes=planes,pitch=pitch,margin=pitch/2 if margin is None else margin,pad_grow=pitch/2,copper_grow=copper_grow,classes=classes,
  net_class=lambda n:net_classes.get(n,'Default'),layer_widths=layer_widths,layer_clearances={l:.152 for l in inner},
  edge_clear=rules['min_copper_edge_clearance'],hole_clear=rules['min_hole_clearance'],hole_to_hole=rules['min_hole_to_hole'],cache='off',spec=0,search='hybrid',budget=50000,field=2)
 assert gr.available()
 neck_refs={'U1','J3','U3','U4','U6','U7','U8','U9','U10'}
 for p in m['pads']:
  cls=net_classes.get(p['net'],'Default')
  if p['layers']!='all' and cls in ('PWR','GND') and p['original_ref'] in neck_refs:
   p['escape_width']=.12 if cls=='PWR' and p['original_ref']=='U3' else .15
 comps={p['ref']:dict(footprint=p['ref']) for p in m['pads']}
 pins={(p['ref'],'1'):p['net'] for p in m['pads'] if p['net']}
 fps={p['ref']:dict(pads=[dict(p,x=0,y=0)],courtyard=[0,0,0,0]) for p in m['pads']}
 bd=board.Board(*m['size'],comps,pins,fps)
 edge=outline_mask(tuple(map(tuple,m['outline'])),bd.ny,bd.nx,pitch,rules['min_copper_edge_clearance'],board.CLEAR)
 def edges(self):board.Board._edges(self);self.occ[:,edge]=board.KEEP
 bd._edges=types.MethodType(edges,bd);bd._edges()
 for p in m['pads']:bd.place(p['ref'],p['x'],p['y'],gap=0)
 for area in m['rule_areas']:
  layers=[l for l in area['layers'] if l in signal]
  if not layers:continue
  for poly in area['polygons']:bd.keepout_polygon(poly,layers=layers,tracks=area['tracks'],vias=area['vias'])
 for g in m.get('copper_graphics',[])+m.get('mask_graphics',[]):bd.keepout_polygon(g['points'],layers=[g['layer']],holes=g['holes'])
 # The project's neck-down rules: PWR/GND tracks touching these courtyards may be
 # 0.15 mm on outer layers, PWR at U3 0.12 mm. Only the part's own side is used.
 if necks:
  power=[n for n,c in net_classes.items() if c in ('PWR','GND')]
  pwr=[n for n in power if net_classes[n]=='PWR']
  for part in m.get('parts',[]):
   if part['ref'] not in neck_refs:continue
   for layer,poly in part.get('courtyards',{}).items():
    if layer not in (signal[0],signal[-1]):continue
    if part['ref']=='U3':
     bd.add_neck(poly,.12,pwr,layers=[layer]);bd.add_neck(poly,.15,[n for n in power if n not in pwr],layers=[layer])
    else:bd.add_neck(poly,.15,power,layers=[layer])
 return bd
