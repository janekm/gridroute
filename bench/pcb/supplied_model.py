"""Conservative routing model for a supplied KiCad project with explicit rules.

The current adapter supports this fixture's layer minimums and permitted outer
pad escape widths. RF's local clearance relaxation is not used: retaining the
larger class clearance is conservative and may leave a route incomplete.
Copper planes are stripped for this cold-routing variant; all six layers remain
available. Native DRC against the unchanged supplied project is the acceptance
oracle. This is not differential-pair or impedance qualification.
"""
import copy,types
import gridroute as gr
from gridroute import board
from model import outline_mask

def build(model,pitch=.05,margin=None,copper_grow=0.):
 m=copy.deepcopy(model);rules=m['design_rules'];classes=m['classes'];net_classes=m['net_class']
 inner=m['layers'][1:-1]
 layer_widths={n:{l:max(r[0],.152) for l in inner} for n,r in classes.items()}
 board.configure(layers=m['layers'],pitch=pitch,margin=pitch/2 if margin is None else margin,pad_grow=pitch/2,copper_grow=copper_grow,classes=classes,
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
  for poly in area['polygons']:bd.keepout_polygon(poly,layers=area['layers'],tracks=area['tracks'],vias=area['vias'])
 for g in m.get('copper_graphics',[])+m.get('mask_graphics',[]):bd.keepout_polygon(g['points'],layers=[g['layer']],holes=g['holes'])
 return bd
