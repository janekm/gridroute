"""Exact polygon/pad adapter for the immutable PCB benchmark models."""
import copy
from functools import lru_cache
import types
import numpy as np
import gridroute as gr
from gridroute import board

def build(model,pitch=0.05,margin=None,pad_grow=None,shifts=None,copper_grow=0.):
    margin = pitch/2 if margin is None else margin
    pad_grow = pitch/2 if pad_grow is None else pad_grow
    model=copy.deepcopy(model)
    if shifts:
        for p in model['pads']:
            dx,dy=shifts.get(p['original_ref'],(0,0));p['x']+=dx;p['y']+=dy
    rules=model['rules']
    net_classes=model.get('net_class',{})
    board.configure(layers=model['layers'],pitch=pitch,margin=margin,pad_grow=pad_grow,classes=model.get('classes',{'Default':tuple(rules)}),
                    net_class=lambda n:net_classes.get(n,'Default'),
                    edge_clear=0.5,hole_clear=0.25,copper_grow=copper_grow,cache='off',spec=0,search='hybrid',budget=50000,field=2)
    assert gr.available(), 'native library failed to load'
    comps={p['ref']:dict(footprint=p['ref']) for p in model['pads']}
    pin_net={(p['ref'],'1'):p['net'] for p in model['pads'] if p['net']}
    fps={p['ref']:dict(pads=[dict(p,x=0,y=0)],courtyard=[0,0,0,0]) for p in model['pads']}
    bd=board.Board(*model['size'],comps,pin_net,fps)
    # Respect the real polygon, including rounded and bevelled corners. The pad/trace
    # clearance dilation supplies the remaining distance to the outline.
    edge_mask=outline_mask(tuple(map(tuple,model['outline'])),bd.ny,bd.nx,pitch,.5,board.CLEAR,
                           tuple(tuple(map(tuple,h)) for h in model.get('outline_holes',[])))
    def edges(self):
        board.Board._edges(self);self.occ[:,edge_mask]=board.KEEP
    bd._edges=types.MethodType(edges,bd);bd._edges()
    for p in model['pads']: bd.place(p['ref'],p['x'],p['y'],gap=0)
    # Unassigned solder-mask apertures must not expose unrelated routed copper.
    # Treat them conservatively as obstacles, never as electrical connectivity.
    for g in model.get('copper_graphics',[])+model.get('mask_graphics',[]):
        if 'circle' in g:bd.keepout_ring(**g['circle'],layers=[g['layer']])
        else:bd.keepout_polygon(g['points'],layers=[g['layer']],holes=g['holes'])
    for g in model.get('edge_graphics',[]):
        if 'circle' in g:bd.keepout_ring(**g['circle'],clearance=.5)
        else:bd.keepout_polygon(g['points'],holes=g['holes'],clearance=.5)
    return bd


@lru_cache(maxsize=12)
def outline_mask(poly,ny,nx,pitch,edge,clearance,holes=()):
    yy,xx=np.mgrid[:ny,:nx].astype(float);xx*=pitch;yy*=pitch
    inside=np.zeros(xx.shape,dtype=bool);dist=np.full(xx.shape,np.inf)
    for ring_index,ring in enumerate((poly,)+holes):
        ring_inside=np.zeros(xx.shape,dtype=bool)
        for a,b in zip(ring,ring[1:]+ring[:1]):
            ax,ay=a;bx,by=b;dx=bx-ax;dy=by-ay
            if dy:
                ring_inside^=((ay>yy)!=(by>yy)) & (xx < (bx-ax)*(yy-ay)/dy+ax)
            den=dx*dx+dy*dy
            t=np.clip(((xx-ax)*dx+(yy-ay)*dy)/den,0,1) if den else 0
            dist=np.minimum(dist,(xx-ax-t*dx)**2+(yy-ay-t*dy)**2)
        if ring_index==0:inside=ring_inside
        else:inside &= ~ring_inside
    edge_mask=(~inside)|(dist<(max(0,edge-clearance)+pitch)**2)
    edge_mask.flags.writeable=False
    return edge_mask
