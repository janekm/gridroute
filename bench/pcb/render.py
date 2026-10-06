#!/usr/bin/env python3
"""Render a routing result per copper layer: pads, tracks, vias and unrouted pads (red rings).

usage: render.py RESULT_DIR [--out PNG] [--scale PX_PER_MM] [--box X0 Y0 X1 Y1] [--layers L ...]
Plane-net vias get a green dot.
RESULT_DIR holds input/model.json, layout.json and routing.json (bench/pcb/run.py output).
"""
import argparse,json,math,pathlib
from PIL import Image,ImageDraw

COLORS=['#f37f72','#8d96ed','#88ca85','#e1b773','#a7a7d8','#64bddd']

def render(out_dir,png=None,scale=24,box=None,only=None):
    out_dir=pathlib.Path(out_dir)
    m=json.loads((out_dir/'input/model.json').read_text());lay=json.loads((out_dir/'layout.json').read_text())
    routing=json.loads((out_dir/'routing.json').read_text()) if (out_dir/'routing.json').exists() else {}
    missing={(r,n) for _,r,n in routing.get('missing',[])}
    planes=set(routing.get('planes',{}))|{z['net'] for z in m.get('plane_intent',[]) if m.get('plane_export')}
    layers=m['layers'];W,H=m['size'];bx,by=0,0
    if box:bx,by=box[0],box[1];W,H=box[2]-box[0],box[3]-box[1]
    shown=[L for L in layers if not only or L in only]
    cols=3 if len(shown)>4 else len(shown);rows=math.ceil(len(shown)/cols)
    pw,ph=int(W*scale)+20,int(H*scale)+40
    im=Image.new('RGB',(pw*cols,ph*rows),'#0d1621');d=ImageDraw.Draw(im)
    canvas=im
    for k,L in enumerate(shown):
        im=Image.new('RGB',(pw,ph),'#0d1621');d=ImageDraw.Draw(im)
        ox,oy=10-bx*scale,30-by*scale
        pt=lambda p:(ox+p[0]*scale,oy+p[1]*scale)
        d.text((10,6),L,fill='#e8f0f6')
        d.polygon([pt(p) for p in m['outline']],fill='#182c34')
        li=layers.index(L);col=COLORS[0] if li==0 else COLORS[-1] if li==len(layers)-1 else COLORS[li%(len(COLORS)-1)]
        for p in m['pads']:
            on=layers if p['layers']=='all' else [layers[0] if p['layers']=='F' else layers[-1]]
            if L not in on:continue
            w,h=p['w'],p['h']
            fill='#e5d8a4' if p['net'] else '#6b6650'
            d.rectangle([pt((p['x']-w/2,p['y']-h/2)),pt((p['x']+w/2,p['y']+h/2))],fill=fill)
        for t in lay['tracks']:
            if t['layer']!=L or len(t['pts'])<2:continue
            d.line([pt(q) for q in t['pts']],fill=col,width=max(1,round(t['width']*scale)),joint='curve')
        for v in lay['vias']:
            x,y=pt((v['x'],v['y']));r=v['d']/2*scale
            d.ellipse((x-r,y-r,x+r,y+r),outline='#e5d8a4',width=max(1,int(scale*.06)))
            if v['net'] in planes:d.ellipse((x-r/3,y-r/3,x+r/3,y+r/3),fill='#5a9')
        for p in m['pads']:
            if (p['ref'],p['num']) in missing:
                x,y=pt((p['x'],p['y']));r=.45*scale
                d.ellipse((x-r,y-r,x+r,y+r),outline='#ff3030',width=2)
        canvas.paste(im,((k%cols)*pw,(k//cols)*ph))
    im=canvas
    png=pathlib.Path(png or out_dir/'copper.png');im.save(png);return png

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('result');a.add_argument('--out');a.add_argument('--scale',type=float,default=24)
    a.add_argument('--box',type=float,nargs=4,help='x0 y0 x1 y1 in mm');a.add_argument('--layers',nargs='*')
    x=a.parse_args();print(render(x.result,x.out,x.scale,x.box,x.layers))
