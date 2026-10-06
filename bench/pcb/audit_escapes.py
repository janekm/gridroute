#!/usr/bin/env python3
"""Escape audit for a fixture: can every connected pin get out of its package (and to a via where its net
must change layer)? Run on the empty board, then after the controller's pre-routing stages.

usage: audit_escapes.py CASE [--options-json '{...}'] [--refs U10 J3]
"""
import argparse,json,pathlib,sys,os
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]/'python'));sys.path.insert(0,str(ROOT))
os.environ.setdefault('GRIDROUTE_CACHE','off')
from model import build
from supplied_model import build as build_supplied
from gridroute.router import Options,RoutingController
from gridroute.escape import audit_escapes

def summary(name,r):
    bad={k:v for k,v in r['packages'].items() if v['exit']<v['pins'] or v['via']<v['via_needed']}
    print(f"{name}: {r['pins']} pins; no exit {len(r['no_exit'])}, no via site {len(r['no_via'])}, "
          f"no joint via site {len(r['joint_unplaced'])}")
    for k,v in sorted(bad.items()):print('  ',k,v)
    for kind in ('no_exit','no_via','joint_unplaced'):
        if r[kind]:print('  ',kind,[tuple(x) for x in r[kind]])

p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('--options-json',default='{}');p.add_argument('--refs',nargs='*')
a=p.parse_args()
model=json.loads((ROOT/'cases'/a.case/'model.json').read_text())
builder=build_supplied if 'design_rules' in model and 'rule_areas' in model else build
bd=builder(model,pitch=.05,margin=.025)
r=audit_escapes(bd,a.refs);summary('empty board',r)
o=json.loads(a.options_json);o.setdefault('fallback_pitches',[])
ctl=RoutingController(lambda pitch:builder(model,pitch=pitch,margin=pitch*.5),Options(**o));ctl._begin()
bd=ctl._fresh(ctl.options.pitch)
# the controller's pre-routing stages, in run() order
if ctl.options.escape_hold:
    from gridroute.escape import dense_packages,plan_dogbones
    n,un=plan_dogbones(bd,dense_packages(bd,ctl.options.fine_pitch_mm),commit=False)
    print(f'held escapes: {n} planned, {len(un)} without one: {[tuple(x) for x in un]}')
    r=audit_escapes(bd,a.refs);summary('after planning held escapes',r)
if ctl.options.escape_reserve_mm>0:ctl._reserve_escapes(bd)
if ctl.options.neck_escapes and bd.__dict__.get('necks'):
    for n in [n for n in ctl._order(bd) if any(q.get('escape_width') for q in bd.pads_of(n))]:ctl._connect(bd,n)
r=audit_escapes(bd,a.refs);summary('after neck-down nets',r)
for n in sorted(bd.config['planes']):ctl._drop_net(bd,n)
r=audit_escapes(bd,a.refs);summary('after plane drops (general routing starts here)',r)
