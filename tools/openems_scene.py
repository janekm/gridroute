#!/usr/bin/env python3
"""Replay a fieldsolve scene in openEMS and compare port impedances.

    ~/opt/openEMS/venv/bin/python tools/openems_scene.py RUN_DIR [--keep]

RUN_DIR is a `fieldsolve run` output directory (scene.json, result.json). The scene is rebuilt with exactly
fieldsolve's mesh lines (no smoothing), copper sheets as zero-thickness metal polygons, dielectric slabs
extruded from the board outline, plated barrels as metal lines on the same mesh nodes fieldsolve uses,
lumped ports/elements on the same boxes, the same Gaussian excitation, PML thickness and end criterion.
Writes RUN_DIR/openems/result.json and prints the impedance comparison.

openEMS is GPLv3 and used here only as an external reference; nothing from it is linked into fieldsolve.
"""
import json
import math
import os
import pathlib
import shutil
import sys

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import EPS0

run = pathlib.Path(sys.argv[1]).resolve()
keep = '--keep' in sys.argv
sc = json.loads((run / 'scene.json').read_text())
fs = json.loads((run / 'result.json').read_text())
out = run / 'openems'
sim_path = out / 'sim'
if sim_path.exists():
    shutil.rmtree(sim_path)
sim_path.mkdir(parents=True)

gx, gy, gz = (np.array(sc['grid'][a]) for a in 'xyz')
fmax = sc['fmax']
FDTD = openEMS(NrTS=int(sc['max_steps']), EndCriteria=10 ** (sc['end_db'] / 10.0))
FDTD.SetGaussExcite(0.0, fmax)
FDTD.SetBoundaryCond([('PML_%d' % n) if n > 0 else 'PEC' for n in sc['pml']])
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(sc['unit'])
mesh.SetLines('x', gx)
mesh.SetLines('y', gy)
mesh.SetLines('z', gz)

# dielectrics
w_loss = 2 * math.pi * sc['f_loss']
outline = sc['outline']
x0, x1, y0, y1 = gx[0], gx[-1], gy[0], gy[-1]
for k, s in enumerate(sc['slabs']):
    m = CSX.AddMaterial('slab%d' % k, epsilon=s['eps'], kappa=w_loss * EPS0 * s['eps'] * s['tand'])
    if outline:
        # loops are outer first; the outer one is the board, the others cut-outs (even-odd)
        for li, loop in enumerate(outline):
            pts = np.array(loop).T
            if li == 0:
                m.AddLinPoly(points=pts, norm_dir='z', elevation=s['z0'], length=s['z1'] - s['z0'], priority=1)
            else:
                air = CSX.AddMaterial('cut%d_%d' % (k, li), epsilon=1.0)
                air.AddLinPoly(points=pts, norm_dir='z', elevation=s['z0'], length=s['z1'] - s['z0'], priority=2)
    else:
        m.AddBox([x0, y0, s['z0']], [x1, y1, s['z1']], priority=1)

metal = CSX.AddMetal('copper')
prio = 10
n_clear = 0
for sh in sc['sheets']:
    for p in sh['prims']:
        prio += 1
        if not p['dark']:
            n_clear += 1
            continue
        pts = np.array(p['pts']).T
        metal.AddPolygon(points=pts, norm_dir='z', elevation=sh['z'], priority=prio)
if n_clear:
    print('openems_scene: %d clear (LPC) primitives ignored' % n_clear, file=sys.stderr)


def snap(lines, v):
    return int(np.argmin(np.abs(lines - v)))


# barrels: z lines at the mesh nodes inside the radius (nearest node if none), as fieldsolve does
for b in sc['barrels']:
    (cx, cy), r = b['c'], b['r']
    ii = [i for i in range(len(gx)) if abs(gx[i] - cx) <= r]
    jj = [j for j in range(len(gy)) if abs(gy[j] - cy) <= r]
    nodes = [(i, j) for i in ii for j in jj if (gx[i] - cx) ** 2 + (gy[j] - cy) ** 2 <= r * r]
    if not nodes:
        nodes = [(snap(gx, cx), snap(gy, cy))]
    za, zb = gz[snap(gz, b['z0'])], gz[snap(gz, b['z1'])]
    for i, j in nodes:
        metal.AddCurve(np.array([[gx[i], gx[i]], [gy[j], gy[j]], [za, zb]]), priority=500)

ports = []
exc = sc['excite']
pnum = 0
for k, l in enumerate(sc['lumped']):
    lo, hi = list(l['lo']), list(l['hi'])
    if l['plus_at_hi']:
        start, stop = lo, hi
    else:
        start, stop = hi, lo
    d = 'xyz'[l['axis']]
    if l['port']:
        ports.append(FDTD.AddLumpedPort(pnum + 1, l['r'], start, stop, d, excite=1.0 if pnum == exc else 0, priority=1000))
        pnum += 1
    elif l['r'] == 0 and l['c'] == 0:
        metal.AddBox(lo, hi, priority=1000)
    else:
        kw = {}
        if l['r'] > 0:
            kw['R'] = l['r']
        if l['c'] > 0:
            kw['C'] = l['c']
        el = CSX.AddLumpedElement('el%d' % k, ny=d, caps=True, **kw)
        el.AddBox(lo, hi, priority=1000)

CSX.Write2XML(str(out / 'scene.xml'))
FDTD.Run(str(sim_path), verbose=0, cleanup=True)

freqs = np.array(fs['ports'][0]['freqs'])
res = {'freqs': freqs.tolist(), 'ports': []}
for k, p in enumerate(ports):
    p.CalcPort(str(sim_path), freqs, ref_impedance=50)
    z = p.uf_tot / p.if_tot
    res['ports'].append({'name': sc['lumped'][k]['name'], 'z_re': z.real.tolist(), 'z_im': z.imag.tolist()})
(out / 'result.json').write_text(json.dumps(res))
if not keep:
    shutil.rmtree(sim_path, ignore_errors=True)

# comparison on the excited port
zf = np.array(fs['ports'][0]['z_re']) + 1j * np.array(fs['ports'][0]['z_im'])
zo = np.array(res['ports'][exc]['z_re']) + 1j * np.array(res['ports'][exc]['z_im'])
band = fs['ports'][0]['l_fit']['band']
sel = (freqs >= band[0]) & (freqs <= band[1])
w = 2 * np.pi * freqs
lf = np.mean(zf.imag[sel] / w[sel])
lo_ = np.mean(zo.imag[sel] / w[sel])
print('port %s: L(fieldsolve) = %.4f nH, L(openEMS) = %.4f nH, diff %+.2f %%' % (fs['ports'][0]['name'], lf * 1e9, lo_ * 1e9, 100 * (lf / lo_ - 1)))
rel = np.abs(zf - zo) / np.maximum(np.abs(zo), 1e-9)
for f in [freqs[len(freqs) // 10], freqs[len(freqs) // 3], freqs[2 * len(freqs) // 3], freqs[-1]]:
    i = int(np.argmin(np.abs(freqs - f)))
    print('  f = %8.3f MHz  Z fs = %10.4f%+10.4fj  openEMS = %10.4f%+10.4fj  |dZ|/|Z| = %.3f' % (freqs[i] * 1e-6, zf[i].real, zf[i].imag, zo[i].real, zo[i].imag, rel[i]))
print('  median |dZ|/|Z| over the band: %.4f' % np.median(rel))
