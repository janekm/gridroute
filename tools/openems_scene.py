#!/usr/bin/env python3
"""Replay a fieldsolve scene in openEMS and compare port impedances and S-parameters.

    ~/opt/openEMS/venv/bin/python tools/openems_scene.py RUN_DIR [--ports all|0,2] [--polygons] [--keep]
                                                         [--threads N] [--out NAME]

RUN_DIR is a `fieldsolve run` output directory (scene.json, result.json, sparams.sNp). The scene is rebuilt
with exactly fieldsolve's mesh lines (no smoothing), the same lumped ports/elements on the same boxes, the
same Gaussian excitation, PML thickness and end criterion. The structure comes from the scene's resolved
raster by default: every perfect-conductor edge fieldsolve used, given to openEMS as wires along mesh lines
(which openEMS maps to exactly those edges), and the dielectric as boxes of whole cells. Both codes then
solve the same discrete structure, so differences are the solvers', not the rasterisers'. `--polygons`
uses the copper polygons, outline and barrels instead (each code rasterises them its own way).

`--ports` selects the excitations (default: the scene's excited port). Each excitation's column is saved
in RUN_DIR/NAME/col_<port>.json and reused when it exists, so a long multi-port run can be resumed. Once
every column exists, RUN_DIR/NAME/sparams.sNp is written and compared with fieldsolve's.

openEMS is GPLv3 and used here only as an external reference; nothing from it is linked into fieldsolve.
"""
import argparse
import json
import math
import pathlib
import shutil
import sys
import time

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import EPS0

ap = argparse.ArgumentParser()
ap.add_argument('run')
ap.add_argument('--ports', default=None)
ap.add_argument('--polygons', action='store_true')
ap.add_argument('--keep', action='store_true')
ap.add_argument('--threads', type=int, default=0)
ap.add_argument('--out', default='openems')
args = ap.parse_args()

run = pathlib.Path(args.run).resolve()
sc = json.loads((run / 'scene.json').read_text())
fs = json.loads((run / 'result.json').read_text())
out = run / args.out
out.mkdir(parents=True, exist_ok=True)
use_raster = 'raster' in sc and not args.polygons
gx, gy, gz = (np.array(sc['grid'][a]) for a in 'xyz')
lines = [gx, gy, gz]
freqs = np.array(sc.get('freqs') or fs['ports'][0]['freqs'])
nports = sum(1 for l in sc['lumped'] if l['port'])
if args.ports in (None, ''):
    excite = [sc['excite'] or 0]
elif args.ports == 'all':
    excite = list(range(nports))
else:
    excite = [int(v) for v in args.ports.split(',')]


def snap(ls, v):
    return int(np.argmin(np.abs(ls - v)))


def build(exc):
    fdtd = openEMS(NrTS=int(sc['max_steps']), EndCriteria=10 ** (sc['end_db'] / 10.0))
    fdtd.SetGaussExcite(0.0, sc['fmax'])
    fdtd.SetBoundaryCond([('PML_%d' % n) if n > 0 else 'PEC' for n in sc['pml']])
    csx = ContinuousStructure()
    fdtd.SetCSX(csx)
    mesh = csx.GetGrid()
    mesh.SetDeltaUnit(sc['unit'])
    for a, ls in zip('xyz', lines):
        mesh.SetLines(a, ls)
    w_loss = 2 * math.pi * sc['f_loss']
    mats = [csx.AddMaterial('slab%d' % k, epsilon=s['eps'], kappa=w_loss * EPS0 * s['eps'] * s['tand']) for k, s in enumerate(sc['slabs'])]
    metal = csx.AddMetal('copper')
    if use_raster:
        r = sc['raster']
        for s, i0, i1, j, k0, k1 in r['cells']:
            mats[s].AddBox([gx[i0], gy[j], gz[k0]], [gx[i1 + 1], gy[j + 1], gz[k1 + 1]], priority=1)
        for a, a0, a1, b, c in r['metal']:
            ib, ic = (a + 1) % 3, (a + 2) % 3
            p0, p1 = [0.0] * 3, [0.0] * 3
            p0[a], p1[a] = lines[a][a0], lines[a][a1 + 1]
            p0[ib] = p1[ib] = lines[ib][b]
            p0[ic] = p1[ic] = lines[ic][c]
            metal.AddCurve(np.array([p0, p1]).T, priority=500)
    else:
        x0, x1, y0, y1 = gx[0], gx[-1], gy[0], gy[-1]
        for k, s in enumerate(sc['slabs']):
            if sc['outline']:
                # loops are outer first; the outer one is the board, the others cut-outs (even-odd)
                for li, loop in enumerate(sc['outline']):
                    pts = np.array(loop).T
                    if li == 0:
                        mats[k].AddLinPoly(points=pts, norm_dir='z', elevation=s['z0'], length=s['z1'] - s['z0'], priority=1)
                    else:
                        air = csx.AddMaterial('cut%d_%d' % (k, li), epsilon=1.0)
                        air.AddLinPoly(points=pts, norm_dir='z', elevation=s['z0'], length=s['z1'] - s['z0'], priority=2)
            else:
                mats[k].AddBox([x0, y0, s['z0']], [x1, y1, s['z1']], priority=1)
        prio = 10
        for sh in sc['sheets']:
            for p in sh['prims']:
                prio += 1
                if p['dark']:
                    metal.AddPolygon(points=np.array(p['pts']).T, norm_dir='z', elevation=sh['z'], priority=prio)
        # barrels: z lines at the mesh nodes inside the radius (nearest node if none), as fieldsolve does
        for b in sc['barrels']:
            (cx, cy), rr = b['c'], b['r']
            nodes = [(i, j) for i in range(len(gx)) for j in range(len(gy)) if (gx[i] - cx) ** 2 + (gy[j] - cy) ** 2 <= rr * rr]
            if not nodes:
                nodes = [(snap(gx, cx), snap(gy, cy))]
            za, zb = gz[snap(gz, b['z0'])], gz[snap(gz, b['z1'])]
            for i, j in nodes:
                metal.AddCurve(np.array([[gx[i], gx[i]], [gy[j], gy[j]], [za, zb]]), priority=500)
    ports = []
    pnum = 0
    for k, l in enumerate(sc['lumped']):
        lo, hi = list(l['lo']), list(l['hi'])
        start, stop = (lo, hi) if l['plus_at_hi'] else (hi, lo)
        d = 'xyz'[l['axis']]
        if l['port']:
            ports.append(fdtd.AddLumpedPort(pnum + 1, l['r'], start, stop, d, excite=1.0 if pnum == exc else 0, priority=1000))
            pnum += 1
        elif l['r'] == 0 and l['c'] == 0:
            if not use_raster:  # the raster already lists shorted elements as metal
                metal.AddBox(lo, hi, priority=1000)
        else:
            kw = {}
            if l['r'] > 0:
                kw['R'] = l['r']
            if l['c'] > 0:
                kw['C'] = l['c']
            el = csx.AddLumpedElement('el%d' % k, ny=d, caps=True, **kw)
            el.AddBox(lo, hi, priority=1000)
    return fdtd, csx, ports


names = [l['name'] for l in sc['lumped'] if l['port']]
print('openems_scene: %s, %d ports, structure from %s' % (run.name, nports, 'the resolved raster' if use_raster else 'polygons'), flush=True)
for exc in excite:
    col = out / ('col_%d.json' % exc)
    if col.exists():
        print('  port %s: reusing %s' % (names[exc], col.name), flush=True)
        continue
    sim_path = out / ('sim_%d' % exc)
    if sim_path.exists():
        shutil.rmtree(sim_path)
    sim_path.mkdir(parents=True)
    fdtd, csx, ports = build(exc)
    csx.Write2XML(str(out / 'scene.xml'))
    t0 = time.time()
    fdtd.Run(str(sim_path), verbose=0, cleanup=True, numThreads=args.threads)
    for p in ports:
        p.CalcPort(str(sim_path), freqs, ref_impedance=50)
    a = ports[exc].uf_inc
    res = {'excite': exc, 'freqs': freqs.tolist(), 'seconds': time.time() - t0,
           's_re': [(p.uf_ref / a).real.tolist() for p in ports], 's_im': [(p.uf_ref / a).imag.tolist() for p in ports],
           'z_re': (ports[exc].uf_tot / ports[exc].if_tot).real.tolist(), 'z_im': (ports[exc].uf_tot / ports[exc].if_tot).imag.tolist()}
    col.write_text(json.dumps(res))
    print('  port %s excited: %.0f s' % (names[exc], res['seconds']), flush=True)
    if not args.keep:
        shutil.rmtree(sim_path, ignore_errors=True)

# excited-port impedance against fieldsolve's (port 0 of result.json is the scene's first excitation)
e0 = sc['excite'] or 0
if (out / ('col_%d.json' % e0)).exists() and fs['ports'][0].get('z_re'):
    c = json.loads((out / ('col_%d.json' % e0)).read_text())
    zf = np.array(fs['ports'][0]['z_re']) + 1j * np.array(fs['ports'][0]['z_im'])
    zo = np.array(c['z_re']) + 1j * np.array(c['z_im'])
    if len(zf) == len(zo):
        band = fs['ports'][0]['l_fit']['band']
        sel = (freqs >= band[0]) & (freqs <= band[1])
        w = 2 * np.pi * freqs
        lf, lo_ = np.mean(zf.imag[sel] / w[sel]), np.mean(zo.imag[sel] / w[sel])
        print('port %s: L(fieldsolve) = %.4f nH, L(openEMS) = %.4f nH, diff %+.2f %%' % (names[e0], lf * 1e9, lo_ * 1e9, 100 * (lf / lo_ - 1)))
        rel = np.abs(zf - zo) / np.maximum(np.abs(zo), 1e-9)
        for f in [freqs[len(freqs) // 10], freqs[len(freqs) // 3], freqs[2 * len(freqs) // 3], freqs[-1]]:
            i = int(np.argmin(np.abs(freqs - f)))
            print('  f = %8.3f MHz  Z fs = %10.4f%+10.4fj  openEMS = %10.4f%+10.4fj  |dZ|/|Z| = %.3f' % (freqs[i] * 1e-6, zf[i].real, zf[i].imag, zo[i].real, zo[i].imag, rel[i]))
        print('  median |dZ|/|Z| over the band: %.4f' % np.median(rel))

# full S-matrix once every column exists
cols = [out / ('col_%d.json' % j) for j in range(nports)]
if nports > 1 and all(c.exists() for c in cols):
    S = np.zeros((len(freqs), nports, nports), complex)
    for j, cp in enumerate(cols):
        c = json.loads(cp.read_text())
        for i in range(nports):
            S[:, i, j] = np.array(c['s_re'][i]) + 1j * np.array(c['s_im'][i])
    with open(out / ('sparams.s%dp' % nports), 'w') as f:
        f.write('! openEMS replay of %s; ports %s\n# HZ S RI R 50\n' % (run.name, ' '.join(names)))
        for k, fr in enumerate(freqs):
            vals = []
            for i in range(nports):
                for j in range(nports):
                    vals += ['%.9e' % S[k, i, j].real, '%.9e' % S[k, i, j].imag]
            # Touchstone: at most 4 pairs per line for N > 2
            f.write('%.6e' % fr)
            for q in range(0, len(vals), 8):
                f.write((' ' if q == 0 else '\n ') + ' '.join(vals[q:q + 8]))
            f.write('\n')
    asym = np.max(np.abs(S - S.transpose(0, 2, 1)))
    print('openEMS S-matrix written (%d ports); max |S - S^T| = %.4f' % (nports, asym))
