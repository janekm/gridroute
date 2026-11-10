"""Calibrate a chip-antenna stand-in on an EM two-port (FEED, LOAD; tools/jti_chip_eval.py): a series C + R on
LOAD (the stand-in's resonance and internal loss), an ideal two-element L-match at 2.44 GHz on FEED, the matched
S11 at the band edges and, through `fieldsolve cosim` with the far field, the radiated share.

    FIELDSOLVE=... python3.12 chip_standin_fit.py RUN_DIR C1,C2,.. R1,R2,.. [tag]

Calibrated on Johanson's evaluation board: C 0.69 pF, R 2.3 ohm (22.5 %, S11 -6.2 dB at the band edges against
the measured 21 % and -5.7 dB). See docs/ebus-fpc-antenna-findings.md.
"""

import json, math, sys
import numpy as np

def read_s2p(path):
    rows = [[float(v) for v in l.split()] for l in open(path) if l.strip() and l[0] not in '!#']
    a = np.array(rows); f = a[:, 0]
    s = np.zeros((len(f), 2, 2), complex)
    s[:, 0, 0] = a[:, 1] + 1j * a[:, 2]; s[:, 1, 0] = a[:, 3] + 1j * a[:, 4]
    s[:, 0, 1] = a[:, 5] + 1j * a[:, 6]; s[:, 1, 1] = a[:, 7] + 1j * a[:, 8]
    return f, s

def s2z(s, z0=50.0):
    I = np.eye(2)
    return np.array([z0 * (I + m) @ np.linalg.inv(I - m) for m in s])

def zin(Z, f, R, L, C):
    w = 2 * np.pi * f
    zc = R + 1j * w * L + 1 / (1j * w * C)
    return Z[:, 0, 0] - Z[:, 0, 1] * Z[:, 1, 0] / (Z[:, 1, 1] + zc), zc

def lmatch(zl, f0, z0=50.0):
    """Ideal two-element L-matches for zl at f0: list of (topology, X_series, B_shunt) solutions."""
    R, X = zl.real, zl.imag
    out = []
    if R < z0:  # series element next to the load, shunt at the source
        q = math.sqrt(R * (z0 - R))
        for sgn in (1, -1):
            Xs = sgn * q - X
            B = sgn * math.sqrt((z0 - R) / R) / z0
            out.append(('series-load', Xs, B))
    if R > 0:  # shunt next to the load, series at the source
        G, Bl = (1 / zl).real, (1 / zl).imag
        if G < 1 / z0:
            for sgn in (1, -1):
                Bs = sgn * math.sqrt(G * (1 / z0 - G)) - Bl
                Xs2 = sgn * math.sqrt((1 / z0 - G) / G) * 1 / (1 / z0) / 1 * (1 / 1)  # placeholder, recomputed below
                out.append(('shunt-load', None, Bs))
    return out

def matched_gamma(zl, f, f0, z0=50.0):
    """|Gamma(f)| for the best ideal L-match designed at f0 (elements scale as ideal L or C)."""
    i0 = np.argmin(abs(f - f0)); zl0 = zl[i0]
    best = None
    w = 2 * np.pi * f; w0 = 2 * np.pi * f0
    def react(X0):  # series reactance X0 at f0 as an ideal L (X0>0) or C
        return X0 * w / w0 if X0 > 0 else X0 * w0 / w
    def susc(B0):
        return B0 * w / w0 if B0 > 0 else B0 * w0 / w
    cands = []
    R, X = zl0.real, zl0.imag
    if R < z0:
        q = math.sqrt(R * (z0 - R))
        for sgn in (1, -1):
            Xs = sgn * q - X; Bp = sgn * math.sqrt((z0 - R) / R) / z0
            z = zl + 1j * react(Xs); y = 1 / z + 1j * susc(Bp); cands.append(abs((1 / y - z0) / (1 / y + z0)))
    yl0 = 1 / zl0; G, B = yl0.real, yl0.imag
    if 0 < G < 1 / z0:
        for sgn in (1, -1):
            Bp = sgn * math.sqrt(G * (1 / z0 - G)) - B
            y1 = 1 / zl0 + 1j * Bp; Xs = -(1 / y1).imag
            y = 1 / zl + 1j * susc(Bp); z = 1 / y + 1j * react(Xs); cands.append(abs((z - z0) / (z + z0)))
    return min(cands, key=lambda g: g[i0] + max(g[np.argmin(abs(f - 2.40e9))], g[np.argmin(abs(f - 2.48e9))]))


import os, subprocess
FS = os.environ.get('FIELDSOLVE', 'fieldsolve')

def zload(f, C, R, L=0.0):
    w = 2 * np.pi * f
    return R + 1j * w * L + 1 / (1j * w * C)

def scan(run, Cs, Rs, f0=2.44e9, tag='scan'):
    f, s = read_s2p(f'{run}/sparams.s2p'); Z = s2z(s); i0 = np.argmin(abs(f - f0))
    var, meta = {}, {}
    for C in Cs:
        for R in Rs:
            zl = zload(f, C * 1e-12, R)
            zi = Z[:, 0, 0] - Z[:, 0, 1] * Z[:, 1, 0] / (Z[:, 1, 1] + zl)
            g = matched_gamma(zi, f, f0)
            e = [20 * math.log10(g[np.argmin(abs(f - x))]) for x in (2.40e9, 2.44e9, 2.48e9)]
            n = f'C{C}_R{R}'
            var[n] = {'LOAD': f'C:{C}p,esr={R}'}
            meta[n] = (C, R, zi[i0], e)
    cfg = {'em': f'{run}/sparams.s2p', 'out': f'{run}_{tag}', 'keep': ['FEED'], 'variants': var, 'farfield': True, 'report': [f0]}
    json.dump(cfg, open(f'{run}_{tag}.json', 'w'), indent=1)
    subprocess.run([FS, 'cosim', f'{run}_{tag}.json'], capture_output=True)
    c = json.load(open(f'{run}_{tag}/cosim.json'))
    res = {}
    for n, (C, R, z, e) in meta.items():
        v = c[n] if n in c else c['variants'][n]
        eff = {round(r['f'] / 1e7): 10 ** ((r['p_rad_db'] - r['p_in_db']) / 10) for r in v['farfield']}
        res[n] = (C, R, z, e, eff)
        print(f'C {C:6.3f} pF R {R:5.2f}: Zin(2.44) {z.real:8.2f}{z.imag:+9.2f}j  matched S11 2.40/2.44/2.48 {e[0]:6.1f}/{e[1]:6.1f}/{e[2]:6.1f} dB  eff 2.40/2.44/2.48 {eff.get(240,0)*100:5.1f}/{eff.get(244,0)*100:5.1f}/{eff.get(248,0)*100:5.1f} %')
    return res

if __name__ == '__main__':
    run = sys.argv[1]
    Cs = [float(x) for x in sys.argv[2].split(',')]; Rs = [float(x) for x in sys.argv[3].split(',')]
    scan(run, Cs, Rs, tag=sys.argv[4] if len(sys.argv) > 4 else 'scan')
