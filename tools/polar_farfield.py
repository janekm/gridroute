"""Polar cuts of a far-field pattern written by `fieldsolve cosim` with "farfield": true.

    python3 tools/polar_farfield.py PATTERN.json [PATTERN.json ...] -o OUT.png [--title T] [--labels A,B]
        [--range 30] [--accepted]

Three cuts in the solver frame (Gerber orientation, z out of the top copper): the board plane (theta 90,
phi 0-360), the x-z plane (phi 0 / 180) and the y-z plane (phi 90 / 270). Each panel shows the total and
the E_theta / E_phi parts. Levels are EIRP per watt incident at the driven port, i.e. realized gain
in dBi for a 50 ohm source; with --accepted they are divided by the accepted power (gain in dBi).
Several patterns (e.g. frequencies) are overlaid, total only, when more than one is given.
"""
import argparse
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('patterns', nargs='+')
ap.add_argument('-o', '--out', required=True)
ap.add_argument('--title', default='')
ap.add_argument('--labels')
ap.add_argument('--range', type=float, default=30.0)
ap.add_argument('--accepted', action='store_true', help='normalise to accepted power (gain, not realized gain)')
a = ap.parse_args()

pats = [json.load(open(p)) for p in a.patterns]
labels = a.labels.split(',') if a.labels else [f"{p['f'] / 1e9:.3f} GHz" for p in pats]


def grids(p):
    th = np.array(p['theta_deg']); ph = np.array(p['phi_deg'])
    sh = (len(th), len(ph))
    off = -p.get('p_in_db', 0.0) if a.accepted else 0.0
    g = {k: np.array(p[k]).reshape(sh) + off for k in ('eirp_per_w_db', 'eirp_theta_db', 'eirp_phi_db') if k in p}
    return th, ph, g


def cut(th, ph, G, kind):
    """angle (rad, polar-plot convention) and level along one cut"""
    if kind == 'xy':
        it = int(np.argmin(abs(th - 90)))
        ang = np.radians(ph); val = G[it, :]
        return np.append(ang, ang[0]), np.append(val, val[0])
    p0, p1 = (0, 180) if kind == 'xz' else (90, 270)
    i0 = int(np.argmin(abs(ph - p0))); i1 = int(np.argmin(abs(ph - p1)))
    # theta from +z: going over the top from the phi0 side (0..180) and back down on the phi1 side
    ang = np.concatenate([np.radians(th), 2 * np.pi - np.radians(th[::-1])])
    val = np.concatenate([G[:, i0], G[::-1, i1]])
    return ang, val


top = max(np.max(grids(p)[2]['eirp_per_w_db']) for p in pats)
top = np.ceil(top / 5) * 5
lo = top - a.range
fig, axs = plt.subplots(1, 3, subplot_kw={'projection': 'polar'}, figsize=(16, 6))
names = {'xy': 'board plane (θ = 90°), φ from +x',
         'xz': 'x-z plane (φ = 0° right, 180° left), θ from +z',
         'yz': 'y-z plane (φ = 90° right, 270° left), θ from +z'}
cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
for ax, kind in zip(axs, ('xy', 'xz', 'yz')):
    for n, (p, lab) in enumerate(zip(pats, labels)):
        th, ph, g = grids(p)
        keys = [('eirp_per_w_db', 'total', '-')]
        if len(pats) == 1:
            keys += [('eirp_theta_db', 'E_θ', '--'), ('eirp_phi_db', 'E_φ', ':')]
        for k, nm, ls in keys:
            if k not in g:
                continue
            ang, val = cut(th, ph, g[k], kind)
            ax.plot(ang, np.clip(val, lo, None), ls, color=cols[n % len(cols)] if len(pats) > 1 else {'total': 'k', 'E_θ': cols[0], 'E_φ': cols[1]}[nm],
                    lw=1.8 if nm == 'total' else 1.2, label=f'{lab} {nm}' if len(pats) > 1 or nm != 'total' else f'{lab} total')
    ax.set_rlim(lo, top)
    ax.set_rticks(np.arange(lo, top + 0.1, 10 if a.range > 25 else 5))
    if kind == 'xy':
        ax.set_theta_zero_location('E'); ax.set_theta_direction(1)
    else:
        ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
    ax.set_title(names[kind], fontsize=10, pad=14)
    ax.grid(True, alpha=0.5)
axs[0].legend(loc='upper left', bbox_to_anchor=(-0.25, 1.18), fontsize=8)
unit = 'gain, dBi (per W accepted)' if a.accepted else 'realized gain, dBi (EIRP per W incident at the driven port)'
fig.suptitle(f'{a.title}\n{unit}', fontsize=12)
fig.savefig(a.out, dpi=110, bbox_inches='tight')
print('wrote', a.out, 'top', top)
