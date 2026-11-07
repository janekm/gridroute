"""Retune a matching ladder by coordinate descent over real part values, scoring each candidate with
`fieldsolve cosim` on an existing EM result. Set up for the MTH R7A PCB PIFA (docs/mth-antenna-findings.md):
edit EM, PORT, CHOICES and BASE for another board.

Run it from a working directory holding the EM run (EM below), Nordic's target network
(nordic_ref_network.s2p) and Murata's S-parameter libraries unpacked as murata_sparams/lib_grm and
murata_sparams/lib_lqg. FIELDSOLVE points at the binary (default: fieldsolve on PATH).

Score per frequency: (1 - |G_t|^2) * P_rad / P_in, where G_t is the power-wave reflection against the chip's
design load (Nordic's network, 35 - j35 ohm) and P_rad / P_in the share of the power entering the chip pin
that is radiated (ladder + antenna losses). The objective is the worst of 2.40 / 2.44 / 2.48 GHz, in dB.

    python ladder_retune.py [rounds]
"""
import glob, json, math, os, random, subprocess, sys

FS = os.environ.get('FIELDSOLVE', 'fieldsolve')
EM = 'mth7_board/sparams.s7p'
PORT = 'CHIP'
GRM = 'murata_sparams/lib_grm'; LQG = 'murata_sparams/lib_lqg'


def cap_file(code):
    for tol in 'BCWDJ':
        f = f'{GRM}/GRM1555C1H{code}{tol}A01.s2p'
        if os.path.exists(f):
            return f
    return None


def ind_file(code):
    for tol in 'BCSJGH':
        f = f'{LQG}/LQG15HS{code}{tol}02.s2p'
        if os.path.exists(f):
            return f
    return None


CAPS = {}
for v in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.6, 6.8, 8.2]:
    code = f'R{int(round(v * 10))}0' if v < 1 else f'{int(v)}R{int(round(v * 10)) % 10}'
    f = cap_file(code)
    if f:
        CAPS[f'{v:g}p'] = 's2p:' + f
INDS = {}
for v in [1.0, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2]:
    code = f'{int(v)}N{int(round(v * 10)) % 10}'
    f = ind_file(code)
    if f:
        INDS[f'{v:g}n'] = 's2p:' + f
small_caps = {k: v for k, v in CAPS.items() if float(k[:-1]) <= 2.0}
CHOICES = {
    'C59': {'open': 'open', **CAPS},
    'L4': {'short': 'short', **INDS},
    'C60': {'open': 'open', **CAPS},
    'L5': {'short': 'short', **INDS},
    'R70': {'short': 'short', **{k: v for k, v in INDS.items() if float(k[:-1]) <= 4.7}, **{k + ' ser': v for k, v in CAPS.items() if float(k[:-1]) >= 0.5}},
    'C20': {'open': 'open', **small_caps},
}
BASE = {'C59': '1p', 'L4': '4.7n', 'C60': '1.2p', 'L5': '2.2n', 'R70': 'short', 'C20': 'open'}
FREQS = (2.40e9, 2.44e9, 2.48e9)
cache = {}
calls = [0]


def evaluate(cands):
    """cands: list of {pos: label}; returns list of (score_db, per_f list, harmonics)"""
    todo = [c for c in cands if tuple(sorted(c.items())) not in cache]
    for start in range(0, len(todo), 120):
        chunk = todo[start:start + 120]
        cfg = {"em": EM, "out": "retune_tmp", "keep": [PORT], "target": "nordic_ref_network.s2p",
               "report": list(FREQS), "plot_ideal": False, "farfield": True,
               "variants": {f"v{i}": {p: CHOICES[p][lab] for p, lab in c.items()} for i, c in enumerate(chunk)}}
        json.dump(cfg, open('retune_tmp.json', 'w'))
        r = subprocess.run([FS, 'cosim', 'retune_tmp.json'], capture_output=True, text=True)
        if r.returncode:
            sys.exit(r.stderr[-2000:])
        calls[0] += 1
        res = json.load(open('retune_tmp/cosim.json'))['variants']
        for i, c in enumerate(chunk):
            v = res[f'v{i}']
            g = {round(row['f'] / 1e7): row.get('gamma_vs_target') for row in v['report']}
            ff = {round(row['f'] / 1e7): row for row in v['farfield']}
            per = []
            for f in FREQS:
                k = round(f / 1e7)
                gk = min(g, key=lambda x: abs(x - k))
                eta = 10 ** ((ff[k]['p_rad_db'] - ff[k]['p_in_db']) / 10)
                mis = 1 - g[gk] ** 2
                per.append((10 * math.log10(max(mis * eta, 1e-12)), g[gk], 10 * math.log10(eta)))
            harm = {round(f / 1e8) / 10: ff[k]['p_rad_db'] for k, f in [(488, 4.88e9), (732, 7.32e9)] if k in ff}
            cache[tuple(sorted(c.items()))] = (min(p[0] for p in per), per, harm)
        for x in glob.glob('retune_tmp/farfield_*.json'):
            os.remove(x)
    return [cache[tuple(sorted(c.items()))] for c in cands]


def descend(start, rounds):
    cur = dict(start)
    best = evaluate([cur])[0]
    for rd in range(rounds):
        improved = False
        for pos in random.sample(list(CHOICES), len(CHOICES)):
            cands = [dict(cur, **{pos: lab}) for lab in CHOICES[pos]]
            scores = evaluate(cands)
            k = max(range(len(cands)), key=lambda i: scores[i][0])
            if scores[k][0] > best[0] + 1e-4:
                cur, best, improved = cands[k], scores[k], True
        if not improved:
            break
    return cur, best


def show(tag, c, s):
    per = ' | '.join(f'{f/1e9:.2f}: {d:+.2f} dB (|G| {g:.2f}, rad/acc {e:+.2f} dB)' for f, (d, g, e) in zip(FREQS, s[1]))
    print(f'{tag:10s} worst {s[0]:+.2f} dB  ' + ' '.join(f'{p}={c[p]}' for p in CHOICES) + f'\n           {per}\n           radiated H2/H3 per W into the chip pin (50 ohm): {s[2]}')


if __name__ == '__main__':
    random.seed(7)
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    b = evaluate([BASE])[0]
    show('as built', BASE, b)
    results = []
    starts = [BASE] + [{p: random.choice(list(CHOICES[p])) for p in CHOICES} for _ in range(5)]
    for i, st in enumerate(starts):
        c, s = descend(st, rounds)
        results.append((s[0], c, s))
        show(f'start {i}', c, s)
    results.sort(key=lambda x: -x[0])
    json.dump([{'score_db': r[0], 'parts': r[1], 'per_f': r[2][1], 'harm': r[2][2]} for r in results], open('ladder_retune_results.json', 'w'), indent=1)
    print(f'{calls[0]} cosim calls, {len(cache)} candidates')
