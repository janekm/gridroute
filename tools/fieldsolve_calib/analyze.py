"""Fast-mode calibration: screening / refined runs against full runs of the same candidates.

Expects runs/<cand>/result.json (the spec as given, `fieldsolve run`) and runs/<cand>_screen/, runs/<cand>_refined/
(`--fast`, `--fidelity refined`) as written by run_all.sh. Prints the table of docs/fieldsolve-fast-mode.md and
writes calib.json.

    python3 analyze.py [--replay runs/base_long]   (replay the objective stop rule on a long full-mesh history)
"""
import json, math, sys
from itertools import combinations

C = ['base', 'trim1', 'trim2', 'ext1', 'ext15', 'feedL05', 'feedL10', 'feedR03', 'gnd+05', 'gnd+10', 'gnd-05', 'trim1_feedL05']
LEVELS = ['screen', 'refined']


def load(path):
    try:
        return json.load(open(path))
    except FileNotFoundError:
        return None


def spearman(a, b):
    def ranks(x):
        r = [0] * len(x)
        for k, i in enumerate(sorted(range(len(x)), key=lambda i: x[i])):
            r[i] = k
        return r
    ra, rb, n = ranks(a), ranks(b), len(a)
    return 1 - 6 * sum((p - q) ** 2 for p, q in zip(ra, rb)) / (n * (n * n - 1))


def pairs(a, b, tie=0.0):
    """Pairs ordered the same way, among pairs whose reference difference exceeds tie."""
    ok = tot = 0
    for i, j in combinations(range(len(a)), 2):
        if abs(a[i] - a[j]) > tie:
            tot += 1
            ok += (a[i] - a[j]) * (b[i] - b[j]) > 0
    return ok, tot


def cost(r):
    return r['cells'] * sum(x['steps'] for x in r['runs'])


rows = []
for c in C:
    v = load(f'runs/{c}/result.json')
    if not v:
        continue
    row = {'cand': c, 'ref': v['objective']['worst_db'], 'ref_per_f': [p['total_db'] for p in v['objective']['per_f']],
           'ref_steps': v['runs'][0]['steps'], 'ref_s': v['runs'][0]['seconds'], 'ref_cost': cost(v)}
    for t in LEVELS:
        r = load(f'runs/{c}_{t}/result.json')
        if r:
            pf = [p['total_db'] for p in r['objective']['per_f']]
            row[t] = {'worst': r['objective']['worst_db'], 'err': r['objective']['worst_db'] - row['ref'],
                      'per_f_err': max(abs(x - y) for x, y in zip(pf, row['ref_per_f'])), 'steps': r['runs'][0]['steps'],
                      's': r['runs'][0]['seconds'], 'cost_ratio': row['ref_cost'] / cost(r), 'stop': r['runs'][0]['stop'],
                      'status': r['fidelity']['status'], 'cells': r['cells'], 'dt': r['dt']}
    rows.append(row)

print(f"{'candidate':<15}{'full':>8}{'steps':>7}" + ''.join(f" | {t:>8}{'err':>7}{'f err':>7}{'steps':>7}{'cost/':>6}" for t in LEVELS))
for r in rows:
    line = f"{r['cand']:<15}{r['ref']:>8.3f}{r['ref_steps']:>7}"
    for t in LEVELS:
        x = r.get(t)
        line += f" | {x['worst']:>8.3f}{x['err']:>+7.3f}{x['per_f_err']:>7.3f}{x['steps']:>7}{x['cost_ratio']:>6.2f}" if x else ' | ' + ' ' * 35
    print(line)

summary = {}
for t in LEVELS:
    xs = [r for r in rows if t in r]
    if len(xs) < 2:
        continue
    a, b = [r['ref'] for r in xs], [r[t]['worst'] for r in xs]
    errs = [r[t]['err'] for r in xs]
    mean = sum(errs) / len(errs)
    sd = math.sqrt(sum((e - mean) ** 2 for e in errs) / (len(errs) - 1))
    ok, tot = pairs(a, b)
    ok2, tot2 = pairs(a, b, 2 * sd)
    costs = sorted(r[t]['cost_ratio'] for r in xs)
    walls = sorted(r['ref_s'] / r[t]['s'] for r in xs)
    s = summary[t] = {'n': len(xs), 'mean_err_db': mean, 'sd_err_db': sd, 'max_abs_err_db': max(map(abs, errs)),
                      'max_abs_per_f_err_db': max(r[t]['per_f_err'] for r in xs), 'spearman': spearman(a, b),
                      'pairs_same_order': [ok, tot], 'pairs_same_order_beyond_2sd': [ok2, tot2],
                      'best_same': max(range(len(a)), key=a.__getitem__) == max(range(len(b)), key=b.__getitem__),
                      'cells_steps_ratio_median': costs[len(costs) // 2], 'cells_steps_ratio_range': [costs[0], costs[-1]],
                      'wall_ratio_median': walls[len(walls) // 2], 'cells_ratio': xs[0]['ref_cost'] / xs[0]['ref_steps'] / xs[0][t]['cells']}
    print(f"\n{t}: n {s['n']}, error mean {mean:+.3f} sd {sd:.3f} max |{s['max_abs_err_db']:.3f}| dB (per f {s['max_abs_per_f_err_db']:.3f}), "
          f"Spearman {s['spearman']:.3f}, pairs {ok}/{tot}, beyond 2 sd {ok2}/{tot2}, best same {s['best_same']}, "
          f"cells x steps /{s['cells_steps_ratio_median']:.2f} ({costs[0]:.2f}-{costs[-1]:.2f}), wall /{s['wall_ratio_median']:.2f}, cells /{s['cells_ratio']:.2f}")

if '--replay' in sys.argv:
    # the objective stop of fdtd.rs (spread over the window <= tol, after the minimum settling time) replayed on
    # a long full-mesh run, against the centre of the radiated power's late beat
    v = load(sys.argv[sys.argv.index('--replay') + 1] + '/result.json')
    run, dt, f = v['runs'][0], v['dt'], v['objective']['band'][0]
    comps, hist = run['objective_components'], run['objective_history']
    rad = [i for i, c in enumerate(comps) if c.startswith('radiated')]
    worst = lambda x: min(x[i] for i in rad)
    late = [worst(x) for s, x in hist if s >= run['steps'] - 50000]
    centre = (max(late) + min(late)) / 2
    std = run.get('standard_rule_step')
    print(f"\nreplay on {run['steps']} steps: late beat centre {centre:.3f} dB, half amplitude {(max(late) - min(late)) / 2:.3f} dB, "
          f"standard rule at {std}")
    rep = {}
    for w in [2, 3, 4, 6, 8]:
        for tol in [0.05, 0.1, 0.2, 0.3]:
            win, ms = max(3000, math.ceil(w / f / dt)), 3600 + math.ceil(8 / f / dt)
            for s, x in hist:
                if s < max(ms, win):
                    continue
                seg = [y for t, y in hist if s - win <= t <= s]
                if max(max(q[c] for q in seg) - min(q[c] for q in seg) for c in range(len(comps))) <= tol:
                    rep[f'{w}/{tol}'] = [s, worst(x) - centre]
                    print(f'  window {w} periods, tol {tol} dB: stop {s} ({s / std:.2f} x standard), error {worst(x) - centre:+.3f} dB')
                    break
    summary['replay'] = rep
json.dump({'rows': rows, 'summary': summary}, open('calib.json', 'w'), indent=1)
