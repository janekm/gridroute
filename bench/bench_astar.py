"""Replay a corpus of recorded A* instances: optimised search vs reference, checking identical paths.

    python3 bench/bench_astar.py CORPUS_DIR [--limit N] [--no-ref]

Each .npz holds one instance as recorded from a router run (see README, "Recording a corpus").
"""
import argparse
import glob
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
import gridroute as gr  # noqa: E402


def load(f):
    z = np.load(f)
    shape = tuple(z['shape'])
    n = int(np.prod(shape))
    unpack = lambda a, m, s: np.unpackbits(a)[:m].astype(bool).reshape(s)
    return dict(blk=unpack(z['blk'], n, shape), vok=unpack(z['vok'], shape[1] * shape[2], shape[1:]),
                tgt=unpack(z['tgt'], n, shape), src=z['src'], mcost=z['mcost'], lay_ok=z['lay_ok'].astype(bool),
                vcost=float(z['par'][0]), turn=float(z['par'][1]), hmul=float(z['par'][2]),
                max_exp=int(z['max_exp']), path=z['path'], t=float(z['t']))


def run(inst, **kw):
    t = time.perf_counter()
    e0 = gr.stats['astar_exp']
    p = gr.astar(inst['blk'], inst['vok'], inst['src'], inst['tgt'], inst['mcost'], inst['lay_ok'], inst['vcost'],
                 inst['turn'], inst['hmul'], inst['max_exp'], **kw)
    return p, time.perf_counter() - t, gr.stats['astar_exp'] - e0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('corpus')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--no-ref', action='store_true')
    a = ap.parse_args()
    files = sorted(glob.glob(os.path.join(a.corpus, '*.npz')))
    if a.limit:
        files = files[:a.limit]
    tot = {'rec': 0.0, 'ref': 0.0, 'new': 0.0}
    nfail = 0
    exps = 0
    for f in files:
        inst = load(f)
        p_new, t_new, e = run(inst)
        exps += e
        rec = inst['path'] if len(inst['path']) else None
        same = (p_new is None and rec is None) or (p_new is not None and rec is not None and np.array_equal(p_new, rec))
        if not same:
            nfail += 1
            print('MISMATCH vs recorded path:', os.path.basename(f))
        tot['new'] += t_new
        tot['rec'] += inst['t']
        if not a.no_ref:
            p_ref, t_ref, _ = run(inst, reference=True)
            tot['ref'] += t_ref
            if not ((p_ref is None and p_new is None) or (p_ref is not None and p_new is not None and np.array_equal(p_ref, p_new))):
                nfail += 1
                print('MISMATCH vs reference:', os.path.basename(f))
    print('%d instances, %d mismatches, %.0fM expansions' % (len(files), nfail, exps / 1e6))
    print('recorded (in run) %.2fs   reference %.2fs   optimised %.2fs   -> %.2fx vs reference' % (
        tot['rec'], tot['ref'], tot['new'], tot['ref'] / tot['new'] if tot['ref'] else float('nan')))
    sys.exit(1 if nfail else 0)


if __name__ == '__main__':
    main()
