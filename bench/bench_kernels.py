"""Micro-benchmark: numpy reference vs gridroute CPU vs Metal for the clearance dilation, by window size.

    python3 bench/bench_kernels.py [--board 170x111] [--grid 0.05]

The source grid is a synthetic board: random pads/tracks of many nets on four layers at roughly the density of a
routed board. Every timed result is also checked against the numpy reference.
"""
import argparse
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
import gridroute as gr  # noqa: E402


def blocked_np(occ, L, nid, win, r, G):
    """The numpy dilation gridroute replaces (reference implementation)."""
    ny, nx = occ.shape[1:]
    i0, j0, i1, j1 = win
    rc = int(math.ceil(r / G))
    pi0, pj0 = max(i0 - rc, 0), max(j0 - rc, 0)
    pi1, pj1 = min(i1 + rc, ny - 1), min(j1 + rc, nx - 1)
    src = occ[L, pi0:pi1 + 1, pj0:pj1 + 1]
    other = (src != 0) & (src != nid)
    H, W = other.shape
    cs = np.zeros((H, W + 1), dtype=np.int32)
    np.cumsum(other, axis=1, out=cs[:, 1:])
    dil = np.zeros_like(other)
    rr = (r / G) ** 2
    jj = np.arange(W)
    for di in range(-rc, rc + 1):
        if di * di >= rr - 1e-9:
            continue
        w = int(math.floor(math.sqrt(rr - di * di) - 1e-9))
        lo, hi = np.clip(jj - w, 0, W), np.clip(jj + w + 1, 0, W)
        hd = (cs[:, hi] - cs[:, lo]) > 0
        if di >= 0:
            dil[:H - di] |= hd[di:]
        else:
            dil[-di:] |= hd[:H + di]
    return dil[i0 - pi0:i1 - pi0 + 1, j0 - pj0:j1 - pj0 + 1]


def synthetic(ny, nx, nl=4, seed=1):
    rng = np.random.default_rng(seed)
    occ = np.zeros((nl, ny, nx), dtype=np.int16)
    for L in range(nl):
        for _ in range(ny * nx // 4000):
            y, x = rng.integers(0, ny - 40), rng.integers(0, nx - 40)
            h, w = rng.integers(2, 30), rng.integers(2, 30)
            occ[L, y:y + h, x:x + w] = rng.integers(1, 400)
    return occ


def best(f, n):
    ts = []
    for _ in range(n):
        t = time.perf_counter()
        r = f()
        ts.append(time.perf_counter() - t)
    return min(ts), r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--board', default='170x111')
    ap.add_argument('--grid', type=float, default=0.05)
    ap.add_argument('--radius', type=float, default=0.45, help='reach in mm (via 0.25 + clearance 0.2)')
    a = ap.parse_args()
    bw, bh = map(float, a.board.split('x'))
    G = a.grid
    ny, nx = int(math.ceil(bh / G)) + 1, int(math.ceil(bw / G)) + 1
    occ = synthetic(ny, nx)
    print('grid %d x %d x 4, reach %.2f mm (%d cells), Metal: %s' % (nx, ny, a.radius, math.ceil(a.radius / G),
                                                                     gr.device_name() or 'no'))
    print('%-14s %10s %10s %10s %10s %8s %8s' % ('window', 'numpy 4L', 'cpu 4L', 'metal 4L', 'cpu OR', 'cpu x', 'metal x'))
    span = gr.disc_span(a.radius / G)
    for mm in (5, 10, 20, 40, 80, 170):
        h, w = min(int(mm / G), ny - 1), min(int(mm / G), nx - 1)
        i0, j0 = (ny - h) // 2, (nx - w) // 2
        win = (i0, j0, i0 + h - 1, j0 + w - 1)
        n = 3 if mm >= 80 else 10
        t_np, ref = best(lambda: np.stack([blocked_np(occ, L, 7, win, a.radius, G) for L in range(4)]), n)
        t_cpu, r_cpu = best(lambda: gr.dilate(occ, [0, 1, 2, 3], win, span, excl=(7,), backend='cpu'), n * 3)
        assert np.array_equal(ref, r_cpu)
        t_or, r_or = best(lambda: gr.dilate(occ, [0, 1, 2, 3], win, span, excl=(7,), reduce_or=True, backend='cpu'), n * 3)
        assert np.array_equal(ref.any(axis=0), r_or)
        if gr.metal_available():
            t_mt, r_mt = best(lambda: gr.dilate(occ, [0, 1, 2, 3], win, span, excl=(7,), backend='metal'), n * 3)
            assert np.array_equal(ref, r_mt)
        else:
            t_mt = float('nan')
        print('%-14s %9.2fms %9.2fms %9.2fms %9.2fms %7.0fx %7.0fx' % (
            '%dx%d mm' % (min(mm, bw), min(mm, bh)), t_np * 1e3, t_cpu * 1e3, t_mt * 1e3, t_or * 1e3,
            t_np / t_cpu, t_np / t_mt))


if __name__ == '__main__':
    main()
