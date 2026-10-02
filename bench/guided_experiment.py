"""Experiment: field-guided A* vs the reference search on a recorded corpus (success, path cost, expansions, time)."""
import ctypes, glob, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridroute as gr
from bench_astar import load
L = gr._LIB
P = ctypes.POINTER
L.gr_astar_guided.argtypes = [ctypes.c_int] * 3 + [P(ctypes.c_uint8)] * 3 + [P(ctypes.c_int32), ctypes.c_int, P(ctypes.c_float), P(ctypes.c_uint8)] + \
    [ctypes.c_float] * 3 + [ctypes.c_int64, P(ctypes.c_int32), ctypes.c_int, P(ctypes.c_int64), P(ctypes.c_double)]
L.gr_astar_guided.restype = ctypes.c_int
DIRS = [(0, 1), (0, -1), (1, 0), (-1, 0), (1, 1), (1, -1), (-1, 1), (-1, -1)]


def cost(path, inst):
    nl, h, w = inst['blk'].shape
    mc = inst['mcost'].reshape(nl, 8)
    c, din = 0.0, 8
    for a, b in zip(path, path[1:]):
        la, ia, ja = a // (h * w), a % (h * w) // w, a % w
        lb, ib, jb = b // (h * w), b % (h * w) // w, b % w
        if la != lb:
            c += inst['vcost']; din = 8; continue
        d = DIRS.index((ib - ia, jb - ja))
        c += mc[la, d] + (inst['turn'] if din != 8 and din != d else 0)
        din = d
    return c


def guided(inst, weight, max_exp):
    nl, h, w = inst['blk'].shape
    u8 = lambda a: np.ascontiguousarray(a).view(np.uint8).ctypes.data_as(P(ctypes.c_uint8))
    src = np.ascontiguousarray(inst['src'], dtype=np.int32)
    mc = np.ascontiguousarray(inst['mcost'], dtype=np.float32).ravel()
    out = np.empty(1 << 22, dtype=np.int32)
    ne, fs = ctypes.c_int64(), ctypes.c_double()
    t = time.perf_counter()
    n = L.gr_astar_guided(nl, h, w, u8(inst['blk']), u8(inst['vok']), u8(inst['tgt']), src.ctypes.data_as(P(ctypes.c_int32)),
                          len(src), mc.ctypes.data_as(P(ctypes.c_float)), u8(inst['lay_ok'].astype(np.uint8)),
                          inst['vcost'], inst['turn'], weight, max_exp, out.ctypes.data_as(P(ctypes.c_int32)), len(out),
                          ctypes.byref(ne), ctypes.byref(fs))
    return (out[:n].copy() if n > 0 else None), time.perf_counter() - t, fs.value, ne.value


files = sorted(glob.glob(os.path.join(sys.argv[1], '*.npz')))
weight = float(sys.argv[2]) if len(sys.argv) > 2 else 1.2
rows = []
for f in files:
    inst = load(f)
    ref = inst['path'] if len(inst['path']) else None
    p, t, tf, ne = guided(inst, weight, inst['max_exp'])
    rows.append((os.path.basename(f), inst['blk'].shape, ref is not None, p is not None, inst['t'], t, tf, ne,
                 cost(ref, inst) if ref is not None else None, cost(p, inst) if p is not None else None))
both = [r for r in rows if r[2] and r[3]]
print('instances %d: ref found %d, guided found %d, both %d' % (len(rows), sum(r[2] for r in rows), sum(r[3] for r in rows), len(both)))
print('time: ref %.1fs  guided %.1fs (field %.1fs)' % (sum(r[4] for r in rows), sum(r[5] for r in rows), sum(r[6] for r in rows)))
if both:
    rel = np.array([r[9] / r[8] for r in both])
    print('path cost guided/ref: mean %.3f  min %.3f  max %.3f' % (rel.mean(), rel.min(), rel.max()))
print('guided expansions total %.1fM' % (sum(r[7] for r in rows) / 1e6))
