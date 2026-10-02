"""GPU vs CPU cost-to-target field on corpus instances: agreement and time."""
import ctypes, glob, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridroute as gr
from bench_astar import load
L = gr._LIB; P = ctypes.POINTER
L.gr_field.argtypes = [ctypes.c_int] * 3 + [P(ctypes.c_uint8)] * 3 + [P(ctypes.c_float), P(ctypes.c_uint8), ctypes.c_float, ctypes.c_int, P(ctypes.c_float), P(ctypes.c_int64)]
def field(inst, backend):
    nl, h, w = inst['blk'].shape
    u8 = lambda a: np.ascontiguousarray(a).view(np.uint8).ctypes.data_as(P(ctypes.c_uint8))
    mc = np.ascontiguousarray(inst['mcost'], dtype=np.float32).ravel()
    out = np.empty(nl * h * w, dtype=np.float32); r = ctypes.c_int64(0)
    t = time.perf_counter()
    b = L.gr_field(nl, h, w, u8(inst['blk']), u8(inst['vok']), u8(inst['tgt']), mc.ctypes.data_as(P(ctypes.c_float)),
                   u8(inst['lay_ok'].astype(np.uint8)), inst['vcost'], backend, out.ctypes.data_as(P(ctypes.c_float)), ctypes.byref(r))
    return out, time.perf_counter() - t, b, r.value
def main():
  files = sorted(glob.glob(os.path.join(sys.argv[1], '*.npz')))
  step = int(sys.argv[2]) if len(sys.argv) > 2 else 1
  tc = tg = 0
  for f in files[::step]:
      inst = load(f)
      a, t1, _, _ = field(inst, 1)
      b, t2, be, rounds = field(inst, 2)
      fin = np.isfinite(a)
      ok = np.array_equal(fin, np.isfinite(b))
      err = ((b[fin] - a[fin]) / np.maximum(a[fin], 1)).max() if fin.any() else 0; err_lo = ((a[fin] - b[fin]) / np.maximum(a[fin], 1)).max() if fin.any() else 0
      tc += t1; tg += t2
      print(os.path.basename(f), inst['blk'].shape, 'cpu %.3fs gpu %.3fs rounds %d  inf-match %s  gpu-high %.2e gpu-low %.2e' % (t1, t2, rounds, ok, err, err_lo), flush=True)
  print('total cpu %.2fs gpu %.2fs' % (tc, tg))

if __name__ == '__main__':
    main()
