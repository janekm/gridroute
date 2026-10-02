"""Time the hybrid search's field phase on a corpus (budget 1 forces it): total field seconds, total time, paths found.

    python3 bench/bench_field_hybrid.py CORPUS [factor] [step]
"""
import glob, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gridroute as gr
from bench_astar import load

files = sorted(glob.glob(os.path.join(sys.argv[1], '*.npz')))[::int(sys.argv[3]) if len(sys.argv) > 3 else 1]
gr.set_field_factor(int(sys.argv[2]) if len(sys.argv) > 2 else 2)
gr.astar(*[np.zeros(1)] * 0) if False else None
t_all, found = 0.0, 0
f0 = gr.stats['field_s']
for f in files:
    i = load(f)
    t = time.perf_counter()
    p = gr.astar(i['blk'], i['vok'], i['src'], i['tgt'], i['mcost'], i['lay_ok'], i['vcost'], i['turn'], i['hmul'],
                 i['max_exp'], budget=1, weight=1.2)
    t_all += time.perf_counter() - t
    found += p is not None
print('%d instances: field %.2fs, total %.2fs, found %d' % (len(files), gr.stats['field_s'] - f0, t_all, found))
