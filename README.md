# gridroute

Fast building blocks for grid-based PCB routers, and a board toolkit built on them.

- **Native kernels** (Rust, C ABI, multi-threaded CPU, Metal on Apple GPUs) with a numpy/ctypes binding: clearance masks, copper islands, multi-layer A*, a GPU cost-to-target field for guided search, content hashing.
- **`gridroute.board`**: a scriptable layout toolkit for "routing by program". Layout scripts place footprints explicitly and route nets in an order and on layers they choose. The toolkit provides the rasters, the router, plane drops, rip-up and reroute, speculative search and a persistent route cache. It saves the result as JSON for a board writer.
- **`gridroute.kicad`**: KiCad XML netlist loader, and a footprint-geometry dump that runs under KiCad's Python.

It was built for a scripted KiCad layout flow, a 4-signal-layer, 483-part, 170 × 111 mm board on a 0.05 mm grid. On it, the full placement-and-routing script went from **390 s to 5.7 s** cold and **1.0 s** for an unchanged re-run. The resulting board passes KiCad DRC.

## Build

```sh
cargo build --release      # target/release/libgridroute.{dylib,so}; the Python binding also builds it on first import
cargo test --release       # optimised, reference, CPU and Metal implementations agree
python3 examples/demo_board.py
```

The Python side needs numpy. Without a GPU the field falls back to an exact CPU Dijkstra. On other platforms the Metal code is compiled out and everything runs on the CPU; this has not been tested outside macOS yet.

## Board toolkit

```python
from gridroute import board
from gridroute.kicad import load_netlist, load_footprints

board.configure(layers=['F.Cu', 'In2.Cu', 'In3.Cu', 'B.Cu'], planes={'GND': 'In1.Cu', '+3V3': 'In4.Cu'},
                classes={'Default': (0.15, 0.15, 0.5, 0.3), 'Power': (0.5, 0.2, 0.6, 0.3)},
                net_class=lambda n: 'Power' if n.startswith('+5V') else 'Default')
bd = board.Board(170, 111, *load_netlist('project.xml'), load_footprints('footprints.json'))
bd.place('U1', 40.0, 30.0, 90)
bd.fanout_all(['U1'])                         # plane drops for the plane-net pads
bd.connect('SDA', margin=5.0, layers=['F.Cu', 'In2.Cu'])
bd.save('layout.json', [170, 111])
```

`examples/demo_board.py` is a complete, self-contained example with inline footprints.

Geometry comes from `python3 -m gridroute.kicad footprints project.xml footprints.json`, run with KiCad's Python. Board writing is left to the project: read `layout.json` in a pcbnew script.

### What it does for you

- **Rasters.** `occ` holds copper grown by the class clearance, `core` the true copper, `thru` the holes, `smd` the SMD pads.
- **Searches.** A route is a multi-layer A* in a window around its pads, with direction costs per layer, a turn penalty and vias.
- **Hybrid search (default).** A search still open after `budget` expansions switches to a GPU-computed cost-to-target field and continues as field-guided A*. On hard searches this cuts expansions by about 700× at the same path cost, and it finds routes that plain A* gives up on. `configure(search='exact')` keeps plain A* and sets `board.REPRODUCE`, so layout scripts can keep reference results byte for byte.
- **Rip-up.** Every rasterising operation is logged. `rip(nets)` clears the region those nets painted and replays the remaining operations there in order, so the rasters are exactly as if that copper had never been added. A blocked net can rip its neighbours, route first, and let them re-route.
- **Speculation.** Upcoming searches can run on worker threads (`speculate()`, `speculate_connect()`). A result is used only if no copper painted since touches its window, so it is always the sequential result.
- **Cache.** Routes, plane drops, clearance checks and cleanup are memoised on disk (`~/Library/Caches/gridroute/`), keyed by the board content around them through incremental tile hashes. A re-run repeats only what an edit touched. Searches are also memoised by their inputs and stay valid while the tiles they read are unchanged.
- **Relaxation.** Routes found early are shaped by copper that later moved, by layer direction preferences and by the weighted, field-guided search, so they take detours. `relax(since=k)` rips each net's copper logged since op index `k` and re-routes it against the finished board with an unweighted search and no direction preferences; the new route stays only if the net is complete and it is shorter (track length plus a via cost), otherwise the old copper goes back. Copper logged before `k` (hand-drawn breakouts, plane drops) is never touched. `pull_tight(since=k)` then string-pulls each track polyline: a run of vertices is replaced by a straight or one-bend 45-degree link when the link is clear and no other copper of the net attaches to the run.
- **Checking.** `accel='verify'` runs the numpy paths and the reference A* alongside the native kernels and asserts identical results. `check()` reports raster clearance clashes; KiCad DRC remains the final word.

Options default from the environment: `GRIDROUTE_SEARCH`, `GRIDROUTE_BUDGET`, `GRIDROUTE_FIELD`, `GRIDROUTE_ACCEL`, `GRIDROUTE_SPEC` and `GRIDROUTE_CACHE`. See `configure()`.

### Rules for robust layout scripts

These are lessons from the board above:

- **Place plane drops before a stage's signal routes.** Draw hand-made copper that is added unchecked (`add_track`, `add_via`) before the drops.
- **Use local rip-up on failure.** Rip the nets around the blocked net's pads, route it first, and promote any ripped net that then fails. Re-running a whole stage in a new order should be the last resort.
- **Print problems loudly.** Clashes, plane pads without a drop and incomplete nets should all reach the output.
- **Relax at the end.** Mark the op index before the router stage (`k = len(bd.ops)`), route, then `bd.relax(since=k)` and `bd.pull_tight(since=k)`. Pass `layers=` to keep a class (for example power) on its layers.
- **Stress the script by varying the search budget.** It is a cheap way to push the router into other configurations.

## Kernels (Python binding)

```python
import gridroute as gr
blk = gr.dilate(occ, layers=[0, 3], win=(i0, j0, i1, j1), span=gr.disc_span(reach / pitch), excl=(net_id,))
isl = gr.island(core, thru, net_id, win, seeds)
path = gr.astar(blk, vok, src_states, tgt, mcost, lay_ok, vcost, turn, hmul, max_exp)                   # exact
path = gr.astar(blk, vok, src_states, tgt, mcost, lay_ok, vcost, turn, hmul, max_exp, budget=50_000)   # hybrid
```

Also `island_crop`, `label_cells`, `mask_scan`, `paint_capsule`, `hash_windows`, `hash_arrays`, `tile_hashes`, `hash_tiles` and `SearchCache`. The C API is in [`include/gridroute.h`](include/gridroute.h).

## Measurements (M3 Max)

Full layout script on the board above:

| Configuration | Time | Result |
|---|---|---|
| numpy masks, Python flood fill, C A* (before) | 390 s | reference board |
| `search='exact'`, cold | 92 s | byte-identical to the reference |
| hybrid search (default), cold | 5.7 s | different board, all nets routed, KiCad DRC clean |
| hybrid, unchanged re-run | 1.0 s | byte-identical to the cold run |
| hybrid, re-run after moving one capacitor 0.5 mm | 3.0 s | byte-identical to an uncached run of the edit |

Per kernel:

- **Clearance dilation:** 4× (5 × 5 mm window) to 136× (whole board) faster than numpy on the CPU. The Metal version gives the same masks but is no faster: the work is memory-bound and unified memory leaves no transfer to save, so `auto` keeps it on the CPU.
- **Exact A\*:** memory-latency bound. It is 1.2× faster than the C original, with the same paths. A concurrent reachability check ends hopeless searches early: 5.9 s to 0.5 s on a board-sized one.
- **Cost-to-target field:** about 9.5× faster on the GPU than CPU Dijkstra, with identical values. At 2× coarser resolution (`set_field_factor(2)`, the hybrid default) it costs about 15 ms for a hard search.
- **Painting and scans** (`paint_capsule`, `mask_scan`): the same arithmetic as numpy, so the results are bit-identical.

### Why nets are not routed in parallel

Measured on the same board: the dozen board-spanning nets overlap one another and must stay sequential, and the other ~290 nets take about 3 ms each, which is about what handing one to a thread costs. An ideal 8-worker schedule would cut the global stage only from 3.55 s to 3.07 s, before thread overhead and the GIL.

What mattered was elsewhere. For example, `connect()` used to ask, pad by pad, whether each pad was already joined to the tree. It now labels the net's copper components at all pads in one native pass (`label_cells`).

## Tools

- `bench/bench_kernels.py`: dilation numpy vs CPU vs Metal, by window size.
- `bench/bench_astar.py CORPUS`: replays recorded A* instances, optimised vs reference, path for path.
- `bench/field_check.py CORPUS`: GPU vs CPU field agreement and time.
- `bench/bench_field_hybrid.py CORPUS`: hybrid search field time.
- `bench/guided_experiment.py CORPUS`: guided vs reference search (success, path cost, expansions).

A corpus is a directory of `.npz` instances: `blk`, `vok` and `tgt` bit-packed, plus `shape`, `src`, `mcost`, `lay_ok`, `par` = (vcost, turn, hmul), `max_exp`, and the recorded `path`. Record one by wrapping `gridroute.astar` in a router run and saving its arguments.

## License

MIT; see [LICENSE](LICENSE).
