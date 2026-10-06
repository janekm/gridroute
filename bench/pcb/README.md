# Native PCB routing benchmark

This harness routes immutable placed boards with gridroute and checks the exported
copper in KiCad. The original JLink-SWD, FT231X and 1Bitsy cases are regression
fixtures. Further PCBench cases and the supplied six-layer AMOLED board exercise
additional geometry, congestion and constraint handling.

## Acceptance

A candidate passes only when native KiCad reports zero unconnected items, no new
error-severity findings relative to the stripped input, and unchanged footprint,
static-board and project-rule content. Warnings are counted separately and are
not silently treated as errors or hidden. Existing source defects remain visible.
This gate does not qualify differential-pair skew, impedance, power integrity or
manufacturing readiness.

The `input/` directory beside each result contains the exact model, unrouted
board, project, custom rules and baseline DRC used by that run, plus SHA-256 hashes.
`routing.json` records options, device, timing, source hashes and recovery events.
`acceptance.json`, `geometry-check.json` and `routed-drc.json` are separate evidence.
Do not infer connectivity or DRC acceptance from a successful process exit.

PCBench profiles use the pinned dataset's class widths as minimum routing widths.
A separate USBee32 profile permits the native board's 0.2 mm pad escapes while
retaining a 0.4 mm preferred routing width. Placement and footprint-correction
experiments have distinct case IDs and cannot count as unchanged-input passes.

## Run

Build the native library with `cargo build --release` from the repository root.
Install NumPy in the routing Python environment. KiCad's `pcbnew` bindings run in
KiCad's own Python environment; the routing process does not require them.

```sh
PYTHONHASHSEED=0 GRIDROUTE_CACHE=off GRIDROUTE_NO_BUILD=1 \
  python3 bench/pcb/suite.py --out /absolute/path/to/results \
  --cases JLink-SWD_JLink-SWD FT231X_breakout_FTDI_FT231XS-U_Breakout 1Bitsy_1bitsy
```

The suite's macOS KiCad defaults can be overridden with `--kicad-python` and
`--kicad-cli`. Check `routing.json` for the actual Metal availability and device.
CPU fallback is valid, but must not be compared with Metal timings as the same
configuration. Runs are serial because Board configuration is process-global.
Use a fresh output directory for every run; existing evidence is never overwritten.

The preferred validation option is `--options-json '{"continuous_check":true}'`.
It runs continuous copper connectivity and standard clearance checks in the routing
process, then tries bounded repairs before a single final KiCad export and DRC.
The checker includes orphan copper islands, rotated/rounded/polygon pads, slotted
holes, copper clearance, hole clearance and hole-to-hole clearance. It uses a
sweep-line broad phase and exact segment/polygon distances, independent of the
routing raster. Static pad/pad defects remain in the source baseline.

`gridroute.continuous.check(board)` is usable without KiCad. Its companion
`repair_connectivity(source, factory, ...)` reroutes at most four affected nets,
tries 0.025 and 0.0125 mm, and shares a 30-second cooperative repair budget. It
retains only candidates with no standard clearance findings and improved
connectivity, or resolved clearance findings without worse connectivity. It
preserves the source board and disables destructive cleanup during repair.
Dense incomplete boards receive diagnostics without expensive additional retries.

This is an explicitly limited DRC subset. It does not evaluate arbitrary custom
rule expressions, solder-mask webs, copper-to-edge/Margin distance, keepout areas,
zone fills or differential-pair constraints. Routing still uses raster obstacles
for these supported input features; final native validation is mandatory. The
checker is presently Python: measured low-millisecond checks do not justify
moving it into Rust before profiling larger boards and incremental updates.

`--native-feedback` retains an optional diagnostic/reference path. It uses native
gaps to generate up to two repair candidates and repeats CAD validation. Prefer
the continuous option above for normal iteration. Per-case wall time includes
all attempts and CAD checks; routing time includes continuous checks and accepted
or rejected in-process repair attempts. Warnings remain separate.

For one candidate:

```sh
python3 bench/pcb/run.py 1Bitsy_1bitsy --out /absolute/path/to/result \
  --options-json '{"fallback_pitches":[0.025],"deadline_seconds":60,"continuous_check":true}'
/path/to/kicad-python bench/pcb/kicad_io.py export 1Bitsy_1bitsy \
  --out /absolute/path/to/result
kicad-cli pcb drc --format json --severity-all --all-track-errors \
  -o /absolute/path/to/result/routed-drc.json /absolute/path/to/result/routed.kicad_pcb
python3 bench/pcb/accept.py 1Bitsy_1bitsy --out /absolute/path/to/result
```

Routing time includes board construction, search and cleanup. It excludes process
startup and the native export/DRC gate. `suite.py` also reports complete wall time.
The in-library deadline is cooperative; the suite applies an outer process timeout.

## Fixtures and provenance

`manifest.json` records source revision, source hashes, physical profile and any
explicit input limitation. PCBench is pinned to
`dec3be75cbdef74787625f9043c7391cd473bb64`. Dataset models, legacy native CAD files
and source metadata are evidence, not executable inputs. Preserve their original
metadata and source attribution when redistributing fixtures.

To prepare a PCBench fixture from a local pinned dataset directory containing
`<case>/processed.kicad_pcb` and `<case>/final.json`:

```sh
/path/to/kicad-python bench/pcb/kicad_io.py prepare CASE --source-dir /path/to/PCBs
kicad-cli pcb drc --format json --severity-all --all-track-errors \
  -o bench/pcb/cases/CASE/baseline-drc.json bench/pcb/cases/CASE/unrouted.kicad_pcb
```

The adapter supports rotated and polygon pads, slotted drills, polygon outlines
with cutouts, pad-local clearance, mask openings, copper graphics, Margin graphics
and copper-layer aliases. Invalid open/disjoint outlines fail explicitly. Legacy
nil object UUIDs receive deterministic IDs because KiCad otherwise resolves DRC
items to arbitrary footprints. This changes identity metadata, not geometry.

The supplied AMOLED model preserves its exact project and custom rules, rule
areas, layer minimums and allowed pad neck-downs. RF's local courtyard clearance
relaxation is conservatively unused. Cold routing removes all copper zones; a
separate plane profile retains the four original zone definitions and refills them
in KiCad. Neither is equivalent to preserving the supplied routed copper.

## Policy and performance

`gridroute.router.RoutingController` starts with the normal grid and retries failed
nets at the configured finer pitches. Physical rules remain constant while the
factory scales numerical raster allowances. Already-complete routes avoid a fine
grid rebuild unless `verify_at_finest=True`. `retain_fine_grid` is optional because
keeping the fine grid can alter subsequent routing order and geometry.

Recovery identifies the actual obstructing tracks and vias, reroutes the target,
then reconnects affected nets transactionally. Default acceptance requires fewer
disconnected pad components. `plateau_budget` permits a bounded number of distinct
equal-score moves. `candidate_paths` can search a few alternative corridors while
preserving disruptive blockers. Both are experimental escalation controls.

Static-board snapshots reduce repeated geometry work. The default raster cache is
bounded to 128 MiB; `cache_static=False` disables it. The initial coarse pass does
not allocate an unused snapshot; templates are created only for retries. A board factory must return
independent, deepcopy-compatible static geometry. Custom methods must operate on
`self` and must not capture another mutable Board in a closure.

Hole masks use local bounding boxes, blocker tests reject distant object/segment
boxes before exact distances, pad flood seeds are compact, and native component
labels avoid repeated full-board floods. Native Rust/Metal kernels are unchanged.
Numerical safety still requires the native CAD gate: a grid result alone is not a
continuous-geometry clearance proof.

```sh
PYTHONPATH=python GRIDROUTE_NO_BUILD=1 python3 -m unittest discover -s tests -p 'test_*.py'
cargo test --release
```
