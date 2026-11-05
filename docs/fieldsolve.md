# fieldsolve: a PCB field solver (Rust, Metal)

`fieldsolve` is a 3-D finite-difference time-domain (FDTD) solver for PCB layouts. It reads the KiCad Gerber
output the router flow already produces (copper with X2 net/pad attributes, Excellon drills, the `.gbrjob`
stack-up), meshes a region of interest and computes port impedances, S-parameters, loop inductance and
field snapshots. It runs on the CPU (rayon) or on Apple GPUs (Metal); both give the same numbers.

Typical checks:

- **Decoupling loop inductance.** Put a port at an IC power pin, short (or load) the decoupling capacitor,
  and read the inductance from Z(f) at low frequency.
- **RF frontends.** Ports at the radio pin and at the u.FL/antenna feed give the S-parameters of the
  matching network and the line between them, including the board's parasitics.
- **Field pictures.** Field snapshots at chosen frequencies export to VTK (ParaView) and load into Blender.

## Choice of method and reference

| Solver | Method | Licence | Fit for us |
|---|---|---|---|
| [openEMS](https://github.com/thliebig/openEMS-Project) | EC-FDTD on a rectilinear non-uniform mesh | GPLv3 | **Reference.** It uses the same method, so it can be fed the exact same mesh. It is the established open PCB solver: [gerber2ems](https://github.com/antmicro/gerber2ems) (Antmicro, Apache-2.0) drives it from KiCad Gerbers and was validated against VNA measurements. |
| [FastHenry2](https://github.com/hpretl/FastHenry2) | PEEC, magneto-quasi-static | free | A good second opinion for loop inductance. It gives L and R only, with no resonances or radiation, and needs a filament model of the copper. |
| Palace (AWS), Elmer | FEM, frequency domain | Apache-2.0 / GPL | Accurate, but meshing PCB copper for FEM (gmsh) is a project of its own. |
| MEEP | FDTD | GPL | Aimed at photonics; it has no lumped ports or PCB input. |

FDTD covers both use cases with one engine. Low frequencies (inductance) and RF (S-parameters, antennas)
come out of the same broadband run. Its stencil is a perfect GPU workload, and a layered board maps
naturally onto a rectilinear mesh. fieldsolve is written from the textbook method (Yee, Taflove; CPML after
Roden and Gedney). openEMS is only run externally as a cross-check, so its GPL does not apply to our code.

## Model

- **Geometry.** Copper is a zero-thickness perfect conductor on the mesh plane of its layer, as in openEMS
  PCB models and gerber2ems. An in-plane E edge is metal when its midpoint lies in copper. Gerber
  primitives (flashes, macros, tracks as capsules, regions, LPC clears) are rasterised in file order.
  Plated holes are solid conductor columns: z edges at the mesh nodes inside the drill, or at the nearest
  node for holes thinner than a cell.
- **Stack-up.** Dielectric thicknesses, εr and tan δ come from the `.gbrjob`, with 4.5 / 0.02 when it gives
  none. The bottom copper is at z = 0 and copper thickness is ignored. tan δ becomes a conductivity at
  `f_loss`. Dielectric fills cells inside the board outline. Solder mask is optional (`solder_mask`).
- **Mesh.** Uniform `cell` in x/y inside the ROI and `cell_z` (at least `min_cells_layer` cells per
  dielectric layer) through the board, growing by at most 1.3× per cell outside, up to λ/20. Every
  copper interface and every port coordinate is a mesh line.
- **Boundaries.** CPML (κ = 1, cubic σ grading, α for low-frequency absorption), `pml` cells on each face.
  A cropped ROI leaves the board's planes cut at the PML, which acts like infinitely extended planes. That
  is right for local loop inductance; for whole-board plane resonances, use the full outline as the ROI.
- **Lumped elements.** Ports, resistors and parallel capacitors sit on a box of edges: k in series and m in
  parallel, R·m/k and C·k/m per edge. They are semi-implicit (folded into the update coefficients) and so
  unconditionally stable. A port is a resistive voltage source. V is the line integral of E from + to −,
  and I = (Vs − V)/R exactly as the update applies it. Port edges override copper they overlap. A `short`
  is a perfect conductor.
- **Excitation and stopping.** The excitation is a Gaussian, openEMS's `SetGaussExcite(0, fmax)`. A run
  stops when the field energy falls `end_db` (default −70 dB) below its peak.

## Usage

```sh
cargo build --release --bin fieldsolve
kicad-cli pcb export gerbers --check-zones -o gbr/ board.kicad_pcb
kicad-cli pcb export drill --format excellon --excellon-separate-th -o gbr/ board.kicad_pcb
target/release/fieldsolve run spec.json            # --backend auto|metal|cpu, --out DIR
target/release/fieldsolve svg gbr/ board.svg       # the imported layout, layers toggleable
target/release/fieldsolve selftest                 # analytic validation cases
```

A decoupling-loop spec for the AMOLED board's ESP32-S3 3V3 pin 46 and its capacitor C112. Coordinates are
in KiCad's frame (`"coords": "kicad"`); pads are `REF.PIN` from the Gerber X2 attributes:

```json
{
  "gerbers": "amoled_gbr", "out": "run_decap",
  "roi": [9.0, 30.0, 18.5, 37.0], "cell": 0.1, "air": 2.0, "air_z": 3.0,
  "fmax": 3e9, "fmin": 10e6, "nf": 300,
  "ports": [{"name": "U1_VDD", "at": "U1.46", "from": "F.Cu", "to": "In1.Cu", "r": 10}],
  "elements": [{"kind": "short", "pos": "C112.1", "neg": "C112.2", "layer": "F.Cu"}],
  "l_band": [50e6, 300e6],
  "fields": {"freqs": [100e6, 1e9]}
}
```

Ports and elements come in two forms:

- **Vertical.** `"at"` (pad or `[x, y]`) with `"from"`/`"to"` layers; the + terminal is on `from`.
  Optionally `"size": [w, h]` spreads the port over an area.
- **In-plane, between two pads.** `"pos"`/`"neg"` with `"layer"`. A copper bar as wide as the smaller pad
  joins the pads, and a one-cell gap at the midpoint holds the element: the body of a two-terminal part.

`kind` is `port` (the default in `ports`), `short`, `R` (`"r"`) or `C` (`"c"`, in parallel). `"excite":
"all"` runs every port and writes the full Touchstone file; by default only the first port is excited.

Outputs in `out/`:

- `result.json`: per excited port, Z(f), |S11| in dB and the R + jωL fit over `l_band`.
- `sparams.sNp`.
- `scene.json`: exactly what was simulated, for the openEMS cross-check.
- `fields/`: see below, including the vector fields for `fieldsolve lines`.

## Validation

`fieldsolve selftest` and `cargo test --release --bin fieldsolve` run closed-form cases. The numbers below
are at the default mesh; Metal and CPU agree to about 1e-6 relative.

| Case | Quantity | fieldsolve | Reference | Error |
|---|---|---|---|---|
| Empty PEC box 30×20×10 mm, graded mesh | TM110 / TM210 / TM120 | 9.010 / 12.496 / 15.790 GHz | analytic | +0.03 / +0.04 / −0.07 % |
| Microstrip w 3, h 1.6 mm, εr 4.3 (probes, 0.2–2 GHz) | Z0, εeff | 49.94 Ω, 3.282 | Hammerstad–Jensen 51.14 Ω, 3.258 | −2.4 %, +0.7 % |
| same, 2× mesh | Z0, εeff | 50.68 Ω, 3.277 | | −0.9 %, +0.6 % |
| Plane pair 40×30 mm, h 0.5 mm, εr 4.3 | f10 / f01 / f11 | 1.770 / 2.350 / 2.975 GHz | cavity model with open-end extension | −0.9 / −1.0 / +0.2 % |
| Patch 10×10 mm over ground, h 0.2 mm | C | 20.86 pF | wide-microstrip C with fringing | +4.1 % |

Lessons from building these, which also apply to board runs:

- **Copper edges.** A zero-thickness copper edge sits best 1/3 of a cell inside the outermost metal mesh
  line (openEMS's "thirds rule"). With the edge exactly on a mesh line, the strip acts about half a cell
  wider and a microstrip's εeff came out 6 % high at 4 cells per substrate. The patch case keeps its edges
  on mesh lines and shows that bias. Board meshes are uniform, so edges fall at random positions; expect a
  few percent from this and check convergence with a finer `cell` (below).
- **Port parasitics.** A port's own vertical current path adds inductance. The microstrip case therefore
  uses voltage probes and current loops far from both ends (`add_vprobe`, `add_iprobe`). Probes within
  ~10 mm of a discontinuity picked up evanescent box modes and biased εeff by 2–6 %.
- **Stopping.** End on field energy, not port voltage: −50 dB of energy is only −25 dB of field, which
  truncated the ring-down.

### Board case: AMOLED ESP32-S3 3V3 decoupling (U1.46 → C112 → In1 GND)

| cell | cells | steps | Metal time | L |
|---|---|---|---|---|
| 0.10 mm | 0.46 M | 37 k | 4.2 s (4.1 Gcell/s) | 0.987 nH |
| 0.05 mm | 1.97 M | 99 k | 72 s (2.7 Gcell/s) | 0.971 nH |

The CPU reference does 200 Mcell/s on the same problem.

openEMS on the identical scene (0.1 mm mesh, same ports, mesh lines, copper, vias and excitation) gives
**L = 1.035 nH against fieldsolve's 0.987 nH (−4.6 %)**. |Z| agrees within 3–4.6 % from 0.3 to 3 GHz, a
steady offset rather than noise. openEMS needed 20 min on the CPU (400 k steps at 150 Mcell/s) where
fieldsolve on Metal needed 4.2 s.

**Open issue: port placement.** The reciprocal set-up should give the same loop inductance: a port across
C112's pads (in-plane, bar with a one-cell gap) and a short from U1.46 to In1. It gives 0.92 nH at 0.1 mm
and 0.85 nH at 0.05 mm, while the original set-up gives 0.987 / 0.971 nH. The gap grows with refinement.
The suspects are the single-node vertical column (a zero-radius thin wire whose inductance depends on the
cell) and the in-plane gap port's current distribution. Until this is resolved, take absolute loop
inductance to about ±10–15 %. Comparisons between layouts on the same mesh and port set-up are tighter.

### Cross-check with openEMS

```sh
~/opt/openEMS/venv/bin/python tools/openems_scene.py RUN_DIR
```

This rebuilds `RUN_DIR/scene.json` in openEMS with fieldsolve's mesh lines (no smoothing):

- copper polygons as zero-thickness metal, in priority order;
- dielectric slabs extruded from the outline, with tan δ as conductivity at the same `f_loss`;
- barrels as metal lines on the same nodes;
- `AddLumpedPort` / `AddLumpedElement` on the same boxes;
- the same Gaussian, PML and end criterion.

It prints both impedances and the inductance fit over the same band. openEMS was built from
`openEMS-Project` with `update_openEMS.sh ~/opt/openEMS --python --disable-GUI`. The Homebrew formula is no
longer maintained.

## Fields and Blender

`fields/` holds, per snapshot frequency:

- `fields_<f>.vtr`: ParaView; |E|, |H| and the real-part vectors on the mesh nodes.
- raw `E_abs_<f>.f32` / `H_abs_<f>.f32` volumes.
- `J_<layer>_<f>.f32`: |n × (H⁺ − H⁻)|, the surface current density on each copper layer.
- `copper_<layer>.obj` and `meta.json`.

```sh
/Applications/Blender.app/Contents/MacOS/Blender --python tools/blender_fields.py -- RUN_DIR/fields \
    --freq 0.1000GHz [--field E|H] [--layers F.Cu,In1.Cu] [--no-copper] [--no-volume] [--top] \
    [--save scene.blend] [--render out.png]
```

The script builds the copper layers, one emissive plane per layer showing log |J| over 40 dB, and the field
magnitude as an OpenVDB volume (resampled to the finest cell, log over 50 dB). It uses only Blender's own
numpy and openvdb. Scene units are millimetres in Gerber orientation, with z up.

### Field lines

```sh
target/release/fieldsolve lines RUN_DIR/fields --freq 0.1000GHz --field H --seeds 300 [--phases 12] \
    [--floor-db -40] [--max-len 60] [--crowd 3]
```

This traces streamlines of the instantaneous field Re[F·e^{jφ}] from the snapshot's complex node vectors
(`E_vec_<f>.f32` / `H_vec_<f>.f32`):

- RK4 on trilinear interpolation over the non-uniform node grid, with a step of half the local cell;
- seeds drawn in proportion to |F|, outside the PML;
- an occupancy grid keeps the line density even;
- closed loops (H around a trace) are detected.

`--phases N` traces N phases over one period. The output is `fields/lines_<F>_<f>.json`.

```sh
Blender --python tools/blender_fields.py -- RUN_DIR/fields --freq 0.1000GHz --lines lines_H_0.1000GHz.json \
    --lines-mode volume|tubes|both [--line-radius MM] [--fog-cut 0.5] [--render out.png]
```

The lines become an edge mesh with a log-magnitude attribute `mag`. A Geometry Nodes modifier turns them
into an emissive fog volume (points → volume, density and glow from the grid, thicker where the field is
stronger, the weakest part dropped below `--fog-cut`) and/or emissive tubes coloured by `mag`. With several
phases there is one object per phase, keyed to show on successive frames: the field animated over a
period. The camera frames the lines, unless you pass `--frame-all`.

Renders use Cycles on the Metal GPU. For interactive viewing, Eevee's volume range is set to the scene size:
its default 0.1–100 units turns millimetre fog into blocks.

![H field lines and surface current, AMOLED U1.46 → C112 decoupling loop at 100 MHz](img/fieldsolve_decap_hlines.jpg)

## Code map

The code is in `src/bin/fieldsolve/`.

| File | Role |
|---|---|
| `gerber.rs` | RS-274X/X2 parser (macros, regions, arcs, polarity, step-repeat), Excellon (slots), `.gbrjob` stack-up, outline chaining, SVG view |
| `geom.rs` | shapes and inside tests |
| `mesh.rs` | graded 1-D line generation with fixed lines |
| `model.rs` | layout → per-edge coefficients, metal, barrels, lumped elements, CPML tables, probes |
| `fdtd.rs` | CPU reference stepping |
| `metal_fdtd.rs` + `shaders/fdtd.metal` | GPU stepping |
| `post.rs` | DFT, Z, power-wave S, R/L fit, Touchstone |
| `sim.rs` | spec → setup → runs → outputs; scene export |
| `export.rs` | VTK, raw volumes, complex vector fields, surface currents, copper OBJ |
| `lines.rs` | field-line tracing |
| `cases.rs` | analytic validation cases |
| `tools/openems_scene.py` | openEMS replay of a scene |
| `tools/blender_fields.py` | Blender import |

## Next steps

The plan for radio frontends (EXT balun, nRF54 filters, u.FL) and antennas (free space, enclosure, pipe) is in
[fieldsolve-rf-roadmap.md](fieldsolve-rf-roadmap.md). Engine items:

- **Port placement accuracy.** Resolve the reciprocal-port discrepancy, for example with a thin-wire model
  for single-node vertical ports or by giving ports a finite footprint, and close the 4.6 % offset from
  openEMS.
- **Capacitor models.** Series R-L-C elements, so a real capacitor (ESR, ESL, C) can replace the short.
  That gives the PDN impedance curve |Z(f)| directly instead of only the loop L.
- **Conductor loss.** A surface-impedance (skin-effect) model for copper, for RF insertion loss and the
  resistive part of loop impedance. Copper is lossless now.
- **Antennas.** A near-to-far-field transform (radiation pattern, efficiency) for chip and trace antennas.
- **Ports and meshing.** Ports from net names (all VDD pins of a part to its GND pins); the thirds rule
  applied in the rasteriser (shrink copper by 1/3 cell at edges) instead of relying on random edge
  positions.
- **Performance.** Fuse the CPML corrections into the main kernels and keep fields in half-precision
  storage where safe; most time on large runs is memory traffic.
- **Router input.** Read gridroute/topor `layout.json` directly, so a routed candidate can be checked
  without a KiCad round trip.
