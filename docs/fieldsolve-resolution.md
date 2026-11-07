# fieldsolve: resolution where the physics needs it

fieldsolve uses a rectilinear (tensor-product) Yee mesh with uniform cells over the copper window, graded
outside it. Every copper edge, gap and port is resolved by making the whole window fine. This page asks how
to put resolution only where it is needed. It reports a measured prototype of dynamic resolution (step 2),
compares the candidate methods (step 3) and recommends one.

Short answer: **keep the rectilinear FDTD and spend the effort on sub-cell accuracy and a-priori
refinement, not on dynamic or unstructured meshes.**

- The sub-cell copper-edge corrections in `edge.rs` bought, at fixed cell size, what 2-4× finer cells buy
  without them (see [fieldsolve.md](fieldsolve.md#sub-cell-copper-edges)). In 3-D that is 8-64× fewer cells
  and 2-4× fewer steps.
- Field-driven dynamic refinement does not pay on PCB problems. Once a pulse has filled a board, the field
  is everywhere: refinement follows it, and the field transfers at each re-mesh add error.
- An octree with geometry-driven refinement would need 1.9-3.2× fewer cells than the present mesh on the
  nRF52840 DK board, plus about 1.3-1.6× from local time stepping. That is a large rewrite (GPU data structures,
  interface kernels, stability) for a moderate gain.

## Step 2: a dynamic-resolution prototype

`amr2d.rs` (CPU, not used by `run`; `fieldsolve amr-check [h_mm] [steps] [absorb_mm]`) solves the 2-D TMz
equations of a plane pair seen from the top (Ez between the planes, H in the plane) on a 2:1-balanced
quadtree. This is the in-plane part of the 3-D problem, the one a dynamic scheme would coarsen.

- **Scheme.** Ez is a node voltage at each leaf's centre with C = ε A, the tangential H on each face is a
  branch current with L = μ d / w (d the distance of the two centres along the face normal, w the face).
  On a uniform grid this is exactly Yee. On a quadtree it is the two-point-flux finite-volume scheme: first
  order at hanging faces, but a passive LC network, so the discrete energy ½ΣCV² + ½ΣLI² is conserved
  exactly between re-meshes and the scheme is stable for dt ≤ 2/√λ_max of the network.
- **Static refinement.** Ports, probes, material boundaries and thin PEC walls (with their tips) get the
  finest cells within two cells. Elsewhere a leaf may grow to λ_min/`ppw` at fmax (and 2:1 balance).
- **Dynamic refinement.** Every 10 steps each leaf's local truncation error is estimated as |mean of the
  neighbours − own value| (h²/4 times the Laplacian, ≈ (hk)²/4 of the amplitude). It is compared with the
  running peak. Busy leaves split; sibling groups that are all quiet merge, down to a ppw-5 floor. Fields
  move through a fine-grid map: Ez as charge-conserving averages, H linear inside a leaf so that a split
  face carries the same current. Thin-wall surface currents are kept per side.
- **Time step.** One global dt, set by the finest cell that may ever appear: refinement during a run must
  not break the CFL limit. Local time stepping (each level at its own step; for example Diaz and Grote,
  SIAM J. Sci. Comput. 31(3), 2009) is not implemented. The run counts what 2:1 local time stepping would
  have needed on the same meshes.

Test: a 40 × 30 mm plane pair, εr 4.3, tan δ 0.02, open (magnetic) edges, a 15 mm slot (thin PEC wall) from
one edge, a 50 Ω port near a corner. Fine cell 0.25 mm. Reference: the same scheme on the uniform fine grid.
|Z| errors are rms over 0.3-5 GHz.

**Closed board** (no absorber, 60 000 steps):

| mesh | leaves | work (cell updates) vs uniform | with 2:1 local time steps | \|Z\| rms error | modes near f10 / f01 / f11 / next (GHz) |
|---|---|---|---|---|---|
| uniform 0.25 mm | 20 480 | 100 % | 100 % | (reference) | 1.53 / 2.28 / 3.25 / 4.19 |
| static, ppw 20 | 3 068 | 15.0 % | 10.0 % | 0.29 dB | 1.53 / 2.31 / 3.25 / 4.19 |
| static, ppw 10 | 2 429 | 11.9 % | 9.1 % | 0.44 dB | 1.53 / 2.31 / 3.25 / 4.19 |
| static, ppw 5 | 2 321 | 11.3 % | 9.0 % | 0.74 dB | 1.56 / 2.31 / 3.25 / 4.19 |
| dynamic, refine at 0.025 | 2 321-17 339 | 26.6 % | 23.3 % | 1.13 dB | 1.53 / 2.31 / 2.90 / 4.19 |
| dynamic, refine at 0.01 | 2 321-20 438 | 52.8 % | 49.6 % | 0.35 dB | 1.53 / 2.31 / 2.87 / 4.19 |

The same board at a 0.125 mm fine cell (81 920 leaves uniform, 120 000 steps, 238 s):

| mesh | leaves | work vs uniform | with 2:1 local time steps | \|Z\| rms error | modes (GHz) |
|---|---|---|---|---|---|
| uniform 0.125 mm | 81 920 | 100 % | 100 % | (reference) | 1.53 / 2.28 / 3.22 / 4.19 |
| static, ppw 20 | 5 312 | 6.5 % | 4.5 % | 0.29 dB | 1.53 / 2.28 / 3.22 / 4.19 |
| static, ppw 5 | 4 565 | 5.6 % | 4.3 % | 0.70 dB | 1.56 / 2.31 / 3.22 / 4.19 |
| dynamic, refine at 0.025 | 4 565-60 434 | 19.0 % | 16.5 % | 1.15 dB | 1.56 / 2.31 / 2.90 / 4.19 |
| dynamic, refine at 0.01 | 4 565-81 572 | 44.2 % | 41.1 % | 0.37 dB | 1.53 / 2.31 / 2.87 / 4.19 |

Halving the fine cell halves the static mesh's relative work again (the fine cells sit along lines: the
slot, the port, the edges), at the same accuracy. The dynamic mesh's errors stay.

**Open board** (10 mm air and an absorbing frame around it, 25 000 steps):

| mesh | leaves | work vs uniform | with 2:1 local time steps | \|Z\| rms error |
|---|---|---|---|---|
| uniform 0.25 mm | 65 536 | 100 % | 100 % | (reference) |
| static, ppw 20 | 8 554 | 13.1 % | 8.2 % | 0.62 dB |
| static, ppw 5 | 5 785 | 8.8 % | 7.0 % | 1.16 dB |
| dynamic, refine at 0.025 | 5 785-6 076 | 8.8 % | 7.0 % | 0.66 dB |
| dynamic, refine at 0.01 | 5 785-27 997 | 16.8 % | 14.1 % | 0.21 dB |

What this shows:

- **Static refinement does most of the work.** Fine cells only at the port, the slot and its tip and the
  board edge, wavelength-sized cells elsewhere: 85-90 % less work for 0.3-0.6 dB in |Z|.
- **Dynamic refinement does not pay on a resonant board.** The pulse fills the board within a few
  nanoseconds and the board rings for tens. From then on the busy region is the whole board: the refined
  mesh is nearly uniform (17 000-20 000 of 20 480 leaves), the work is 2-4× that of the static mesh, and
  the accuracy is no better. The f11 resonance came out 11-12 % low: each merge and split adds an error, and
  in a resonator these errors accumulate over the ring-down.
- **It helps a little where energy leaves.** On the open board the field passes through and is absorbed.
  There the dynamic mesh (refine at 0.01) got 3× closer to the reference than static ppw 20 for 1.3× the
  work. PCB S-parameter runs are the resonant case (they run until the ports have rung down). Open
  antenna runs are partly the other case, but the antenna and its feed are what rings.
- **Energy at re-mesh events.** A merge loses the energy of the field structure below the new cell: 0.1-4.4 %
  of the injected energy over a run, as an extra, unphysical loss. Splits add none (fields are
  interpolated). The scheme stays stable throughout: between re-meshes it is a passive network.
- **Coarse/fine interfaces** are first order (hanging faces, where the centre-to-centre line is not normal
  to the face). On the static meshes the modes shift by up to 1.3 % (closed) and 2.7 % (open), both at f01, the mode the
  slot perturbs most. That is the
  same kind of first-order error as staircasing.
- **Time step.** With one global step the coarse cells do far more updates than they need. 2:1 local time
  stepping would cut the work by a further 1.3-1.6× on these meshes, not more: most of the leaves of a
  well-graded mesh are already near the finest level.
- **Ports and the absorber.** The port leaf is kept at the finest level, so its own impedance does not move
  with the mesh. The absorbing frame ran on cells of up to 32 fine cells. Its loss grading is evaluated at
  each leaf's centre, so a coarse leaf sees a step instead of a ramp (not measured separately here). A
  real PML (stretched coordinates) would need its profile per level.

Measured on the M3 Max CPU, single-threaded. The uniform 0.25 mm run takes 16 s for 60 000 steps. In 3-D
the same structure holds, with the interface problems (late-time instability of non-symmetric
interpolation, reflections, PML on mixed levels) are known to be harder.

## Step 3: candidate methods

Each candidate is weighed against what limits fieldsolve's PCB runs today. In-plane, copper features
(0.1-0.3 mm gaps and traces, pads, antipads) set the cell, not the wavelength: at 8 GHz in FR4, λ/20 is
0.9 mm, while the RF path needs 0.05 mm. Through the board, thin layers (75 µm prepreg) set it. The
rectilinear mesh makes the whole window as fine as its finest feature, and the time step follows the
smallest cell.

| candidate | what it saves | accuracy | implementation in fieldsolve | risk | fit |
|---|---|---|---|---|---|
| **A. Static a-priori refinement on the rectilinear mesh**, driven by geometry (RF nets, copper edges, ports) and a coarse pre-run's error estimate, plus the sub-cell edge corrections | The corrections: a 2-4× coarser cell for the same edge accuracy (8-64× cells in 3-D). Refinement regions: fine lines only where copper needs them, instead of the whole window | 1st order, but the edge corrections remove the first-order edge error of zero-thickness copper | small: `mesh.rs` already takes refinement regions; the corrections are in | low: same kernels, same stability | excellent |
| **B. Octree / true subgridding with local time stepping** (stable FIT-type interfaces: Thoma and Weiland 1996; Chilton and Lee, IEEE TAP 55(9), 2007; Bérenger's Huygens subgridding, IEEE TAP 54(12), 2006) | 1.9-3.2× fewer cells on the DK board (estimate below), ×1.3-1.6 from local time stepping | 1st order at interfaces; reflections −30 to −40 dB are typical of published schemes | large: new data layout, interface and LTS kernels on Metal, CPML and lumped ports per level, conformal edges at interfaces | high: late-time instability is the known failure of non-symmetric schemes; symmetric ones cost accuracy | poor: everything (material table, ports, far-field box, exports) assumes one rectilinear grid |
| **C. Unstructured meshes: tetrahedral DGTD or FEM-TD** | cells only where needed; conformal copper | high order possible, conformal geometry | very large: a PCB tet mesher (zero-thickness copper, thousands of vias), a new solver, local time stepping (the smallest tet sets dt) | very high | none: a separate solver |
| **D. Higher-order FDTD stencils** (FDTD(2,4) and kin) | fewer points per wavelength (6-10 instead of 20) | better dispersion in bulk; worse at PEC and material boundaries without special treatment | moderate | moderate: wider stencils at every boundary | **does not address the limit**: board meshes are geometry-bound, not dispersion-bound (0.05 mm vs λ/20 = 0.9 mm) |
| **E. Frequency-domain FEM** (Palace, Elmer) or **planar MoM / surface integral** for the copper layers (layered-medium Green's functions: Michalski and Mosig, IEEE TAP 45(3), 1997), with a coarse volume model | unknowns on the copper surface only (MoM); adaptive mesh refinement (FEM) | MoM: very good for planar RF circuits (edge-singular basis functions); FEM: high order, adaptive | large (MoM: Green's functions, singular integrals, finite ground planes and cut-outs break the layered-medium assumption); FEM: meshing again | high | an independent cross-check, not a replacement: the time-domain outputs (TDR, transient fields, broadband in one run) go away |

### Evidence from this codebase

- **Edge corrections vs refinement.** On the microstrip (h = 1.6 mm, edges at seven positions in a cell),
  with both corrections, 4 cells across the substrate give εeff within 0.73 % and C within 2.3 % rms. Without
  them, 4 cells are 18 % out in εeff, and even 16 cells are 4.9 % out. See
  [fieldsolve.md](fieldsolve.md#sub-cell-copper-edges) for the tables.
- **On a real board.** On the nRF52840 DK the corrected 0.05 mm mesh gives the S-parameters the plain mesh
  approaches as the cell shrinks: extrapolated from 0.05 and 0.025 mm, the plain model lands on the corrected
  one (the fitted extra shunt capacitance per node: −0.33, −0.28, extrapolated −0.24 pF; corrected −0.24 and
  −0.25 pF). The corrected 0.05 and 0.025 mm runs agree to |ΔS| 0.016 below 3 GHz.
- **Where the corrections stop helping.** A 0.3 mm trace on 75 µm prepreg at 0.05 mm cells (1.5 cells
  across the substrate, as on the nRF52840 DK) keeps ±8 % in C after the corrections: the strip is only 6
  cells wide and the field under it 1.5 cells high. There the answer is finer cells around that trace,
  which is exactly what refinement by net (below) gives without refining the whole window. A 45° trace keeps
  a staircased current path at 0.2 mm (εeff +5 %), which again asks for local refinement, not a different
  solver.
- **Where an octree would save.** `fieldsolve octree-count spec.json` estimates the cells of a
  geometry-driven octree on a board's mesh. In-plane, cell size grows as h₀ + (r − 1)·d with the distance d
  from the nearest copper edge. In z it grows with the distance from the nearest copper plane, never finer
  than the existing mesh there.

  | nRF52840 DK, 0.05 mm (18.78 M cells) | grading r 1.3 | 1.5 | 2.0 |
  |---|---|---|---|
  | octree estimate | 9.85 M | 7.53 M | 5.92 M |
  | fewer cells | 1.9× | 2.5× | 3.2× |

  The saving is modest because a PCB has copper edges nearly everywhere: planes with antipads, pours,
  traces. The rectilinear mesh is already graded outside the region of interest.
- **Dynamic refinement**: no gain on resonant PCB runs (step 2).

### Recommendation: A

Static a-priori refinement on the existing rectilinear mesh, together with the sub-cell edge corrections.
In order:

1. Use `"edges"` on RF runs, and coarsen the copper window where the corrections allow it (validate per
   board with a 2-mesh check).
2. Refinement regions by net. Fine lines only around the RF path's copper (and what couples to it), with
   graded cells over planes and pours. `mesh.rs` already takes several `Refine` regions. What is missing is
   choosing them from the layout (nets, pads, gaps under a size).
3. An error-estimate pre-run. Run the board coarse and accumulate an indicator (local field energy ×
   (cell size × curvature)², as in the 2-D prototype) over the run. Then add mesh lines where it is
   large. Being static, this keeps the global time step but avoids refining the whole window.
4. Revisit B only if a board's runs are still too slow after 1-3. The octree's 2-3× would then come on top.

Not recommended: dynamic (field-driven) refinement (step 2), higher-order stencils (they address
dispersion, which is not the limit), and unstructured or frequency-domain solvers as a replacement. A planar
MoM would be a valuable independent cross-check of matching networks (like openEMS today), not a
successor.

## Implemented and measured: static refinement and an octree (Nov 2026)

Branch `claude/fieldsolve-refine-octree` (on top of the edge corrections and the fast-mode branch's
`"refine"` boxes) adds both and compares them.

### Static refinement on the rectilinear mesh (`meshpol.rs`, spec `"mesh"`)

```json
"mesh": {"no_slivers": true, "z_grade": 3.0, "copper_bands": {"layers": ["F.Cu"], "cell": 0.05, "margin": 0.1}}
```

- `no_slivers` fills the interval between two fixed lines with ceil(n − 0.35) cells instead of ceil(n).
  Fixed lines (port and part coordinates) a little over a cell apart used to leave two half cells, and the
  smallest cell sets the time step for the whole mesh: on the DK the finest x cell was 0.025 mm and y
  0.035 mm in a 0.05 mm mesh.
- `z_grade` keeps `cell_z` within two cells of every copper plane and dielectric interface and lets cells
  grow (1.3 per cell, up to 3 · `cell_z`) through thick cores.
- `copper_bands` makes the spec's `cell` the base of the copper window and adds bands of a finer cell
  around the listed layers' copper edges (found on a supersampled raster: x bands where copper changes
  along x, y bands along y).

`fieldsolve octree-count spec.json` now also prints the cells per axis by size. The DK at 0.05 mm: more than
half of the cells along each axis are coarse (≥ 8× the finest), so only about 5 % of the 18.8 M cells are
fine in all three directions. That is the tensor-product waste an octree recovers.

nRF52840 DK, configuration 1, all five excitations, with `"edges"`. "Work" is cells / dt, the cost per
simulated nanosecond, relative to the 0.05 mm mesh. The reference is the 0.025 mm run with corrections:

| mesh | cells | dt | work | max \|ΔS\| 0.1-3 / 3-6 / 6-8 GHz | best ΔC (edges / rms) |
|---|---|---|---|---|---|
| 0.05 mm (as before) | 18.8 M | 56.7 fs | 1.00 | 0.016 / 0.163 / 0.094 | −0.24 / −0.22 pF |
| + no slivers | 18.3 M | 71.4 fs | 0.77 | 0.019 / 0.045 / 0.113 | −0.23 / −0.22 pF |
| + z grading | 16.4 M | 71.0 fs | 0.69 | 0.017 / 0.041 / 0.104 | −0.23 / −0.22 pF |
| + 0.1 mm base, 0.05 mm copper bands | 16.0 M | 78.5 fs | 0.61 | 0.017 / 0.085 / 0.071 | −0.24 / −0.24 pF |
| 0.1 mm base, no bands | 10.6 M | 79.6 fs | 0.40 | 0.025 / 0.210 / 0.300 | −0.24 / −0.24 pF |

- **1.6× less work at the same in-band accuracy**, almost all of it from removing slivers (larger dt) and
  grading z. The copper bands add little: the RF path has copper edges at nearly every x and y, so the
  bands cover most of the window.
- A 0.1 mm window everywhere is 2.5× cheaper and still fine in band (with the edge corrections), but loses
  accuracy above 3 GHz.
- **The number of steps is not a reliable cost measure on this board.** Runs stop when every port's signal
  has fallen 60 dB. A weakly coupled board resonance (late-time Q 45-65) sometimes stays above that line for
  a port and sometimes not, so the same board took 63 k steps on one mesh and 250 k on another. The table
  therefore compares cost per simulated time. (The 3-6 GHz column is affected the same way: where a run, or
  the reference, stopped early, that resonance is truncated.)

### An octree (`octree.rs`, CPU prototype, `fieldsolve octree-check`)

A 3-D FIT scheme on a 2:1-balanced octree (faces and edges), with box cells (the finest cell may differ per
axis). Edge voltages and face fluxes on the leaf edges and faces split at hanging nodes; hanging edges are
constrained to the coarser leaf's edges (halves 1/2, edges inside a coarse face 1/4 + 1/4), Ampère acts
through Pᵀ, and the mass is the lumped Pᵀ M_ε P. The scheme is symmetric in the energy inner product, so
it conserves energy and is stable at the Yee time step of the finest cell (power iteration: dt²λ = 3.39
at 0.4 mm, limit 3.61; `refined_octree_is_stable_and_conserves_energy`). Strip current is exact in this
form (the curl term of the strip's metal edges).

Microstrip (w ≈ 3 mm, h = 1.6 mm, PEC box 14 mm), plain staircase (the edge corrections are not ported to
the octree), leaves sized by distance from the strip's edges with grading r:

| finest cell | mesh | cells / leaves | Z0 vs uniform | εeff vs uniform |
|---|---|---|---|---|
| 0.4 mm, 4 across h | uniform octree (= Yee) | 793 600 | | |
| | octree, r 1.3 | 44 000-49 000 (16-18×) | −0.16 to +0.09 % | −0.15 to −0.14 % |
| 0.2 mm, 8 across h | uniform octree (= Yee) | 2 571 000 | | |
| | octree, r 1.3 | 79 000 (32×) | −0.30 % | −0.03 % |
| | octree, r 1.5 | 49 000 (52×) | +0.68 % | −0.08 % |
| | octree, r 2.0 | 36 000 (71×) | +1.16 % | −0.14 % |
| | graded rectilinear (fieldsolve) | 250 000 | | |

- **The octree's interfaces cost little accuracy** (under 0.3 % at r 1.3, about 1 % at r 2), and against
  the graded rectilinear mesh it needs 3.2× (r 1.3) to 6.8× (r 2) fewer cells. With local time stepping
  it would do 77 % of the updates at r 1.3 (most leaves are near the finest level around the strip).
- **But it is not faster per unknown.** On the CPU, 79 000 leaves (472 000 unknowns) took the same 22 s as
  the rectilinear mesh's 250 000 cells (1.5 M unknowns): gathering through edge and face lists costs about
  3× more per unknown than the structured stencil. On the GPU the structured kernels run at 2 000 Mcell/s
  because neighbours are implicit; an octree kernel is gather-bound, so the 3-7× cell saving would turn into
  perhaps 1-2× in time, after a rewrite of the GPU backend, ports, CPML, far field and exports.
- Against Hammerstad-Jensen the octree runs (uniform and refined alike) are within −3.8 to +3.1 % in Z0,
  against −5.6 to −0.3 % for the rectilinear mesh at the same insets. The setups differ: a closed PEC box
  (doubling it moved Z0 by 0.7 %) against absorbing walls. The comparison above is therefore octree against
  uniform octree, on the same setup.

### Conclusion

Static refinement on the existing mesh is the cheap win: about 1.6× on the DK with nothing new to trust,
and it combines with the edge corrections (which allow coarser cells in the first place). The octree
delivers the predicted cell savings (3-7× against a graded mesh) with small interface errors, but per
unknown it is 3× more expensive on the CPU and would be gather-bound on the GPU. It pays only on boards
where the tensor-product waste is extreme (sparse copper in a large domain), and only after a full GPU
implementation. The recommendation stands: rectilinear mesh, edge corrections, static refinement;
`no_slivers` and `z_grade` are worth making the default once the fast-mode calibration is re-run.

## Edge corrections and antenna efficiency: the MTH R7A check (Nov 2026)

A PIFA run (MTH R7A, C20 one-port, `"extent": "board"`, 0.05 mm ROI with cells graded to 1.48 mm around
it, wideband FR4 loss, `"mesh": {"no_slivers", "z_grade": 3}`) lost about 8 points of radiation efficiency
at 2.40 GHz with `"edges"`, and a 0.035 mm mesh took it back. The discriminating runs (radiation
efficiency at 2.40 / 2.44 / 2.48 GHz, port stop at -60 dB):

| run (0.05 mm unless noted) | 2.40 | 2.44 | 2.48 GHz |
|---|---|---|---|
| no edges, same mesh policy | 89.5 | 89.1 | 88.0 % |
| no edges, without the `loss.copper` stand-in | 90.6 | 90.3 | 89.2 % |
| edges, without `loss.copper` | 84.0 | 85.4 | 87.6 % |
| conformal only | 83.9 | 86.4 | 88.5 % |
| singular only | 89.1 | 89.0 | 88.3 % |
| edges, limited to the fine mesh (`max_cell`) | 89.9 | | |
| edges, limited, 0.035 mm | 91.8 | | |
| no edges, 0.035 mm | 90.4 | 89.2 | 87.2 % |

The mesh policy and the copper-loss term are not involved; the conformal fractions are. Of the
permittivity × volume they added, 98 % sat on graded 1.48 mm in-plane edges where ground pours end, just
outside the ROI, at factors up to 1/min_open = 6.7. The corrections are now limited to edges no longer
than 2.5 fine cells (`"edges"` / `max_cell`, `fieldsolve.md`).

The efficiency spectrum explains the size of the effect. With 21 far-field frequencies over 2-3 GHz and
the energy stop at -50 dB (3.3× the steps of the port stop), each mesh shows one narrow lossy board
resonance (late-time Q 70-80): uncorrected at 2.25 GHz (81 % against 92 % around it), with the unlimited
corrections at 2.40 GHz (77 %). Away from it the port and energy stops agree within 1-2 points; the
energy stop only sharpens the dip. With `"conductor": true` (surface-impedance copper) instead of the
stand-in, limited corrections and none agree within 1.5 points over 2.40-2.50 GHz (84.6 / 82.1 / 82.4 %
and 84.8 / 83.5 / 81.0 % at 2.40 / 2.45 / 2.50), both with the dip at 2.25 GHz. The resonance lives on the
coarse board mesh, so its frequency is uncertain by several percent either way: single-frequency
efficiencies near a band edge are better read from a dense frequency list, and a run whose field energy
still decays slower than a few dB/ns at the stop has such a mode.
