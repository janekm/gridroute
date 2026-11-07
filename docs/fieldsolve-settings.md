# fieldsolve settings: what to use, and the traps

This is the page to read before writing a `fieldsolve run` spec. [fieldsolve.md](fieldsolve.md) describes
the model behind each key. This page says which keys to set for which job, and which mistakes have already
cost us results.

## Start from a preset

```json
"preset": "antenna"
```

A preset fills in the recommended model settings that the spec leaves out. A key the spec sets itself
always wins, and `result.json` lists what the preset filled in (`"settings"`).

| key | `"rf"` | `"antenna"` | why |
|---|---|---|---|
| `"loss"` | `{"model": "wideband"}` | same | FR4's tan δ is flat with frequency. Without this, tan δ is right only at `f_loss`, and resonances ring 4-8× too long at 4-8 GHz. |
| `"conductor"` | `true` | same | Copper gets a surface impedance (skin effect, thickness, 1 µm roughness). Without it, copper is lossless. |
| `"edges"` | `true` | same | Sub-cell copper edges (conformal and edge-singularity corrections): a 2-4× coarser cell gives the same accuracy. |
| `"mesh"` | `{"no_slivers": true, "z_grade": 3.0}` | same | No half-size cells next to fixed lines (they set the time step for the whole board), and z cells grade through thick cores. About 1.3-1.5× less work at the same accuracy. |
| `"extent"` | (cropped to `roi` + `air`) | `"board"` | On a small board the whole ground plane radiates, so an antenna run needs all of it. |

- **`"rf"`** is for conducted problems in a cropped window: matching networks, transitions, filters,
  decoupling.
- **`"antenna"`** is for anything that radiates or whose far field you want.

Neither preset sets the geometry keys (`roi`, `cell`, `air`), the band, the ports or the stop rule. Those
are covered below.

## Read the model line and the warnings

Every run prints the model it actually uses before it starts, for example:

```
fieldsolve: model: preset antenna | dielectric: wideband tan d (...) | copper: surface impedance (sigma 5.80e7 S/m, roughness 1 um) | edges: conformal + singular | mesh: {...} | fine region [...] at 0.05 mm, extent board
```

It then prints `WARNING:` lines for known traps. The same text is in `result.json` under `"settings"` and
`"checks"`. Do a 20-step dry run first (`"max_steps": 20`, `"excite": "first"`), read these lines, and only
then start the long run.

| warning | what it means | fix |
|---|---|---|
| `"loss": {"copper": true} is a plane-pair stand-in` | The copper stays lossless. Instead, about δ_s/h is added to the dielectric's tan δ. That is a rough stand-in for loss between two planes, and wrong for traces and antennas. | `"conductor": true`, and drop `"copper"` from `"loss"` |
| `far-field run with lossless copper` | Radiation efficiency and gain come out optimistic. | `"conductor": true`, or `"preset": "antenna"` |
| `far-field run without "loss": {"model": "wideband"}` | Substrate loss is right at one frequency only. | `"loss": {"model": "wideband"}` |
| `far-field run with the domain cropped to "roi"` | The board's ground beyond the window is missing, and on a small board it is part of the antenna. | `"extent": "board"` |
| `far-field run without "edges"` | Copper edges are staircased. | `"edges": true` |
| `narrow copper ... in cells coarser than ...` | Traces, meanders or antenna arms outside the fine region sit in cells several times larger than `cell`. Their length and coupling are then wrong. | Extend `roi` over anything that carries RF current. After you have checked the rest, list it in `"coarse_copper_ok"`. |
| `a high-Q mode is still ringing at the stop` (far-field runs) | When the run stopped, the field energy was still falling slower than 5 dB/ns, so a board mode with a Q above about 30 is still ringing. The warning gives the mode's frequency and Q. Efficiency within about f/Q of that frequency depends on the mesh. | List far-field frequencies 0.02-0.05 GHz apart across the band, quote band averages rather than single frequencies, and refine the mesh (`max_cell`, a `"refine"` box) where the mode lives. A field export at its frequency shows where that is. |
| `port X reached max_steps before any stop rule` | The run was cut off before it converged, so the S-parameters are an estimate. | Raise `max_steps`, or loosen `end_db` |
| `power balance ... check the port's placement` (on a full run) | The port's energy does not add up. The port bridge may miss a pad, or short the port. | Check the printed nets under each port, and `fieldsolve raster` |

## Geometry: the fine region

- **`roi` must cover every piece of copper that carries RF current.** That is the antenna (including the
  open end of the radiator), the feed, the matching parts and the chip's RF pin. Outside `roi`, cells grade
  up to `max_cell`. A 0.4 mm meander in 0.4 mm cells is roughly the right length, but not the right shape.
  - **MTH R7A, Nov 2026:** the ROI started 3 mm short of the PIFA's open end. The resonance came out
    about 3 % low, the antenna looked inductive, and the recommended ladder fix was wrong. The narrow-copper
    warning catches this now.
- **Check the ROI against the board drawing** (`fieldsolve svg gbr/ board.svg`), not against component
  positions.
- **`"refine"` boxes** add fine cells in more places without widening the whole ROI. Use them for a remote
  feed, a shorting stub or a gap.
- **`cell`:** 0.05 mm with `"edges"` is accurate for 0402 parts and 0.1-0.2 mm traces on 0.1-1.6 mm
  substrates at 2.4 GHz.
  - Exception: thin prepreg under a trace (75 µm on the nRF52840 DK) still needs finer cells across that
    trace.
  - For a number you will publish or act on, do a two-mesh check: the run again at about 0.7 × `cell`.
- **`air`:** 15 mm (about λ/8 at 2.4 GHz) is what the MTH and DK far-field runs used. The far-field box
  sits `farfield.gap` cells inside the air margin.
- **`"extent": "board"`** for antennas (the `"antenna"` preset sets it).

## The model

- **Dielectric:** use `"dielectrics"` when the `.gbrjob` gives no stack-up. Write down where εr came from:
  an assumed εr shifts the resonance as its square root.
- **Copper loss:** use `"conductor"` and never `"loss": {"copper": true}`. The surface-impedance model
  assumes current on one face. That is right for a trace over a plane, and up to 2× too lossy for an
  isolated strip carrying current on both faces, such as an antenna arm over cleared ground.
  `"roughness_um"` defaults to 1 µm, which is standard electrodeposited foil.
- **Solder mask** (`"solder_mask": true`) lowers an antenna's resonance. Leave it out only when comparing
  with a model that also leaves it out.
- **Decoupling capacitors near the RF section:** model them as `"elements"` shorts. Without them, the
  supply copper floats at RF.
- **Tracks to unmodelled ICs:** a long track that ends at a chip pin is open-ended in the model, and so a
  resonant stub (a quarter wave is about 17 mm at 2.4 GHz on FR4 microstrip). Near an antenna it can take
  efficiency in a narrow band, as on MTH (see Known open issues). Load such pins with `"elements"` (a few
  pF to ground) when the tracks are longer than about 10 mm and lie within a few cm of the antenna.

## Ports, band and stopping

- **Ports:** one per part position you want to vary later (`"excite": "all"`), so that `fieldsolve cosim`
  and `tune` can place any part without a new EM run. Check the printed nets under each port. A
  rotated-part port gets a staircase bridge, and both pads must be reached.
- **Band:** `fmax` at least 3.2× the operating frequency if you want H2 and H3. `nf` 300-400.
- **Far field:** list frequencies on a dense grid across the band, 0.02-0.05 GHz apart, plus the
  harmonics. A board can carry a narrow, lossy resonance about 0.05 GHz wide that removes 5-12 points of
  efficiency at one frequency and moves between meshes. Three spot frequencies (2.40 / 2.44 / 2.48 GHz) can
  land on such a dip, or miss it, by chance. Quote band averages.
- **Stopping:** `"end_db": -60` with the default port stop is the usual choice. A weakly coupled board
  resonance can make one port wait 2-4× longer: that is expected, not a fault. Stopping on field energy
  instead (`"stop": "energy"`) took 3.3× the steps on the MTH board, but moved broadband efficiency by only
  0.3-0.5 points. It sharpens a narrow mode's dip rather than removing it, so it isn't the default. For searches,
  `"objective": {..., "extrapolate": true}` with `"stop": "objective"` stops as soon as the extrapolated
  objective is stable ([fieldsolve-fast-mode.md](fieldsolve-fast-mode.md)).
- **Fidelity levels:** `--fast` (screening) is for searches only. Report final numbers from a plain or
  `--fidelity validated` run.

## After the run

- `result.json` `"settings"` is the model that was used, and `"checks"` the warnings. Quote both when
  writing up results.
- **Port power balance** should be about 0 per port over the band.
- **Far field:** the box flux check should agree with the far-field integral to about 1 %.
- **Timing:** two jobs sharing the GPU each run 1.5-2× slower. Compare speed-ups in cells × steps, or time
  them on an idle GPU.

## Known open issues

- **Open-ended tracks at unmodelled ICs make resonant stubs.**
  - On the MTH board, two signal tracks run side by side from the nRF52840 to header H5: P1_04 (0.15 mm
    wide, 24.7 mm long) and SWDCLK (about 32 mm). The chip isn't modelled, so both are open at the U12 end,
    and together they resonate near 2.25-2.40 GHz with a Q of about 70-80. The resonance takes 5-12 points
    of antenna efficiency in a narrow band.
  - The tracks run along y where the y cells are 1.48 mm, so their frequency also depends on the mesh.
  - On the real board, the pin capacitance moves the resonance, so the dip is partly a model artefact.
  - Not yet automated. For an antenna run, load the pins of unmodelled ICs that long tracks near the
    antenna lead to, with `"elements"` (a few pF to ground, or an R), or check where the ringing
    warning's mode lives with a field export at its frequency.

## Fixed issues (keep in mind when reading older results)

- **`"edges"` on coarse cells (fixed Nov 2026, 58f602d).** Conformal fractions were applied on the coarse
  graded cells outside the ROI. On the MTH PIFA at 0.05 mm, they moved the track resonance above onto
  2.40 GHz and cut efficiency there by about 8 points. Edge factors now apply only to edges no longer than 2.5×
  the fine cell (`"edges": {"fine": 2.5}`, or an absolute `{"max_cell": mm}`). `"edges": true` needs no
  change.
- **Copper loss through the plane-pair stand-in.** On the MTH PIFA, real copper (`"conductor": true`)
  costs 4-5 points of radiation efficiency compared with `"loss": {"copper": true}`. Every far-field
  result from before Nov 2026 that used the stand-in is that much optimistic.
