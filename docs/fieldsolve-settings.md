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
| `a high-Q mode is still ringing at the stop` (far-field runs) | When the run stopped, a mode with a Q above 30 was still ringing (at 2.4 GHz, field energy falling slower than about 2.2 dB/ns). An antenna's own loaded ringing, with a Q of 10-20, is not flagged. The warning gives the mode's frequency and Q. Efficiency within about f/Q of that frequency depends on the mesh. | List far-field frequencies 0.02-0.05 GHz apart across the band, quote band averages rather than single frequencies, and refine the mesh (`max_cell`, a `"refine"` box) where the mode lives. A field export at its frequency shows where that is. |
| `N signal nets ... carry 10 mm or more of track with nothing modelled on them` (far-field runs) | Long signal tracks near the fine region end at a chip that isn't modelled. Their far ends are often unconnected pins too, so they are resonant stubs. | `"pin_loads"` on the chip they run to (see The model). Nets you have checked go in `"open_nets_ok"`. |
| `port X reached max_steps before any stop rule` | The run was cut off before it converged, so the S-parameters are an estimate. | Raise `max_steps`, or loosen `end_db` |
| `power balance ... check the port's placement` (on a full run) | The port's energy does not add up. The port bridge may miss a pad, or short the port. | Check the printed nets under each port, and `fieldsolve raster` |
| `edit ... (net X) overlaps copper of net Y` | An `"add"` edit carries a net and lies on another net's copper. The net-aware raster cuts copper of two nets apart where they touch, so the added copper is **not** connected to it. | Leave `"net"` out of edits that must join existing copper (an inverted-F's short into the ground pour). Net-less copper never conflicts. |

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
  - With an explicit `"dielectrics"` stack-up, the mask used to be silently dropped. It now adds a layer
    over each outer copper layer: `"solder_mask_t"` (mm, default 0.02), `"solder_mask_er"` (3.8) and
    `"solder_mask_tand"` (0.025).
  - A 20 µm layer forces fine z cells and costs about 5× the run time. An equivalent 40 µm layer with
    `"min_cells_layer": 1` costs about 2.3×. Keep (εr − 1)·t and εr·tan δ·t equal: for 20 µm of 3.8 /
    0.025, that is 40 µm of 2.4 / 0.0198.
  - On the MTH PIFA, the mask mostly detunes the antenna (2.44 GHz feed impedance 49 → 40 + j8 Ω); its own
    loss is small.
- **Decoupling capacitors near the RF section:** model them as `"elements"` shorts. Without them, the
  supply copper floats at RF. Short only parts the production BOM fits: a footprint with a value in
  KiCad may be unfitted. On MTH R7A, C6, C7 and C23 have values but aren't fitted, and the first MTH specs
  shorted them anyway.
- **Tracks to unmodelled ICs:** a track ending at a chip pin is open-ended in the model, which makes it a
  resonant stub (a quarter wave is about 17 mm at 2.4 GHz on FR4 microstrip). Near an antenna, such a
  stub can take efficiency in a narrow band, as on MTH. Load the chip's pins:

  ```json
  "pin_loads": [{"ref": "U12", "c": 3e-12}]
  ```

  - Each pin of `ref` gets a lumped capacitor from its pad to the next copper layer below. You can also
    give `"r"` in ohms, alone or in parallel with the capacitor.
  - Pins on GND or VSS nets are skipped, and so are pins already used by a port or element.
  - **Pins whose net is already tied to ground by a modelled short are skipped too,** for example a DEC
    pin or a supply pin with its decoupling capacitor modelled as a short. A load there closes a tiny
    lossless loop (pad, short, vias, load) that rings on its own. On MTH, those loops rang at 2.15 GHz
    with a Q of 20-35 and inflated in-band efficiency by 3.5 points. `"load_decoupled": true` loads them
    anyway.
  - `"pins"` / `"except"` restrict the set, and `"to"` names the layer to load against.
  - The loads sit on the nearest mesh node, so they add no mesh lines.
  - The run prints how many pins were loaded, and how many land on no ground copper (those loads float).
  - Use the chip's pad capacitance from its datasheet: about 3 pF for nRF52840 GPIOs.
- **Enclosures and off-board antennas** (`"bodies"`, `"layers"`; see [fieldsolve.md](fieldsolve.md#usage)):
  - Set `"f_loss"` to the band centre (e.g. 2.44e9). Bodies take a constant conductivity from `tand` at
    `f_loss`, not the wideband model.
  - Use the enclosure's STL in the board's frame and check the frame first: with `"y_up"`, the CAD default,
    STL y is −(KiCad y). Each body prints how many cells it fills, so a zero or a surprising count means
    a wrong frame or offset. Plot a section through the model before a long run.
  - Write down where εr and tan δ came from. Plastic datasheets usually quote 1 MHz values.
  - Model the board under an off-board antenna fully: decoupling shorts, pin loads and acknowledged open
    nets, as for a PCB antenna. Run the antenna in free space, with the enclosure only, and with the board
    only, as well as installed, so the cause of any detuning or loss is visible.
  - An antenna's cable matters. Model the shield's outer surface (posts plus a strip on an extra layer)
    and compare with a run without it.

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
  - **Now handled:** far-field runs warn about such nets, and `"pin_loads"` loads the chip's pins. On
    MTH, loading U12's 59 signal pins with 3 pF (the 12 decoupled pins skipped) removed the dip: at
    2.25 GHz the antenna alone went from 82.4 % to 88.9 %. In band (2.40-2.48 GHz), the mean moved from
    83.8 % to 82.4 %, so it stayed within about 1.5 points. The feed impedance moved by under 2 Ω, and the
    run cost nothing extra.
  - **What still rings after that:** the antenna's own decay (Q about 13-19 at 2.43 GHz). There is also a
    weaker component (Q about 25 at 2.47 GHz) around the VCC3V3 via beside C30, a decoupling capacitor
    modelled as an ideal short. It is probably current circulating in that short's loop; not confirmed.
    Neither shows up as a dip in the efficiency spectrum at 0.05 GHz spacing. Turning the supply nets'
    capacitors outside the ROI into shorts (`"decouple"`) changed nothing.

## Fixed issues (keep in mind when reading older results)

- **An antenna short drawn with its own net (warned since Nov 2026).** Two eBUS inverted-F studies drew the
  whole antenna, short included, as `"net": "ANT"` over the GND pour. The raster cut the short off, so the
  antenna was an open-ended monopole: capacitive below resonance, with 3-9 Ω at a series resonance. An
  inverted-F with its short connected is inductive below resonance and shows a high-resistance parallel
  resonance. The edit now warns, and the eBUS results were re-run.
- **`"edges"` on coarse cells (fixed Nov 2026, 58f602d).** Conformal fractions were applied on the coarse
  graded cells outside the ROI. On the MTH PIFA at 0.05 mm, they moved the track resonance above onto
  2.40 GHz and cut efficiency there by about 8 points. Edge factors now apply only to edges no longer than 2.5×
  the fine cell (`"edges": {"fine": 2.5}`, or an absolute `{"max_cell": mm}`). `"edges": true` needs no
  change.
- **Copper loss through the plane-pair stand-in.** On the MTH PIFA, real copper (`"conductor": true`)
  costs 4-5 points of radiation efficiency compared with `"loss": {"copper": true}`. Every far-field
  result from before Nov 2026 that used the stand-in is that much optimistic.
