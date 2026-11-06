# fieldsolve: what radio frontends and antennas need

This is the plan for taking fieldsolve from "loop inductance and S-parameters of bare copper" to:

- **RF frontends:** full simulation of our radio frontends, such as the ESP32 balun on EXT, the nRF54
  frontend filters with their layout, and the u.FL footprint;
- **antennas:** PCB antennas in free space, in an enclosure, and next to a pipe (1BP).

It builds on what exists today (see [fieldsolve.md](fieldsolve.md)): Gerber import, layered model, lumped
ports and R/C/short elements, CPML, CPU + Metal stepping, Z/S/Touchstone output, field snapshots and the
openEMS cross-check.

## 1. Radio frontends

A frontend is mostly small parts on short traces: a balun or IPD, LC filters, matching components, ESD, and
the connector. Its performance is set by the parts and by the copper between them. That copper includes the
pads' capacitance to the planes, the via and return-path inductance, coupling between neighbouring
components, and the connector footprint. The plan is EM–circuit co-simulation, as Sonnet/ADS-style flows do
it:

1. One fieldsolve run of the board with a port at every component pad and at the chip and connector
   reference planes. The parts are removed and the copper is left as drawn. This gives the board's N-port
   S-matrix, including every parasitic and coupling path.
2. Vendor S-parameter models of the actual parts are then connected to that N-port in post-processing. This
   covers the balun (e.g. a Johanson 2450BM.. or an ST/Murata IPD .s3p), the filters, and Murata/TDK C and L
   models, which carry real ESR, ESL and self-resonance. Connecting networks is linear algebra (sub-network
   connection), so it takes milliseconds per frequency.

What this buys:

- **Correct physics.** The vendor model carries the part's parasitics and the EM run carries the
  layout's.
- **Instant tuning.** Swapping a matching capacitor from 1.0 to 1.2 pF, or comparing two balun part numbers,
  needs no new EM run. Matching values can be optimised over E-series against the vendor models.
- **The same flow for every radio.** ESP32 (balun on EXT) and nRF54 (filter or matching network) use the
  same machinery; only the parts list and reference impedances differ.

### Work items (frontends)

| # | Item | Why | Size |
|---|---|---|---|
| F1 | **Touchstone reader and network connection**: N-port S with complex reference impedances, connect/terminate, power-wave renormalisation, mixed-mode (differential) conversion. *Done (Nov 2026), in `cosim.rs` / `fieldsolve cosim`:* Touchstone reader, two-terminal parts as port terminations (ideal, parasitic, `s2p`/`s1p`), ideal-ladder comparison, complex reference renormalisation. *Still to do:* multi-terminal parts (a balun `.s3p` connected port to port) and mixed-mode conversion. | Co-simulation core. The chip's RF pin impedance from the datasheet becomes a complex reference. Differential balun ports need mixed-mode. | M |
| F2 | **Ports from parts**: `"components": ["L3", "C7", "FL1", "U4"]` puts a port on every signal pad, referenced to the nearest GND pad or the plane underneath. Two-pad parts get an in-plane gap port; balun and filter GND pins are shorted to their pads/vias. Ports on RF nets are found the way `bench/pcb/topor_rf.py` already finds them. | Writing 20 port definitions by hand does not scale. | M |
| F3 | **Port accuracy**: resolve the reciprocal-port discrepancy (single-node vertical columns act as zero-radius thin wires; the gap port's current spreads differently). Plan: ports with a finite footprint matching the pad, a thin-wire correction, and an openEMS cross-check per port type. | At 2.4 GHz, 0.1 nH is 1.5 Ω of reactance: the matching network sees it. | M |
| F4 | **RF meshing**: on RF nets, mesh lines follow the copper edges with the thirds rule (the microstrip case showed 6 % εeff error without it); optional two-plane thick copper (35 µm) for gaps ≤ 0.15 mm; a sub-cell solder-mask layer (εr ≈ 3.5 over the traces, folded into edge averaging, no extra mesh plane) | Without these, line impedance and coupling are off by a few percent. CPWG and filter gaps are sensitive to all three. | M |
| F5 | **Conductor loss**: a surface-impedance (skin-effect) sheet model for copper, fitted over the band, plus roughness factor | Filter insertion loss, frontend loss budget, antenna efficiency. Copper is lossless now. | M |
| F6 | **TDR output**: an impedance profile versus distance from the S11 step response | The clearest way to see the u.FL pad's excess capacitance, a via stub, or a neck-down | S |
| F7 | **u.FL connector**: (a) footprint only. Port at the pad, with the plane cut-out under the pad as drawn; the pad's excess C appears in S11 and TDR. That already answers "does the void under the pad help?" (b) The connector body. 3-D metal from STL (shell ring, centre contact) and a **radial coaxial port** at the mating plane: lumped edges arranged around the centre pin, 50 Ω coax reference. Optionally cascade Hirose's mated-connector S-parameters. | The footprint parasitic is the main u.FL effect. The body matters above about 3 GHz and for the harmonics. | (a) S, (b) M |
| F8 | **Aggressor coupling**: ports on DC-DC switch nodes, crystal and fast digital nets; report S21 into the frontend and antenna feed | Desense and spurs: the frontend's other job | S once F2 exists |
| F9 | **Frequency range for filters**: run to the 3rd harmonic (7.5 GHz), with the mesh at ≤ λ/20 in FR4 there (≈ 0.6 mm; the frontend ROI is finer anyway) | Harmonic rejection of the nRF54 filters and the ESP32 balun | config |

Outputs per frontend:

- S-parameters chip → u.FL, at the real reference planes, with parts;
- insertion and return loss over the band, harmonic rejection at 2f0 and 3f0;
- the optimised matching values with their sensitivity;
- a TDR profile of the RF path;
- field snapshots and field lines (already available) to show where coupling happens.

## 2. PCB antennas

### What the engine needs

| # | Item | Why | Size |
|---|---|---|---|
| A1 | **Near-to-far-field transform**: record the tangential E/H running DFT on a closed box just inside the PML (surface only, on Metal), then compute the far field. Outputs: 3-D pattern, directivity, gain, realised gain, radiated power, radiation efficiency, total efficiency, polarisation/axial ratio. Writes pattern files for plotting and for Blender (a coloured 3-D pattern mesh). | Antenna performance means pattern and efficiency, not just S11 | M–L |
| A2 | **3-D objects**: STL (and STEP via a converter) and primitives (box, cylinder, sphere) as metal or dielectric. Voxelised with fill-fraction ε averaging, so curved housings and pipes detune correctly without a fine mesh. Mesh hints come from object bounds. | Enclosure, pipe, battery, display, shield cans, the u.FL body, other boards | M |
| A3 | **Material library**: FR4 grades, solder mask, ABS/PC/PA housings (εr 2.6–3.2, tan δ 0.005–0.02), glass, water and wet soil, metals (σ) | Enclosure and pipe studies are only as good as these numbers | S |
| A4 | **Lossy and dispersive materials**: σ (done), plus Debye/Lorentz dispersion for water (εr ≈ 77 − j10 at 2.4 GHz, rising loss with frequency) and optionally a human-tissue phantom | A water-filled pipe, or a hand | M |
| A5 | **Larger domains**: per-edge material indices instead of four f32 coefficient arrays (48 → ~26 B/cell), CPML fused into the main kernels, f16 storage where safe | Free-space domains are 10–30 M cells; enclosures and pipes are bigger | M |
| A6 | **Antenna validation set**: half-wave dipole (D = 1.64, R ≈ 73 Ω), a probe-fed patch, and a vendor reference antenna with published efficiency (Nordic/Espressif reference designs), all cross-checked against openEMS's nf2ff on identical scenes | Trust in pattern and efficiency numbers | S–M |

### Scenarios

- **Free space.** The whole board outline with its ground plane (the board is part of the antenna).
  Matching parts are lumped or come from vendor models through F1. Air margin ≥ λ/4 (≈ 30 mm) to the NF2FF
  box, then PML. A rough size: 100 × 100 × 70 mm, 0.1 mm cells near the antenna graded to 3 mm, about
  10–15 M cells, 60–100 k steps, about 5–10 min on Metal. Matching-network optimisation needs no re-run
  (antenna port plus ports at the matching pads, F1/F2).
- **In an enclosure.** Add the housing STL (A2/A3), battery and display as metal, screws, and the other PCB.
  The usual effects are a resonance shift of 2–8 % downward, efficiency loss from housing tan δ and nearby
  metal, and pattern changes. Deliverables are S11, total efficiency and pattern with and without the
  housing, and re-tuned matching values.
- **Against a pipe (1BP).** The pipe as a cylinder (A2):
  - a metal pipe as a conductor (σ of steel or copper);
  - a plastic pipe as a dielectric shell around water (A4), where loss and detuning dominate.

  The pipe runs through the PML, which CPML handles. Sweep the device–pipe distance and orientation and the
  pipe diameter. Deliverables are detuning, efficiency loss and pattern blockage versus distance, which tell
  us whether the antenna needs a ground-plane spacer, a different orientation, or a matching network tuned
  for the installed case.

## 3. 1BP: what we know (Nov 2026)

- **Board:** 1BP_mini R6B is the current 1BP ("1 Pipe Boiler"). 1BP R6B has the same RF section on a longer
  board.
- **Radio and RF path:** nRF52840 (U1, AQFN-73), RF pin H23 (net ANT). All parts are 0402 on top, a ladder
  to an IPEX ANB02001 connector J1 (pin 4):
  1. C21 1.0 pF shunt
  2. L4 4.7 nH series
  3. C18 1.2 pF shunt
  4. L3 2.2 nH series

  There is no balun or ESD part.
- **Antenna:** an external FPC (YL1761-2.4, ≤ 1.19 dBi) on a 40 mm 1.13 mm coax, stuck inside the plastic
  shell. The device straps to a heating pipe.
- **Layout (KiCad coordinates, mm):**
  - parts: U1.H23 (149.42, 92.87), C21 at y 91.1, L4 at y 90.1–89.1, C18 at y 88.05, L3 at y 87.0–86.0,
    J1.4 (149.60, 84.72) with GND pads at y 83.2 and 81.68;
  - the feed is a ≈ 0.79 mm top-layer fill;
  - both inner layers are kept out under the matching network, so the feed is referenced to the
    bottom/top GND;
  - In1 is the `VCC2V4` power plane and In2 is GND. Reference ports to GND, not to In1 (see the results
    below).
- **Stack-up:** JLC 4-layer 0.8 mm:

  | layer | material | thickness | εr | tan δ |
  |---|---|---|---|---|
  | top prepreg | 7628 | 0.21 mm | 4.45 | 0.015 |
  | core | | 0.25 mm | 4.1 | |
  | bottom prepreg | 7628 | 0.21 mm | 4.45 | 0.015 |

  Copper is 35/15/15/35 µm. Pass the dielectrics in the spec as
  `"dielectrics": [{"t":0.21,"er":4.45,"tand":0.015},{"t":0.25,"er":4.1,"tand":0.015},{"t":0.21,"er":4.45,"tand":0.015}]`;
  the KiCad conversion says 1.6 mm, which is wrong.
- **Gerbers:** X2 with nets and pads from the KiCad conversion,
  `p-hw/easyeda2kicad/6B/1BP_mini/fabrication/gerbers/`; they load in fieldsolve. The EasyEDA production
  Gerbers carry no net attributes.
- **Validation data:** VNA plots of a separate chip-antenna test PCB, "1BP v260427"
  (`p-hw/docs/antenna/C293767/`). There is no 1BP enclosure CAD and there are no range tests.

First 1BP runs, in order:

1. **Ladder as copper only.** *Done (Nov 2026), results below.* The EM run has six ports:
   - CHIP: H23 to the GND ball F23;
   - C21, L4, C18 and L3 at their pads;
   - J1.4 to B.Cu.

   The parts are connected afterwards with `fieldsolve cosim`, using vendor S-parameters: Murata
   GRM1555C1H1R0CA01D, LQG15HS4N7B02D and LQG15HS2N2B02D, and Yageo CC0402BRNPO9BN1R2, which was
   measured as a shunt. The result is S chip → J1, checked against the ideal ladder.
2. **u.FL pad study.** The excess capacitance of J1.4 with and without the inner-layer keep-out, shown as
   TDR and as S11 at 2.44 GHz.
3. **Antenna on the pipe** (A1–A4): the FPC antenna's own model is not in our layout. Simulate it as a
   generic 2.4 GHz FPC dipole, or ask the vendor for its geometry, plus the board, coax and the pipe.

### 1BP results: matching ladder chip → J1 (Nov 2026)

Run settings:
- window x 146.5–153.5, y 80.5–94.5, plus a 3 mm margin; cell 0.05 mm; 6.2 M cells;
- six ports at 250 k steps each, 88 min in all on an M3 Max.

The decoupling capacitors in the window, C4, C8, C9, C3, C11, C13, C14, C17 and C24, are `short`
elements. They tie the In1 power plane and the DEC nets to ground.

Chip-side impedance looking towards J1, with J1 at 50 Ω:

| case | 2.44 GHz Zin | S21 at 2.44 GHz | S21 at 4.87 GHz | S21 at 7.33 GHz |
|---|---|---|---|---|
| ladder alone, nominal values | 65.9 − j30.4 Ω | −0.37 dB | −11.7 dB | −27.6 dB |
| ladder alone, vendor models | 65.8 − j30.3 Ω | −0.77 dB | −38.2 dB | −95 dB |
| layout + nominal parts | 53.7 − j36.3 Ω | −1.07 dB | −23.0 dB | −35.6 dB |
| layout + vendor parts | 42.5 − j47.6 Ω | −2.38 dB | −42.0 dB | −38.7 dB |
| layout, no parts | 48.0 + j1.6 Ω (S11 −31.7 dB) | −0.40 dB | −0.55 dB | −0.63 dB |

- **The feed copper on its own is a good 50 Ω line.** With L4 and L3 shorted and the capacitors open,
  S11 is −31.7 dB at 2.44 GHz.
- **Pads and vias pull the impedance the chip sees.** With the vendor parts it moves from 66 − j30 Ω
  to 42 − j48 Ω.
  - Against the ladder-alone impedance that would be |Γ| ≈ 0.27. But that is the wrong target: Nordic's
    measured network presents 35 − j35 Ω; see the re-tune below. Against that, the fitted design is at
    |Γ| ≈ 0.16–0.20 over 2.40–2.48 GHz, about 0.17 dB of mismatch.
  - Dissipated loss, from the power the network accepts to the power reaching J1, is 1.33 dB with the
    vendor parts. About 0.37 dB of that is the port-calibration floor; see the last point.
- **The third-harmonic rejection is limited by the layout.** C18's self-resonance near 7.4 GHz would
  give −95 dB as a bare part. On the board it gives −39 dB: the pad and via inductance move the series
  resonance, and the two ends of the ladder couple directly.
- **Port-calibration floor.** Every port loses a constant 1–7 % of power, the same from 100 MHz to
  8 GHz: 6.7 % at the floating L4 port, 1.8 % at CHIP. This is a port-extraction error, not physics. It
  shows as about 0.37 dB of "loss" on bare copper at 300 MHz. The largest eigenvalue of S^H S is 1.016,
  so the matrix is slightly non-passive. Fixing this is a follow-up.
- **First attempt, kept as a warning.** The chip port first ran vertically from H23 to In1, the
  `VCC2V4` plane, and the decoupling caps were not modelled. The floating plane rang as a patch:
  absorption notches at 2.32, 2.95 and 7.55 GHz, and S21 −1.4 dB on bare copper. `fieldsolve run` now
  prints each port terminal's net and warns about vertical ports with no GND end.


### 1BP re-tune against Nordic's measured network (Nov 2026)

**Target.** Nordic's matching-network data for the nRF52840 QIAAC0, dated 2018-10-03, contains
`nRF52840 QIAAC0 matching network.s2p`. It is their reference network measured as a two-port on their
board, port 1 on the chip side. Loaded with 50 Ω it presents:
- 38.1 − j34.9 Ω at 2.40 GHz;
- 35.2 − j35.5 Ω at 2.44 GHz;
- 32.4 − j35.7 Ω at 2.48 GHz.

This is the load the chip is designed for, Nordic's layout parasitics included. Their network's own
rejection is −32.5 dB at H2 and −28 dB at H3, and it dissipates 0.7 dB.

The same PDF shows impedances measured with the chip transmitting, e.g. 115.8 − j33.4 Ω at 0 dBm. Those
are hot impedances seen from the antenna connector into the network plus the running chip, not the chip
pin's load. A forum post reads 115.8 − j33.4 Ω as the pin impedance. Its "check" uses the nRF52832's
two-part network (0.8 pF shunt, 3.9 nH series: 111.8 − j32.9 Ω), a different chip at a different
reference plane.

**Search.** `tune_ladder.py` in the session scratchpad went through every combination of Murata
GRM1555C1H capacitors (0.5–2.5 pF) and LQG15HS inductors (1.0–6.8 nH) on the four pad sets, 194 k in
all, in 1.4 s. It reuses the one EM run, and every part uses Murata's measured S-parameters. Candidates
were scored on:
- the worst |Γ| against the target over 2.40–2.48 GHz;
- H2 and H3 rejection;
- the worst |Γ| over the 3⁴ tolerance corners (C ±0.1 pF, L ±0.1 nH, B tolerance).

| C21 / L4 / C18 / L3 | nominal \|Γ\| | worst-corner \|Γ\| | H2 / H3 (worst in band) |
|---|---|---|---|
| as fitted: 1.0 / 4.7 / 1.2 (Yageo) / 2.2 | 0.20 | 0.32 (Murata 1.2 pF) | −41 / −37 dB |
| C18 → 1.0 pF | 0.04 | 0.24 | −39 / −35 dB |
| **1.2 / 4.7 / 1.2 / 3.0** | **0.04** | **0.20** | **−39 / −55 dB** |
| same, Yageo C18 kept | 0.10 | | −41 / −65 dB |

**Recommendation: C21 1.0 → 1.2 pF (GRM1555C1H1R2BA01D) and L3 2.2 → 3.0 nH (LQG15HS3N0B02D). L4 and C18
stay.**
- Matching gains little: about 0.15 dB, because the fitted design is already close.
- The third harmonic drops from about −37 dB to about −55 dB in this model. This is where Nordic says
  the four-part network matters, at 5 dBm or more.
- The H3 figure rests on the EM model's coupling between the two ends of the ladder. Check it with a
  spectrum analyser before relying on it for approval margins.
- The calibration floor (about 0.37 dB) still sits in every loss figure. Port-model accuracy and the
  0.05 mm mesh limit the absolute numbers. The ranking between candidates, which share the same EM run,
  is more reliable than the absolute values.

## 4. Inputs we need

- **KiCad sources of EXT and 1BP** (Revision_6B in p-hw). Their prepared inputs were lost on 2026-11-04 and
  have to be re-exported. The nRF54 board's KiCad file is also needed.
- **BOM part numbers on the RF path:** balun, filters, matching L/C, ESD, u.FL and antenna parts, plus
  vendor S-parameter files (Johanson, Murata, TDK, ST and Hirose all publish them).
- **Radio pin reference impedances** from the ESP32 and nRF54 datasheets or hardware design guides
  (complex, per band).
- **Mechanical data:** enclosure STEP/STL and material; 1BP pipe diameter, wall, material and contents;
  mounting distance.
- **Measurements** if we have them: VNA S11 of built boards, conducted power at the u.FL, any chamber/OTA
  data. These calibrate the models. Antmicro's gerber2ems examples come with VNA data that serves as a
  public check of the whole flow.

## 5. Suggested order

1. **Port accuracy and RF meshing** (F3, F4) with openEMS cross-checks. Everything else depends on ports
   that are right to ~2 %.
2. **The co-simulation core** (F1, F2, F6). First deliverable: the EXT ESP32 frontend with vendor balun and
   matching models, S chip → u.FL, plus the u.FL footprint study (F7a).
3. **Conductor loss** (F5), then the nRF54 filters with harmonics (F9).
4. **NF2FF and the antenna validation set** (A1, A6), then free-space antenna runs.
5. **3-D objects and materials** (A2–A4), then enclosure and pipe studies. Memory work (A5) as the domains
   grow.
6. The u.FL body with a coaxial port (F7b) and aggressor coupling (F8) when the above is solid.
