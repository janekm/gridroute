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
| F1 | **Touchstone reader and network connection** in `post.rs`: N-port S with complex reference impedances, connect/terminate, power-wave renormalisation, mixed-mode (differential) conversion | Co-simulation core. The chip's RF pin impedance from the datasheet becomes a complex reference. Differential balun ports need mixed-mode. | M |
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
    bottom/top GND.
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

1. **Ladder as copper only.** Ports at C21 (ANT–GND), L4, C18, L3 and J1.4–GND: a 5-port EM run. Then
   connect the Murata LQG15HS/GRM1555 S-parameters (needs F1) and get S chip → J1, checked against the
   ideal ladder.
2. **u.FL pad study.** The excess capacitance of J1.4 with and without the inner-layer keep-out, shown as
   TDR and as S11 at 2.44 GHz.
3. **Antenna on the pipe** (A1–A4): the FPC antenna's own model is not in our layout. Simulate it as a
   generic 2.4 GHz FPC dipole, or ask the vendor for its geometry, plus the board, coax and the pipe.

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
