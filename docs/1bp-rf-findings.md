# 1BP RF frontend: findings from fieldsolve (Nov 2026)

1BP_mini R6B: nRF52840 (aQFN73), RF pin H23, a four-part ladder C21 1.0 pF shunt, L4 4.7 nH series, C18
1.2 pF shunt, L3 2.2 nH series, then the u.FL J1. These are exactly the values of Nordic's
four-component reference, configuration 7. The question was whether the 1BP layout of that network is
faulty compared with Nordic's.

## Summary

- **One real deviation from Nordic's layout: how C21 is grounded.** Nordic grounds the first shunt
  capacitor (their C3) only through its own trace, back to the chip's VSS_PA ball, alongside the ANT
  trace. Nothing else connects to its ground pad. 1BP floods C21's ground pad into the top ground pour.
- **That costs about 12-15 dB of 2nd-harmonic suppression in the model.** With identical parts and
  identical modelling, 1BP as built conducts and radiates about 20 dB more H2 than Nordic's reference.
  With Nordic's C21 routing transplanted from their Gerbers, the gap drops to 4-7 dB. The rest is the
  remaining layout (L4/C18/L3, the u.FL, a 0.8 mm board against 1.6 mm).
- **The H2 conclusion is robust.** It holds whether the chip's harmonic output is taken as a 50 ohm
  source, a voltage source or a current source.
- **The H3 numbers are not reliable for the tight-loop layouts, and are probably pessimistic** (see
  below). Don't use them to choose between layouts.
- **The in-band match cannot be judged from the model yet.** The model reproduces the shape of Nordic's
  measured network but puts its low-pass cutoff about 10-15 % low (see Validation).
- **The inner-layer keep-out under the ladder is not a fault.** Nordic's reference does the same.
- **Recommendation for the next revision:** copy Nordic's C21 routing exactly: the isolated ground pad
  and the return trace to F23 alongside the ANT trace (figure below). Then check on hardware:
  conducted H2/H3 at J1 with a spectrum analyser at the TX power you use (Nordic measured about -58 to
  -61 dBc at H2), and the network from the ANT pad with a VNA, chip removed, as Nordic measured theirs.

## What was compared

All on the whole board, from Gerbers. Every case uses the same fitted parts as measured Touchstone models:
Murata GRM1555C1H1R0 (1.0 pF), LQG15HS4N7B02 (4.7 nH), LQG15HS2N2B02 (2.2 nH) and Yageo
CC0402BRNPO9BN1R2 (1.2 pF). Nordic's BOM names only values. The cases:

1. **Nordic reference, configuration 7**, from `nrf52840-qiaa-reference-layout-1_3.zip`. Altium Gerbers
   with no nets: ports placed by coordinate, the 13 decoupling capacitors as 30 mOhm resistors, the
   stack-up taken from the `.PcbDoc`.
2. **1BP as built**: C21's ground flooded into the pour.
3. **1BP, hand-made tie**: C21 mirrored, its ground pad isolated, and a hand-drawn 0.15 mm trace to F23.
4. **1BP with Nordic's C3 routing transplanted.** Nordic's top copper from the ANT ball to C3 is rotated
   onto 1BP. The two boards share the package and the feed handedness, so Nordic's VSS_PA ball lands
   exactly on 1BP's F23. Over the area from the chip edge to just past C3, 1BP's pour and ground vias are
   removed (Nordic has none there). 1BP's other nets (DEC3/DEC4/DEC5, the DEC4 via) are kept.

![Transplant (left) against Nordic's copper rotated into 1BP's frame (right)](img/1bp_nordic_c3_transplant.jpg)

Model settings, the same for every case:
- 0.05 mm mesh along the RF path, graded outside it;
- wideband FR4 (Debye fit, tan delta 0.015, flat over 0.3-10 GHz);
- far field from a Huygens box in 12-15 mm of air;
- each part a measured model terminating its own EM port (`fieldsolve cosim`);
- harmonics given per watt into the chip port.

## Results

| per W into the chip port | Nordic ref | 1BP as built | hand tie | Nordic routing transplanted |
|---|---|---|---|---|
| H2 radiated / to connector | -49.3 / -63.8 dB | -29.7 / -41.8 dB | -37.7 / -53.3 dB | -41.9 / -58.1 dB |
| H3 radiated / to connector | -18.6 / -16.7 dB | -25.2 / -39.4 dB | -15.6 / -24.5 dB | -14.8 / -22.1 dB |
| chip-side impedance at H3 | 3.1 + j3.6 ohm | 5.5 + j91 ohm | 8.4 + j0.6 ohm | 11.2 + j1.5 ohm |
| chip-side impedance at 2.44 GHz | 12.8 - j45.5 ohm | 45.6 - j52.8 ohm | 22.9 - j50.7 ohm | 23.4 - j55.4 ohm |

H2 relative to the reference holds under each assumption about the chip's harmonic output. Radiated
per unit pin voltage: as built +19 dB, transplant +6 dB. Per unit pin current: +20 dB and +8 dB.

![H2 (top) and H3 (bottom), shared scale per row: Nordic reference, 1BP as built, hand tie, transplant](img/1bp_h2h3_four_layouts.jpg)

### Why the H3 numbers are unreliable

The reference, the hand tie and the transplant all have a short C3 ground loop. Its series resonance
falls near 7.3 GHz in the model, so the chip sees close to a short at H3. Above that frequency the
shunt capacitors are past their self-resonance, and the ladder's transfer rises steeply:

| chip -> connector, dB | 6.6 GHz | 6.9 | 7.1 | 7.32 | 7.5 | 7.7 | 7.9 |
|---|---|---|---|---|---|---|---|
| Nordic reference | -45.1 | -29.3 | -22.4 | -16.6 | -12.3 | -9.8 | -6.6 |
| 1BP as built | -46.0 | -41.7 | -39.9 | -39.4 | -40.4 | -42.1 | -36.0 |
| transplant | -42.9 | -31.5 | -24.9 | -22.0 | -20.7 | -17.5 | -15.2 |

- **That is 4-5 dB per 100 MHz at H3.** A frequency-scale error of a few percent moves the H3 result by
  tens of dB.
- **The model's frequencies run 10-15 % low** against Nordic's measurement (Validation). If the real
  features sit that much higher, the board at 7.32 GHz behaves like the model near 6.4-6.6 GHz: 25-30 dB
  better.
- **Losses are missing.** There is no copper loss in the traces, and the reference run had the plane
  copper-loss term off. Real losses blunt the resonance further.
- **"Per W into a 50 ohm port" exaggerates a near-short load.** The chip's real behaviour at H3 is
  unknown.
- **The measured data agree.** Nordic measured -48 to -51 dBc of H3 at their connector, which is far
  lower than the model's reference transfer would imply.

1BP as built does not sit on such an edge (flat around -40 dB), so its H3 figure is the trustworthy one.

## Validation against Nordic's measurement

Nordic's matching-network data (`nRF52840 QIAAC0 matching network`, 2018) holds a measured two-port, chip
side to antenna side. Its harmonic rejection (H2 -32.5 dB, H3 -28 dB) points to the three-part network,
so configuration 1 was simulated with both of Nordic's value sets and compared on magnitudes. Magnitudes
do not depend on the measurement's offset reference planes.

![Simulated configuration 1 (v1.0, v1.1) against Nordic's measured network](img/nordic_config1_vs_measured.png)

- **The shape matches:** low-pass, the H2 null, total reflection at the harmonics.
- **The simulated cutoff sits about 10-15 % low.** At 2.44 GHz the model gives S11 -3 dB and S21 -3.7 dB,
  against -7.8 dB and -1.5 dB measured.
- **Likely causes:**
  - the measured board is not exactly this layout: the measurement carries about 400 ps of extra cable or
    trace delay;
  - zero-thickness copper edges on the mesh make conductors act slightly wider, which narrows Nordic's
    0.15-0.25 mm coplanar gaps and adds capacitance.
- **Next step:** a mesh-refinement run, 0.025 mm on configuration 1, would separate the two causes.

## Other findings along the way

- **In1 is the VCC2V4 plane.** Referencing a port to it, without its decoupling capacitors modelled, made
  the plane ring as a patch: spurious notches at 2.32, 2.95 and 7.55 GHz.
- **Nordic's 115.8 - j33.4 ohm is not the pin load.** It is a hot impedance measured from the connector
  with the chip transmitting. The load the chip is designed for is in Nordic's `.s2p` file (35 - j35 ohm
  at 2.44 GHz), within the uncertainty of their offset reference plane.
- **The board radiates more harmonic than reaches J1,** in every case, by 9-16 dB at H2. A radiated
  pre-scan is worth doing as well as the conducted check.

## Reproducing

The inputs and scripts were in the session scratchpad. To rebuild:

- **Models:** specs in the style of `1bp_board_ff_w2.json`, run with `fieldsolve run`, then
  `fieldsolve cosim` with `"farfield": true`.
- **Nordic reference:** the Gerbers plus a hand-written `.gbrjob` carrying the `.PcbDoc` stack-up; ports
  by coordinate.
- **Transplant:** Nordic's top layer dumped with `fieldsolve prims`, rotated by
  (u, v) -> (149.419 - v, -92.8695 + u) from Nordic's ANT ball at (2.751, 1.001), clipped, and appended
  to 1BP's top-layer Gerber as regions with X2 nets.
- **Figures:** `tools/blender_farfield.py` and `tools/render_farfield_panels.py`.
