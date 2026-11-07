# MTH R7A PCB antenna: findings from fieldsolve (Nov 2026)

MTH R7A (`p-hw/Revision_7/MTH`, KiCad conversion `p-hw/easyeda2kicad/7A/MTH`) is a 19.5 × 60 mm, 2-layer,
0.8 mm board. Its antenna, ANT1, is a meandered PIFA ("WCH_24G_PIFA_10mm_mod") across the top 6 mm, with
ground cleared on both layers below it. The antenna was left deliberately long, so it can be trimmed once
it is in the enclosure.

RF path, as populated:

- nRF52840 (U12, aQFN73, rotated 45°), ANT ball H23;
- C59 1.0 pF shunt, L4 4.7 nH series, C60 1.2 pF shunt, L5 2.2 nH series (Nordic's configuration 7
  values, the same as 1BP);
- R70 0 Ω, then the PIFA feed.
- **C20** is a shunt pad at the feed, not fitted.
- **The IPEX connector U24** hangs off R69, which isn't fitted, so it isn't in the RF path.
- **Parts:** the capacitors are Fenghua 0402CG and the inductors Murata LQG15HS.

## Summary

- **The antenna is good.** In free space it resonates just below 2.40 GHz, as intended for a design that
  will be trimmed, and radiates 82-88 % of the power it accepts across the band.
- **The ladder isn't matched to it.** The configuration-7 values were designed for a 50 Ω connector. On
  the PIFA (28 + j16 Ω at 2.44 GHz), the chip sees 94 − j21 Ω instead of its 35 − j35 Ω design load. The
  ladder dissipates 1.0-2.0 dB on top of the antenna's own loss. The worst in-band case is −4.6 dB, at
  2.48 GHz.
- **The cheapest fix: L5 → 0 Ω and C20 = 0.5 pF.** The worst in-band case improves from −4.6 to about
  −2.1 dB, and band centre from −3.1 to −1.4 dB, with H2/H3 unchanged. Full retunes gain only about 0.5 dB
  more.
- **H2 suppression is weak.** The radiated H2 is −25.5 dB per watt into the chip pin, against −29.7 dB for
  1BP as built and −49.3 dB for Nordic's reference layout. C59's ground pad is flooded into the top pour,
  the same deviation found on 1BP ([1bp-rf-findings.md](1bp-rf-findings.md)). The same fix should apply:
  an isolated return to VSS_PA. That fix hasn't been simulated here.
- **Retune again after trimming.** The trimmed antenna in the enclosure will present a different feed
  impedance. Re-run the EM model on the trimmed copper, plus the enclosure if there is data, and repeat
  the retune. The retune takes seconds; the EM run about 45 minutes.

![RF section, top copper (left) and bottom (right): the diagonal ladder from the chip to R70, the PIFA feed and short above](img/mth7_rf_copper.jpg)

## Model

- **The whole board** in free space, from the KiCad conversion's X2 Gerbers. 0.05 mm cells over the RF
  section and the antenna, graded outside, with far field from a Huygens box in 15 mm of air.
- **Stack-up:** one 0.73 mm FR4 core (εr 4.5, tan δ 0.015, wideband Debye fit) between 35 µm copper.
  The order notes give the finished thickness but not the laminate, so εr is assumed. Copper loss is
  included.
- **Ports:** the chip pin (H23 to the F23 ground ball), each ladder part, R70 and the C20 pads. Twelve
  decoupling capacitors near the chip are modelled as shorts. `fieldsolve cosim` places measured parts
  on the ports: Murata GRM1555 for the Fenghua capacitors, which have no published data, and the exact
  LQG15HS inductors.
- **Not modelled:** solder mask, enclosure, battery. All three pull the resonance down.
- **A solver bug was found and fixed along the way.** Every part here is rotated 45°, the chip included.
  Ports between pads that don't share a mesh line used to get a straight bridge that missed the far pad
  and shorted the port: the first run showed the chip seeing 0.05 + j26 Ω. Such ports now get a
  staircase bridge, with the gap in the copper-free stretch between the pads (commit 2918126, tests in
  `review_tests.rs`). No earlier board had rotated parts.

## The antenna alone

Seen at the feed (the C20 pads), with R70 open:

| | 2.40 GHz | 2.44 GHz | 2.48 GHz |
|---|---|---|---|
| impedance | 37.5 + j5 Ω | 28.5 + j16 Ω | 21.9 + j26 Ω |
| S11 (50 Ω) | −16.3 dB | −9.5 dB | −6.0 dB |
| radiation efficiency | 88 % | 85 % | 82 % |

## As built, seen from the chip pin

The score is the power the chip delivers into its design load (1 − |Γ|² against Nordic's measured network,
35 − j35 Ω at 2.44 GHz), times the share of the power entering the pin that is radiated (ladder plus
antenna loss).

| | 2.40 GHz | 2.44 GHz | 2.48 GHz |
|---|---|---|---|
| chip sees | 68.5 − j19 Ω | 94 − j21 Ω | 111 − j67 Ω |
| \|Γ\| against the design load | 0.31 | 0.46 | 0.58 |
| radiated / accepted | 70 % (−1.5 dB) | 62 % (−2.1 dB) | 52 % (−2.9 dB) |
| score | −2.0 dB | −3.1 dB | −4.6 dB |

Treating R70 as a real 0 Ω jumper (about 0.4 nH) instead of an ideal short costs another 0.3-0.5 dB.

## Radiation pattern

At 2.44 GHz the peak realized gain is +1.0 dBi, for a 50 Ω source at the chip pin. The 60 mm board is
about half a wavelength long, so its ground currents radiate like a dipole along the board's long axis:

- a near-omnidirectional ring of about −1 dBi (E_φ) in the plane across the board (x-z);
- a figure-of-eight of E_θ in the y-z plane, with nulls off the board's ends;
- the antenna arm's own x-directed current fills those nulls only partly.

![Polar cuts at 2.44 GHz: board plane, x-z and y-z planes, total and both polarisations](img/mth7_polar_2440.jpg)

![Polar cuts at 2.40, 2.44 and 2.48 GHz](img/mth7_polar_band.jpg)

![Far field at 2.44 GHz over the board (radius and colour both EIRP per W into the chip pin)](img/mth7_farfield_2440.jpg)

![Far field at 2.44 GHz and at H2 and H3, on a shared scale per row](img/mth7_farfield_h1h2h3.jpg)

Radiated harmonics, per watt into the chip pin (50 Ω): **H2 −25.5 dB, H3 −21.6 dB.**

## Retune

Every ladder position was varied over real part values:

- capacitors: Murata GRM1555, 0.2-8.2 pF, or open;
- inductors: Murata LQG15HS, 1.0-8.2 nH, or short;
- R70: a short, an inductor up to 4.7 nH, or a series capacitor;
- C20: open or up to 2 pF.

Each candidate was scored with `fieldsolve cosim` on the existing 7-port EM result, by the worst score
of 2.40 / 2.44 / 2.48 GHz. The search was coordinate descent from 60 random starts plus every one- and
two-part change from as built: about 25,000 candidates, taking seconds.

| | 2.40 GHz | 2.44 GHz | 2.48 GHz | radiated H2 / H3 |
|---|---|---|---|---|
| as built | −2.0 dB | −3.1 dB | −4.6 dB | −25.5 / −21.6 dB |
| **L5 → 0 Ω, C20 = 0.5 pF** | −1.6 | −1.4 | −1.6 | −25.9 / −21.6 |
| same, with real 0 Ω jumpers (0.4 nH each) | −1.3 | −1.4 | −2.1 | −25.4 / −21.6 |
| L5 → 0 Ω only | −1.3 | −1.4 | −1.8 | −24.6 / −21.6 |
| best with harmonics no worse: 1.0p / 3.9n / 1.0p / 1.6n, R70 → 2.4 pF series, C20 1.1 pF | −1.6 | −1.4 | −1.6 | −25.1 / −21.1 |
| best overall: 0.4p / short / 0.4p / short, R70 → 1.6 nH, C20 0.7 pF | −1.2 | −1.1 | −1.2 | −22.1 / −18.4 |

- **Why it works:** at 2.44 GHz the PIFA's feed is already inductive, and L5 added more inductance in
  series. Removing it fixes most of the match, and C20 trims it across the band.
- **It's robust:** C20 anywhere from 0.4 to 0.6 pF lands within 0.1-0.2 dB.
- **What's left:** after the fix, about 1.2-1.4 dB is lost between the chip pin and radiation. About
  0.7 dB of that is the antenna and board.
- **Further retuning isn't worth it:** retunes beyond about −1.6 dB worst case trade H2 for ~0.4 dB.

## Caveats

- **Free-space values only.** No enclosure, solder mask or battery; the antenna will also be trimmed.
- **The Fenghua capacitors are modelled with Murata data.**
- **The chip's optimum load** is taken as Nordic's measured target. How far off it is in band is
  examined in [1bp-rf-findings.md](1bp-rf-findings.md).
- **The harmonic levels** are per watt into a 50 Ω port. The chip's real harmonic source impedance is
  unknown, and the model's levels near the nulls are mesh-sensitive by a few dB.

## Reproducing

- **EM:** a `fieldsolve run` spec with `"gerbers"` the KiCad conversion's `fabrication/gerbers`,
  `"coords": "kicad"`, `"extent": "board"`, ROI `[145, 70, 159.6, 87.5]`, cell 0.05, air 15 mm,
  `"dielectrics": [{"t": 0.73, "er": 4.5, "tand": 0.015}]`, and wideband loss with copper. Ports
  `CHIP` U12.H23→U12.F23, `C59` .2→.1, `L4` .1→.2, `C60` .2→.1, `L5` .1→.2, `R70` .2→.1, `C20` .1→.2.
  Decoupling shorts on C5-C8, C12, C15, C16, C23, C24, C29, C30, C40. Far field at 2.40/2.44/2.48/4.88/
  7.32 GHz.
- **Cosim:** `keep` `["CHIP"]` with `"target": "nordic_ref_network.s2p"` and `"farfield": true`. For the
  antenna alone, keep `["C20"]` with R70 and CHIP open.
- **Retune:** `tools/ladder_retune.py` (the MTH version: coordinate descent over part libraries,
  scored as above).
- **Plots:** `tools/polar_farfield.py` and `tools/render_farfield_panels.py`.
