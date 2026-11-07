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

> **Correction (7 Nov 2026).** The EM run behind the first version of this note meshed the antenna's open
> end coarsely. Its fine-mesh region started at x 145 mm, but the radiator arm runs from x 141.9 to 144.35,
> so the open end fell in graded cells of 0.35-0.5 mm. A coarse open end acts electrically longer, so the
> resonance came out too low.
>
> Re-run with the fine region from x 141.5 (one port at C20, in the fast-mode work), the antenna resonates
> near 2.43 GHz, not below 2.40 GHz. The old ROI reproduced the 7-port run to within about 1 Ω, so the ROI
> alone explains the difference. The antenna and ladder results below have been corrected; superseded
> values are marked.
>
> The corrected values are still provisional:
>
> - **No new 7-port EM run yet.** The ladder numbers come from the old 7-port result with the
>   better-resolved antenna swapped in at the feed. See [Retune](#retune) for the method and its error.
> - **The mesh is not converged at the open end.** Going from about 0.4 mm cells to 0.1 mm to 0.05 mm
>   moved the antenna's worst-case score from −2.09 to −1.32 to −0.91 dB, and the trend is still moving.
>
> A full 7-port re-run, with the open end resolved and a convergence check, is planned once the mesh
> refinement work lands.

## Summary

- **The antenna is good.** In free space it resonates near 2.43 GHz, mid-band, and radiates about 87-91 %
  of the power it accepts. (Superseded: "just below 2.40 GHz, 82-88 %".) The enclosure, solder mask and
  battery will all pull the resonance down, so there is less length in hand for trimming than first
  thought.
- **The ladder isn't matched to it, but by less than first reported.** The configuration-7 values were
  designed for a 50 Ω connector. On the PIFA (about 44 + j2 Ω at 2.44 GHz), the chip sees about 60 − j36 Ω
  instead of its 35 − j35 Ω design load. The worst in-band case is about −2.5 to −3.1 dB, at 2.48 GHz.
  (Superseded: 94 − j21 Ω, −4.6 dB.)
- **The earlier fix, L5 → 0 Ω and C20 = 0.5 pF, is no longer the right one.** On the corrected antenna it
  still helps at 2.48 GHz, but it pulls the chip impedance down to about 22-25 Ω. That makes it worse than
  as built at 2.40 and 2.44 GHz, and its worst case is only −1.8 to −2.2 dB.
- **Provisional cheapest fix: C59 and C60 → 0.6 pF, with L4, L5, R70 and C20 unchanged.** This is a
  two-part value change with no new footprint. The chip sees about 61 − j37 / 43 − j36 / 30 − j28 Ω across
  the band, close to its design load. The worst case improves from about −2.5 to −3.1 dB to about −1.0 to
  −1.4 dB, and it beats the old fix at every frequency. Check it against the 7-port re-run before changing
  the BOM.
- **H2 suppression is weak.** The radiated H2 is −25.5 dB per watt into the chip pin (from the coarse
  open-end run, not yet re-evaluated), against −29.7 dB for
  1BP as built and −49.3 dB for Nordic's reference layout. C59's ground pad is flooded into the top pour,
  the same deviation found on 1BP ([1bp-rf-findings.md](1bp-rf-findings.md)). The same fix should apply:
  an isolated return to VSS_PA. That fix hasn't been simulated here.
- **Retune again after trimming.** The trimmed antenna in the enclosure will present a different feed
  impedance. Re-run the EM model on the trimmed copper, plus the enclosure if there is data, and repeat
  the retune. The retune takes seconds; the EM run about 45 minutes.

![RF section, top copper (left) and bottom (right): the diagonal ladder from the chip to R70, the PIFA feed and short above](img/mth7_rf_copper.jpg)

## Model

- **The whole board** in free space, from the KiCad conversion's X2 Gerbers. 0.05 mm cells over the ROI
  [145, 70, 159.6, 87.5], graded outside. That ROI misses the radiator's open end, x 141.9-144.35 (see
  the correction at the top), with far field from a Huygens box in 15 mm of air.
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
| impedance | 59.8 − j8.8 Ω | 44.3 + j2.3 Ω | 33.4 + j13.2 Ω |
| S11 (50 Ω) | −18.5 dB | −23.8 dB | −12.0 dB |
| radiation efficiency | 91 % | 89 % | 87 % |
| *superseded (open end coarse):* impedance | *37.5 + j5 Ω* | *28.5 + j16 Ω* | *21.9 + j26 Ω* |
| *superseded:* S11 (50 Ω) | *−16.3 dB* | *−9.5 dB* | *−6.0 dB* |
| *superseded:* radiation efficiency | *88 %* | *85 %* | *82 %* |

The corrected values are from a one-port run at C20 (R70, CHIP and the ladder open), with the fine region
from x 141.5. It is otherwise identical to the 7-port run.

## As built, seen from the chip pin

The score is the power the chip delivers into its design load (1 − |Γ|² against Nordic's measured network,
35 − j35 Ω at 2.44 GHz), times the share of the power entering the pin that is radiated (ladder plus
antenna loss).

| | 2.40 GHz | 2.44 GHz | 2.48 GHz |
|---|---|---|---|
| chip sees | 55 − j45 Ω | 60 − j36 Ω | 75 − j41 Ω |
| \|Γ\| against the design load | 0.21 | 0.26 | 0.40 |
| score | −1.1 dB | −1.5 to −1.7 dB | −2.5 to −3.1 dB |
| *superseded:* chip sees | *68.5 − j19 Ω* | *94 − j21 Ω* | *111 − j67 Ω* |
| *superseded:* score | *−2.0 dB* | *−3.1 dB* | *−4.6 dB* |

The ranges are the two antenna-swap estimates described under [Retune](#retune).

Treating R70 as a real 0 Ω jumper (about 0.4 nH) instead of an ideal short cost another 0.3-0.5 dB in the
superseded model. This hasn't been re-evaluated.

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

The pattern, gain and harmonic figures in this section are from the run with the coarse open end. The
dipole-like board mode should not change much, but the numbers will move a little.

## Retune

### Corrected (provisional)

There is no 7-port EM run with the open end resolved yet. Instead, the old 7-port result (whose ladder
region is meshed finely) was re-scored with the better-resolved antenna swapped in at the feed. The swap
uses the change in antenna impedance between the old-ROI and widened-ROI one-port runs, and scales the
antenna's radiated/accepted share the same way.

There are two ways to make the swap, and neither is exact, because the antenna and ladder copper couple:

- **series swap at the R70 port:** checked against the old far-field scores, it reads 0.3-0.5 dB
  optimistic;
- **shunt-admittance swap at the C20 node:** it reads 0-0.6 dB pessimistic.

The two agree within 0.1 dB at 2.40 and 2.44 GHz and differ by up to 0.6 dB at 2.48 GHz. Both give the chip
impedance exactly on the old model. The table gives series / shunt.

| | chip sees, 2.40 / 2.44 / 2.48 GHz | 2.40 GHz | 2.44 GHz | 2.48 GHz |
|---|---|---|---|---|
| as built | 55 − j45 / 60 − j36 / 75 − j41 Ω | −1.1 / −1.1 dB | −1.5 / −1.7 | −2.5 / −3.1 |
| **C59, C60 → 0.6 pF** | 61 − j37 / 43 − j36 / 30 − j28 Ω | −0.7 / −0.6 | −0.8 / −0.8 | −1.0 / −1.4 |
| C59, C60 → 0.5 pF | 71 − j46 / 50 − j50 / 33 − j43 Ω | −0.9 / −0.8 | −1.0 / −1.0 | −0.9 / −1.3 |
| C59 0.6 pF, C60 0.5 pF | 55 − j44 / 38 − j44 / 26 − j37 Ω | −0.6 / −0.6 | −0.8 / −0.8 | −1.0 / −1.4 |
| 0.5p / 2.7n / 0.5p / 2.0n (four parts) | 56 − j24 / 45 − j37 / 29 − j38 Ω | −0.6 / −0.4 | −0.7 / −0.7 | −0.9 / −1.2 |
| R70 → 1.8 pF series | 35 − j52 / 32 − j40 / 36 − j29 Ω | −1.3 / −1.2 | −1.3 / −1.2 | −1.5 / −1.8 |
| L5 → 0 Ω | 34 − j59 / 29 − j48 / 30 − j37 Ω | −1.5 / −1.5 | −1.3 / −1.3 | −1.2 / −1.5 |
| old fix: L5 → 0 Ω, C20 = 0.5 pF | 23 − j60 / 22 − j50 / 25 − j42 Ω | −2.2 / −1.8 | −1.8 / −1.5 | −1.5 / −1.5 |

- **Why the fix changed:** the resolved antenna is close to 50 Ω and nearly real in band, not inductive.
  The configuration-7 ladder transforms it too high, about 60 Ω real. Lighter shunt capacitors bring the
  real part down to the design load's 35 Ω. Removing L5 now overshoots.
- **Search:** a coordinate descent over the same part libraries (series swap) found worst cases of about
  −0.8 to −0.9 dB, but only by changing three or four parts. Given the swap's error at 2.48 GHz, that is
  not clearly better than the two-capacitor change.
- **Harmonics:** not re-evaluated. The H2/H3 effect of lighter shunt capacitors needs the 7-port re-run
  with far field.

### Superseded: retune on the coarse open end

The rest of this section is the original retune. It was fitted to the coarse antenna (resonance below
2.40 GHz) and is kept for reference only.

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
- **The open end's mesh isn't converged** (see the correction at the top), and the corrected ladder
  numbers use an antenna swap rather than a 7-port run.
- **The Fenghua capacitors are modelled with Murata data.**
- **The chip's optimum load** is taken as Nordic's measured target. How far off it is in band is
  examined in [1bp-rf-findings.md](1bp-rf-findings.md).
- **The harmonic levels** are per watt into a 50 Ω port. The chip's real harmonic source impedance is
  unknown, and the model's levels near the nulls are mesh-sensitive by a few dB.

## Reproducing

- **EM:** a `fieldsolve run` spec with `"gerbers"` the KiCad conversion's `fabrication/gerbers`,
  `"coords": "kicad"`, `"extent": "board"`, ROI `[141.5, 70, 159.6, 87.5]` (the first run used x0 = 145,
  which leaves the antenna's open end coarse), cell 0.05, air 15 mm,
  `"dielectrics": [{"t": 0.73, "er": 4.5, "tand": 0.015}]`, and wideband loss with copper. Ports
  `CHIP` U12.H23→U12.F23, `C59` .2→.1, `L4` .1→.2, `C60` .2→.1, `L5` .1→.2, `R70` .2→.1, `C20` .1→.2.
  Decoupling shorts on C5-C8, C12, C15, C16, C23, C24, C29, C30, C40. Far field at 2.40/2.44/2.48/4.88/
  7.32 GHz.
- **Cosim:** `keep` `["CHIP"]` with `"target": "nordic_ref_network.s2p"` and `"farfield": true`. For the
  antenna alone, keep `["C20"]` with R70 and CHIP open.
- **Retune:** `tools/ladder_retune.py` (the MTH version: coordinate descent over part libraries,
  scored as above).
- **Plots:** `tools/polar_farfield.py` and `tools/render_farfield_panels.py`.
