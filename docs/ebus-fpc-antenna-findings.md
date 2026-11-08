# eBUS_IF inline: FPC antenna in the lid pocket (fieldsolve, Nov 2026)

The eBUS_IF inline adapter (`p-hw` branch `claude/ebus-if-inline`, `Revision_7/eBUS_IF_inline`, draft A.4 at
81d5260c) has no PCB antenna. The nRF52840's configuration-7 front end ends in a u.FL (J3). The moulded
enclosure (`enclosure-moulded/`) puts an adhesive FPC antenna in a 0.15 mm pocket in the cover ceiling,
centred at KiCad x 20.0, y 16.4. The README specifies it as "34.9 x 9.0 x 0.1 mm, Molex 146153 class". Its
cable plugs into J3.

## Summary

- **With the FORCE YL1439-BT-V1, the FPC used on other boards, it is worse:** 26-45 % radiation efficiency
  in the pocket, and −7.5 to −11.6 dB in total from 50 Ω, depending on the cable. In free space on 1.05 mm of
  PC/ABS it is 85-89 %. The part appears tuned for about 0.7 mm of plastic. See
  [the YL1439 section](#the-force-yl1439-bt-v1-the-fpc-used-on-other-boards).
- The rest of this summary is about the stand-in dipole first modelled for the "Molex 146153 class" part.
- **The pocket puts the FPC about 8.3 mm above the board's ground, parallel to it, and that wrecks it.**
  - The board is 50 × 28 mm, a 4-layer 0.8 mm board with ground on the inner layers, and it fills the
    enclosure. The FPC's copper sits at z 9.0, and the board's top copper at z 0.67.
  - In free space on a 1.05 mm PC/ABS plate, the stand-in FPC used here is 35 + j4 Ω at 2.44 GHz and
    radiates 97 % of what it accepts.
  - In the pocket, it is 9 + j46 Ω (S11 −1.8 dB) and radiates 55-62 % across 2.40-2.50 GHz.
  - Fed from 50 Ω, the total in band is −7.1 to −7.6 dB, against −0.3 dB in free space.
- **The board causes it; the plastic hardly does.**
  - With the enclosure but no board, the FPC detunes by about 60-80 MHz and still radiates 96-97 %.
  - With the board but no enclosure, its resistance falls to 5-10 Ω and its efficiency to 69-79 %.
  - It behaves like a low patch: its own image in the ground plane cancels most of its radiation
    resistance. The currents that remain are large, so the loss in the board's FR4 (tan δ 0.02) and
    copper, and in the plastic, takes a big share.
  - Its directivity rises from 1.9 dBi to 6.5 dBi, beamed up through the lid.
- **Rematching can't fix it.** A matching network retuned for 9 + j46 Ω would remove most of the mismatch.
  The radiation efficiency would still cap the total at about −2.1 to −2.6 dB, over a narrow band
  (R ≈ 9 Ω, so the matched bandwidth is small), and the mismatch through the 1.13 mm coax would add cable
  loss on top.
- **Placement is the lever.** An FPC needs clearance from ground planes; vendors typically ask for 10-20 mm
  or more, or for mounting perpendicular to the board. In this enclosure the board fills the plan, so a
  ceiling pocket has no ground-free area to use. Options to simulate next:
  - an end or side wall, perpendicular to the board;
  - a taller cover;
  - a PCB antenna in a ground-free strip at the board's end.

## Model

![Model: a section at x 15.3 through the feed (STL enclosure outlines, board, FPC and carrier, cable-shield riser, screws), and the plan (FPC copper, cable shield, screws, port)](img/ebus_fpc_model.jpg)

- **Board:** the A.4 Gerbers (`exports/fab/gerbers`), with the 4-layer stack-up from the job file. It uses
  `"preset": "antenna"` (surface-impedance copper, wideband FR4 loss, sub-cell edges, mesh refinement,
  the whole board).
  - Decoupling: every capacitor of 47 pF or more is a short (22 parts). The 12 pF crystal loads and the
    1-1.2 pF match are left open.
  - Chips (U4 nRF52840, U2, U5, U6, U1) have 3 pF on each signal pin.
  - The eBUS lines (EBUS_A, EBUS_B, TX_B) and CC_S are acknowledged as open. In use they run on into the bus
    cable, which isn't modelled.
- **Enclosure:** `build/asm-base.stl` and `asm-cover.stl` as dielectric bodies. The material is FR PC/ABS,
  with εr 2.9 and tan δ 0.008 assumed at 2.4 GHz (datasheets give 1 MHz values).
- **Screws:** the three ST2.2 screws (H1, H2, H3) are metal posts of radius 1.1 mm from z −2.8 to 8.3.
- **FPC:** a stand-in. Molex doesn't publish the copper pattern, so this is an off-centre-fed dipole of the
  same outline (34.9 × 9.0 mm).
  - A 13 × 8 mm shield-side patch at the left.
  - A 1.5 mm radiator arm along the edge nearer J3, folded up the right end and back 8.5 mm along the top.
    That return length tunes it to about 2.44 GHz on a 1.05 mm PC/ABS plate in free space, the vendor-style
    condition.
  - Copper at z 9.00 on a 0.1 mm polyimide carrier (εr 3.4, z 8.95-9.05) in the pocket.
  - Port: 50 Ω across the 0.5 mm feed gap.
- **Cable:** the shield of the 1.13 mm coax is modelled by its outer surface only. It rises from the FPC's
  shield pad to z 3.0, runs straight to the mated u.FL plug on J3, and drops to J3's ground pads. Without
  it, the installed results change by under 0.3 dB (see the table).
- **Mesh:** 0.1 mm over the FPC, and at most 0.4 mm over the enclosure. A finer-mesh check at 0.07 mm is
  below.
- **Solver features:** these runs use fieldsolve's new `"bodies"` (STL and box dielectrics, metal posts) and
  `"layers"` (copper off the board) keys; see [fieldsolve.md](fieldsolve.md#usage).

## Results

| Case | Z at 2.44 GHz | S11 at 2.44 | radiation efficiency 2.40 / 2.44 / 2.48 / 2.50 GHz | total from 50 Ω at 2.44 | peak directivity |
|---|---|---|---|---|---|
| free space, 1.05 mm PC/ABS plate | 35 + j4 Ω | −14.8 dB | 96.6 / 96.9 / 97.2 / 97.3 % | −0.28 dB | 1.9 dBi |
| enclosure, no board | 41 + j23 Ω | −11.8 dB | 96.3 / 96.6 / 96.8 / 96.8 % | −0.45 dB | 2.0 dBi |
| board, no enclosure | 5 − j21 Ω | −1.6 dB | 68.9 / 76.8 / 79.1 / 78.1 % | −6.3 dB | 6.6 dBi |
| **installed** (enclosure, board, screws, cable shield) | **9 + j46 Ω** | **−1.8 dB** | **55.5 / 55.5 / 58.8 / 61.8 %** | **−7.35 dB** | **6.5 dBi** |
| installed, no cable shield | 9 + j38 Ω | −2.0 dB | 52.3 / 52.7 / 56.1 / 59.3 % | −7.05 dB | 6.4 dBi |

![|S11| (left) and radiation efficiency with the total from a 50 Ω feed (right) for the four cases](img/ebus_fpc_cases.jpg)

**Two-mesh check.** The installed case at 0.07 mm cells (18.3 M cells) against 0.1 mm (10.0 M):

| | 2.40 | 2.44 | 2.48 | 2.50 GHz |
|---|---|---|---|---|
| Z, 0.1 mm | 8 + j35 | 9 + j46 | 10 + j57 | 11 + j63 Ω |
| Z, 0.07 mm | 8 + j31 | 9 + j41 | 10 + j52 | 10 + j58 Ω |
| radiation efficiency, 0.1 / 0.07 mm | 55.5 / 52.7 | 55.5 / 56.3 | 58.8 / 61.2 | 61.8 / 64.3 % |
| total from 50 Ω, 0.1 / 0.07 mm | −7.12 / −7.12 | −7.35 / −7.20 | −7.55 / −7.32 | −7.61 / −7.38 dB |

The conclusion holds on both meshes. The efficiency is mesh-sensitive by about ±3 points, and the
reactance by about 5 Ω.

## The FORCE YL1439-BT-V1 (the FPC used on other boards)

ShenZhen Powerful Technology's FPCs are marked "FORCE". Two are on file:
- the YL1439 2.4/5 GHz Wi-Fi V2, 46 × 10 mm and dual band, specified for Modbus_IF in
  `Revision_6/Modbus_IF/docs/manufacturing/original/KNA-019286R0_YL1439_antenna_spec.pdf`;
- the single-band **YL1439-BT-V1**, 36 × 10 mm, whose vendor AutoCAD drawing (Oct 2022) shows the copper
  outline. That drawing is `YL1439-BT天线-FPC.pdf` in Downloads.

The BT part is the logical one for a 2.4 GHz-only design. It is 0.5 × 0.4 mm larger than the pocket drawn for
the Molex part, so the pocket would grow to about 36.5 × 10.4 mm.

![YL1439-BT-V1: copper as traced (orange) over the vendor drawing; feed port in red](img/ebus_yl1439_trace.jpg)

- **The pattern is a printed inverted-F.**
  - A 1.25 mm arm runs about 28 mm from an open end at the left (with a tab) to a feed pad.
  - A 0.5 mm bridge over a diagonal slot shorts the arm to the block where the cable lands.
  - A 1.4 mm ground strip runs back from that block.
  - The 1.13 mm cable (63 mm) leaves along the FPC's axis. The model's port spans the slot from the feed pad to
    the block.
- **No measurement of this part is on file**, so the model isn't validated. The Wi-Fi V2's report (62-81 %
  efficiency, 5.0-5.8 dBi peak gain) is a different pattern.
- **It is tuned for a thin plastic carrier.**
  - In free air (no cable) it sits at 2.85 GHz.
  - On a PC/ABS plate, it centres on the band with **about 0.7 mm** of plastic behind it: best match at
    2.46 GHz, VSWR ≤ 2 from 2.42 to 2.50 GHz.
  - Under the eBus lid's 1.05 mm (the 1.2 mm wall less the 0.15 mm pocket) it sits at 2.40 GHz. Each further
    0.5 mm lowers it about 40 MHz.
  - On the 1.05 mm plate it radiates 85-89 %, with a peak gain of about 1.5 dBi.
  - These free-space runs agree at 0.1 and 0.07 mm cells: the same best-match frequency (2.397 GHz), with
    the impedance within 1-3 Ω.
- **Installed, it radiates only 26-33 %.** These are converged numbers: 0.07 and 0.05 mm agree.

  | 0.05 mm mesh | 2.40 | 2.44 | 2.48 | 2.50 GHz |
  |---|---|---|---|---|
  | Z | 7 + j59 | 15 + j79 | 72 + j119 | 169 + j85 Ω |
  | S11 | −1.0 | −1.5 | −3.0 | −4.1 dB |
  | radiation efficiency | 32.9 | 31.6 | 28.3 | 27.3 % |
  | total from 50 Ω | −11.6 | −10.4 | −8.5 | −7.8 dB |

  - The antenna's own resonance moves up to about 2.5-2.55 GHz.
  - A sharp mode at 1.93-1.96 GHz (Q about 32) is still ringing when the runs stop, and the 0.05 mm run hit
    its 200 000-step limit.
  - The 0.1 mm installed run put the resonance near 2.15 GHz instead. **It isn't converged.** The 0.1 mm
    shortened-arm runs from the same mesh are void too.
  - The stand-in dipole above did agree between 0.1 and 0.07 mm.

- **The 1.93-1.96 GHz mode belongs to the cable.**
  - Without the cable shield (0.07 mm), it is gone. The antenna then resonates near 2.3 GHz and sees
    10 − j24 / 5 − j7 / 4 + j1 / 3 + j4 Ω at 2.40 / 2.44 / 2.48 / 2.50 GHz, radiating 36 / 40 / 40 / 45 %
    (total −7.5 to −10.1 dB).
  - That run agrees with the 0.1 mm one (9 − j18 / 5 − j5 Ω at 2.40 / 2.44 GHz), so the cable model is what
    needed the finer mesh: its strip runs 0.4 mm under the FPC.
  - The real cable's route sits between those two cases. Either way the YL1439 does worse here than the
    stand-in dipole (55-62 %). The pocket over the board's ground is still the problem; a thinner lid would
    only retune it.

![YL1439-BT-V1: |S11| against the plastic thickness in free space (left), and installed on three meshes against free space (right)](img/ebus_yl1439_cases.jpg)

## Caveats

- **The FPC is a stand-in.** A real Molex 146153 has its own pattern and tuning, and its measured
  impedance and efficiency on a ground plane would differ. The mechanism doesn't depend on the pattern,
  though: any small horizontal antenna 8 mm over a ground plane its own size loses most of its radiation
  resistance. The real part's datasheet clearance rules apply.
- **The cable is crude.** Only the shield's outer surface is modelled, along one route. A real cable's
  route, length and any ferrite change the common-mode currents.
- **Not modelled:** the eBUS cable on J1, the plug's metal contacts, and the parts' bodies on the board
  (the inductors L1-L3, the electrolytics C7/C8, SW1). Only their pads are there.
- **The enclosure material is estimated:** εr 2.9, tan δ 0.008.
- **The matching network isn't included.** These are antenna-alone results at the FPC's feed.

## Reproducing

`tools/ebus_fpc_specs.py` writes the specs (run it from a directory holding the board's Gerbers as `gbr/`
and the air "dummy" board as `dummy_gbr/`; see its docstring):

- `python3 ebus_fpc_specs.py free fpc_free 8.5 ff`
- `python3 ebus_fpc_specs.py inst fpc_inst 8.5 ff` (`nocable` for the run without the shield)
- The enclosure-only and board-only cases drop the board or the `base`/`cover` bodies from those specs.
