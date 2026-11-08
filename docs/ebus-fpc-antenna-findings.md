# eBUS_IF inline: FPC antenna in the lid pocket (fieldsolve, Nov 2026)

The eBUS_IF inline adapter (`p-hw` branch `claude/ebus-if-inline`, `Revision_7/eBUS_IF_inline`, draft A.4 at
81d5260c) has no PCB antenna. The nRF52840's configuration-7 front end ends in a u.FL (J3). The moulded
enclosure (`enclosure-moulded/`) puts an adhesive FPC antenna in a 0.15 mm pocket in the cover ceiling,
centred at KiCad x 20.0, y 16.4. The README specifies it as "34.9 x 9.0 x 0.1 mm, Molex 146153 class". Its
cable plugs into J3.

## Summary

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
