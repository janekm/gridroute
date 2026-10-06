# Routing results and plan (6 October 2026)

Uniform policy: `suite.py --deadline 30 --options-json '{"continuous_check":true}'`, manifest options
per case, Apple M3 Max / Metal, serial. Baseline: the handover's full frozen-source replay
(`final-suite-v3`, same settings). Acceptance is native KiCad: zero unconnected items, no new errors,
unchanged footprints/rules.

| | Handover replay | This branch |
|---|---:|---:|
| Profiles passing | 24 / 39 | 24 / 39 |
| Native unconnected items (all profiles) | 337 | 259 |
| New DRC errors (all profiles) | 7 | 0 |

Every profile that changed improved; nothing regressed (passing boards keep their times).

| Profile | Handover | Now |
|---|---:|---:|
| genz_amoled18_planes | 59 | **31** |
| genz_amoled18 (cold) | 62 | 56 |
| TMC261-stepstick | 5 | **1** |
| miniboard-stm32f0 | 9 + 6 errors | **3 + 0** |
| ESP32-board | 6 + 1 error | 5 + 0 |
| Hardware_Playground rfm69w | 17 | 12 |
| Hornbill (rule barrier) | 31 | 26 |
| OtterPill / C7 shifts | 35 / 35 / 35 | 27 / 28 / 27 |

AMOLED planes, 31 remaining: U10/J3 (display level shifter and FPC, 0.4 mm pitch) 9, U3 (PMIC QFN-40) 3,
the rest scattered. 256 vias, 0 continuous findings, 23 dangling-track warnings from failed nets' stubs.
The cold profile has to route GND (105 pads) and VCC3V3 (42 pads) as tracks; the source design uses
In1/In4 planes, which the `_planes` profile reproduces.

## What changed

- **Plane nets**: via drops instead of routing; plane-aware scoring and continuous check.
- **Neck-down** (DRU courtyard width exceptions): thin search regions and straight escape stubs to where
  the class width fits; such nets route first. The source design itself breaks its 0.35 mm PWR minimum
  (199 width errors); we keep the rule.
- **Escape reservations** for fine-pitch pins during the first pass (TraceMaker idea).
- **Component joining**: a boxed-in first pad no longer fails a whole multi-pad net (cold VCC3V3 went
  from 41 missing pads to one island).
- **Clearance fixer**: continuous-check findings are rerouted with a wider margin or dropped.
- **Kernel**: optional per-state cost in A* (`gr_astar_hybrid_cost`); negotiated rip-up mode built on it.
- **Tools**: `portfolio.py` (parallel variants, best by in-process check), `render.py`.

## TraceMaker comparison (github.com/DingoOz/TraceMaker)

C++20/CUDA, octilinear lattice A* with exact-geometry legality, PathFinder history, an 8-variant
portfolio and lossless KiCad I/O. They report PCBench tiers A/B/C/D at 100 / 67.5 / 60 / 50 % KiCad-clean
(120 s, 8 variants). Overlapping boards are scored differently (connections routed, not native
unconnected), so not directly comparable: USBee32-S2 159/159 (ours: strict 0.4 mm profile 6 open;
native-minimum profile passes in 4.4 s), OtterPillG 103-104/109, 1Bitsy 149/160 on human placement
(ours passes in 2.6 s).

Adopted: escape reservations, history-cost search (kernel), portfolio, via-cost polish
(`polish_via_cost_mm`, via `Board.relax`), and their raster-margin lesson (0.71 x pitch). Not yet:
exact-geometry legality per step, reachability proofs to skip dead connections, hardest-first restarts,
via LNS.

## Plan

1. **Fine-pitch escape planning** (U10/J3/U3 = 12 of AMOLED's 31): assign each pin of a dense package
   an exit (same-layer outward, or a via site between rows / staggered outward) as one matching
   problem, then reserve those sites, instead of first-come via placement.
2. **Exact legality at the boundary**: check emitted segments against continuous geometry while
   searching near obstacles, so the 0.5 x pitch margin can stay (0.75 costs ~14 connections on AMOLED).
3. **Negotiation that wins**: history costs plus hardest-first restarts keeping the best snapshot;
   the current negotiate mode rips too many nets per move (cost tuning, per-victim rip caps).
4. **Human-like quality**: enable the via polish by default once measured on passing boards; trim
   failed nets' stubs; tree-style power (trunk then branches) for PWR nets; loop-area metric for
   decoupling paths.
5. **Cold AMOLED**: emulate the source planes with dedicated GND/VCC3V3 track layers, or accept that
   this profile is a planes problem.
6. **Determinism under load**: expansion budgets instead of wall-clock deadlines inside portfolios.
