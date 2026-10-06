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

## Round 2: escape planning, pad entry, stack-up (6 October 2026)

AMOLED planes profile, 90 s, four seeds each (net-order jitter), in-process deficit (matches native).
Seed noise is large (base 31-50), so single runs are not evidence.

| Configuration | Deficits (seeds 0-3) | Mean |
|---|---|---:|
| Base (fanout off, late fanout, escape reservations) | 31, 40, 50, 42 | 40.8 |
| + local buses only | 45, 42, 45, 50 | 45.5 |
| + escape plan for U10/J3 | 40, 48, 50, 42 | 45.0 |
| + escape plan, all dense packages | 47, 54, 49, 52 | 50.5 |
| Base, outer-layer cost x2 / x3 | means | 39.0 / 41.0 |
| **Base + plane ties (new default)** | 31, 30, 41, 31 | **33.3** |
| + escape plan U10/J3 | 35, 38, 43, 38 | 38.5 |
| + axial pad entry (0.1 mm halo) | 37, 42, 37, 41 | 39.3 |
| + axial pad entry (one-cell halo) | 37, 39, 38, 35 | 37.3 |
| 8 layers (In5/In6 signal), plane ties | 16, 25, 38, 26 | 26.3 |

Native KiCad, best candidates, all with zero new errors: 6 layers 30, 8 layers 16, 6 layers with
axial entry 35.

Findings:

- Escape planning works mechanically (dog-bones with staggered and between-row vias, crossing-free
  buses), but completion does not improve: U10/J3 pins fail later, when their long runs to U1 are
  sealed off, not at the package. Kept as `escape_plan`, off.
- The router already puts a lot of length on inner layers (base 434 mm inner vs reference 701 mm);
  segment counts on F/B were misleading. A higher outer-layer cost moves length inward without
  helping completion.
- The biggest win was plane pads without room for a via: they are now tied to a neighbouring pad's
  drop (`_tie_to_dropped`), as a designer would.
- Two more signal layers help (mean 33 -> 26, best 30 -> 16) but do not finish the board: the
  remaining failures are pin access at U10, U3, J3 and J1 plus the I2C nets.
- Axial pad entry costs about 4 connections on this board. The continuous checker flags the RF net at
  U1 (class clearance), which the project's courtyard rule permits natively.
