# fieldsolve fast iteration: component tuning and fidelity policies

November 2026. This implements the first two stages of the
[fast-mode proposal](../reports/pcb-field-solver-review-2097498/fast-mode/proposal.md) (outside the repository,
in `reports/pcb-field-solver-review-2097498/fast-mode/`):

1. **Component changes reuse one EM extraction.** `fieldsolve tune` scores part combinations on a saved
   N-port result, with no new EM run.
2. **Geometry changes get cheaper EM runs with a measured accuracy budget.** `fieldsolve run SPEC --fast`
   resolves to a recorded fidelity policy (`screening`, `refined` or `validated`).

Stages 3 and 4 (surrogate-assisted outer search, adjoint or other solver backends) are not implemented; see
[Follow-ups](#follow-ups).

All timings are on the M3 Max (12 P + 4 E cores, 40-core GPU).

## Stage 1: component tuning without new EM runs

### What it does

```sh
target/release/fieldsolve tune tune.json        # [--restarts N] [--rounds N] [--exhaustive 0|1|2] [--top N] [--no-outputs]
target/release/fieldsolve radmatrix run/farfield.json [candidates]   # radiation matrix vs combine + p_rad
```

`tune` loads `sparams.sNp`, `result.json` and `farfield.json` of an all-port EM run once. It then:

- **caches every allowed part's reflection** at the objective frequencies (one solve point per band frequency
  for the match, one for the radiation, one per harmonic; the frequency conventions are exactly `cosim`'s);
- **reduces the patterns to radiation matrices.** For embedded patterns E_i and incident waves a_i,
  `R_ij = ∮ (conj Eθ_i Eθ_j + conj Eφ_i Eφ_j) / (2 η0) dΩ` and `P_rad = Re(aᴴ R a)`, on the far field's
  own quadrature. This is algebraically `Pattern::combine(..).p_rad()`; the off-diagonal terms carry the
  ports' interference, so adding port powers would be wrong;
- **scores a candidate** with one small network solve per point (the terminated ports' `(I − S_tt Γ) x = S_te`,
  allocation-free on the stack) and a quadratic form. The score is `tools/ladder_retune.py`'s: per band
  frequency `(1 − |Γ_target|²) × P_rad / P_in`, worst of the band in dB, with optional harmonic limits
  (`"max_db": -20` or `"base"`, no worse than as built). Infeasible candidates rank below all feasible ones;
- **searches** by coordinate descent from as built and from random starts (positions in random order, best
  choice per position, until a round gains nothing), plus every one- and two-part change from as built;
- **writes detail only for the selected candidates.** The best `top` candidates (and as built) become a
  `cosim` config, which `cosim` runs for Touchstone files and combined patterns. The evaluator's scores
  are then checked against `cosim`'s output, and the largest difference is reported as
  `tune.json` `cosim_check`;
- **`"compare_em"`** re-scores the best candidates on another extraction of the same board (for example a
  screening run against the validated one) and reports score differences and rank agreement.

Choices are part specs (`"open"`, `"C:1p"`, `"s2p:file"`), named sets built from a vendor file template
(`"template": "lib/GRM1555C1H{code}{tol}A01.s2p"`, RKM codes `R50`, `1R2`, `4N7`, the first tolerance letter with
a file), filtered by value (`{"set": "caps", "min": "0.5p", "suffix": " ser"}`). The MTH config, the
equivalent of `ladder_retune.py`, is:

```json
{ "em": "mth7_board/sparams.s7p", "out": "tune_mth7", "drive": "CHIP", "target": "nordic_ref_network.s2p",
  "band": [2.40e9, 2.44e9, 2.48e9], "harmonics": [{"f": 4.88e9, "max_db": "base"}, {"f": 7.32e9, "max_db": "base"}],
  "sets": {"caps": {"template": "murata_sparams/lib_grm/GRM1555C1H{code}{tol}A01.s2p", "tols": "BCWDJ", "unit": "p",
                    "values": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.6, 6.8, 8.2]},
           "inds": {"template": "murata_sparams/lib_lqg/LQG15HS{code}{tol}02.s2p", "tols": "BCSJGH", "unit": "n",
                    "values": [1.0, 1.2, 1.3, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.7, 3.0, 3.3, 3.6, 3.9, 4.3, 4.7, 5.1, 5.6, 6.2, 6.8, 7.5, 8.2]}},
  "positions": {"C59": ["open", "@caps"], "L4": ["short", "@inds"], "C60": ["open", "@caps"], "L5": ["short", "@inds"],
                "R70": ["short", {"set": "inds", "max": "4.7n"}, {"set": "caps", "min": "0.5p", "suffix": " ser"}],
                "C20": ["open", {"set": "caps", "max": "2p"}]},
  "base": {"C59": "1p", "L4": "4.7n", "C60": "1.2p", "L5": "2.2n", "R70": "short", "C20": "open"},
  "search": {"restarts": 60, "rounds": 8, "exhaustive": 2, "seed": 7}, "top": 8 }
```

### Radiation matrix on the real board patterns

`fieldsolve radmatrix` on the MTH R7A 7-port far field (5° grid, 2664 directions), 2000 random complex
weight vectors per frequency, plus, per port, the weights of the other six that cancel that port's
radiation as deeply as possible (`a_rest = −R_rr⁻¹ R_ri`):

| f | max \|ΔP\| (W) | ΔP / incoherent sum | deepest cancellation | incoherent sum / true power, worst | matrix set-up | scoring speed-up |
|---|---|---|---|---|---|---|
| 2.40 GHz | 3.3e-15 | 3.1e-15 | L5 −62.9 dB | 3.1e6 | 0.2 ms | 782× |
| 2.44 GHz | 7.2e-15 | 1.9e-15 | L5 −61.2 dB | 2.3e6 | 0.2 ms | 777× |
| 2.48 GHz | 2.8e-14 | 5.6e-15 | L5 −61.0 dB | 2.3e6 | 0.7 ms | 1075× |
| 4.88 GHz | 1.5e-16 | 3.6e-15 | C20 −10.7 dB … L5 −31.6 dB | 2343 | 0.4 ms | 769× |
| 7.32 GHz | 1.9e-16 | 3.3e-15 | L5 −27.3 dB | 1178 | 0.4 ms | 800× |

Errors are relative to the incoherent sum `Σ |a_k|² R_kk`, the scale before cancellation. The two forms agree
to rounding, including where the combination cancels to −63 dB, and summing port powers instead would be
wrong by up to six orders of magnitude. The speed-up covers only the radiated-power scoring (`combine` +
`p_rad` against the quadratic form), as in the proposal's synthetic benchmark (775×).
`radiation_matrix_matches_combine_on_board_patterns` repeats this on a committed slice of the same run
(`src/bin/fieldsolve/testdata/mth7_mini`, 30° grid).

### Against cosim and the Python prototype

The evaluator reproduces `cosim` exactly: on the MTH extraction the 9-10 selected candidates of each search
differ from `cosim`'s output by at most 2.5e-14 dB, and `evaluator_matches_cosim` checks four part sets
(ideal, vendor-file and parasitic models, one per port) on the fixture to 1e-9 dB in the score, |Γ|, P_in,
P_rad and H2. As built scores −1.98 / −3.13 / −4.59 dB, the Python prototype's numbers, and the cheap fix
from [mth-antenna-findings.md](mth-antenna-findings.md) (L5 → 0 Ω, C20 = 0.5 pF) scores −1.64 / −1.43 /
−1.62 dB with H2 −25.9 dB, as documented there.

Benchmark on the MTH data, Python first (`tools/ladder_retune.py`: as built plus 5 random starts, 4 rounds,
each batch of up to 120 candidates through one `fieldsolve cosim` process), then `tune`:

| | candidates | wall time | candidates / s | best worst-in-band |
|---|---|---|---|---|
| `ladder_retune.py` (6 starts, 4 rounds) | 2111 unique, 85 `cosim` calls | 28.6 s | 74 | −1.39 dB |
| `tune`, the same search (6 starts, 4 rounds, no exhaustive pass) | 3592 | 0.004 s search, 0.11 s total | ≈ 9e5 | −1.31 dB |
| `tune`, default search (61 starts, 8 rounds, all 1- and 2-part changes) | 44 210 | 0.02-0.07 s search, 0.27 s total with outputs | 6e5-2.2e6 | −1.15 dB without harmonic limits, −1.36 dB with |
| `tune`, 2001 starts, 1 thread | 1.1 M | 4.6-5.7 s | 1.9-2.4e5 | |
| `tune`, 2001 starts, 2 threads | 1.1 M | 3.0-3.2 s | 3.4-3.7e5 | |
| `tune`, 2001 starts, 4 / 8 / 12 / 16 threads | 1.1 M | 1.8-2.4 s | 4.6-6.1e5 | |

Two runs per row, on the 7-port MTH config without harmonic limits. Other sessions shared the machine (load
average about 11), so the multi-thread rows are a lower bound. Throughput levels off at about 4 threads,
and in an earlier test, eight independent single-thread processes levelled off the same way, at about
8e4 /s each, so the shared machine sets the ceiling, not the search's synchronisation. The figure that
matters is 1.1 M candidates in 2-6 s, where the Python search needed half a minute for 2000.

Per candidate that is about 13.6 ms in Python (process start-up, JSON, re-reading the 7.5 MB far field and
writing a pattern file per candidate) against 4-6 µs single-threaded, a factor of about 3000; the search
itself is no longer the bottleneck. Stage timings of the default run (`tune.json` `timings_s`): load
0.045 s (Touchstone, `result.json`, the 7.5 MB `farfield.json`), parts 0.05 s (163 vendor files parsed and
their reflections cached), radiation matrices 0.002 s, search 0.03 s, outputs 0.09 s (`cosim` on the
selected candidates).

### MTH results

| | 2.40 GHz | 2.44 GHz | 2.48 GHz | radiated H2 / H3 per W into the chip pin |
|---|---|---|---|---|
| as built | −1.98 | −3.13 | −4.59 | −25.5 / −21.6 dB |
| best, H2 and H3 no worse than as built: C59 open, L4 1.0 nH, C60 0.8 pF, L5 1.5 nH, R70 0 Ω, C20 open | −1.30 | −1.19 | −1.36 | −26.1 / −27.6 dB |
| best ignoring harmonics: C59 0.7 pF, L4 1.0 nH, C60 0.5 pF, L5 0 Ω, R70 0 Ω, C20 open | −1.10 | −0.99 | −1.14 | −21.3 / −20.2 dB |

The constrained optimum is about 0.25 dB better in the band than the result the Python search reported
(−1.6 dB worst) and lowers H3 by 6 dB, with three changed parts. As before, all of these are free-space
values for the untrimmed antenna.

## Stage 2: fidelity policies for geometry changes

> **Provisional.** The calibration below is without the edge-singularity and conformal-metal corrections
> that the `fieldsolve-adaptive-mesh` branch adds to the model. Most of the screening error is staircase error
> on coarse cells, which those corrections are meant to reduce, so the policies and these tables need a
> re-run with the corrections on once that branch is merged (`tools/fieldsolve_calib/`, below). A level
> switches them on through `"fidelity": {"set": {"screening": {"<spec key>": value}}}`, which is passed into
> the resolved spec and recorded.

### What it does

```sh
target/release/fieldsolve run spec.json --fast                  # = --fidelity screening
target/release/fieldsolve run spec.json --fidelity refined      # or validated; a plain run is level "spec"
```

A policy changes the numerical effort of the spec and nothing of its physics. The whole board stays in the
domain, since the ground plane is part of the antenna, and materials, loss, ports, PML and air are untouched.

| level | fine cell | protected features | outputs | stop |
|---|---|---|---|---|
| `spec` (plain `run`) | as given | | as given | as given (energy / port decay) |
| `validated` | as given | | as given | objective spread ≤ 0.1 dB over 8 periods |
| `refined` | × 1.5 | spec cell | no fields, frames or scene; far field only in the band and harmonics | objective spread ≤ 0.1 dB over 6 periods |
| `screening` (`--fast`) | × 2 | spec cell | as `refined` | objective spread ≤ 0.2 dB over 4 periods |

- **Protected features** keep the spec's cell size through `"refine"` boxes. These cover every port (0.25 mm
  margin) and the spec's `"fidelity": {"protect": [...]}` list (pad names, boxes), for the feed, gaps and
  shorting stubs. `"refine"` is a new spec key and is usable without a policy.
- **The objective** is the spec's `"objective": {"port", "band", "harmonics"}`, with the band defaulting to
  the far-field frequencies. Per band frequency, a run records the driven port's accepted power, its coupling
  to the other ports and the radiated box flux, each in dB, at every energy check. `result.json`
  `"objective"` holds the per-frequency totals: radiated power per watt incident, with the other ports
  matched (the radiated part from the far field), and the worst of the band.
- **The objective stop** ends an excitation once every component's spread (max − min) over the last
  `window_periods` periods of the band's lowest frequency is at most `tol_db`, and the source has ended,
  followed by `min_settle_periods` periods (8 by default). The standard energy/port rule is then only
  recorded (`standard_rule_step`). A run that reaches `max_steps` makes the record's `status` `incomplete`.
- **Everything is recorded** in `result.json` `"fidelity"`:
  - the level and policy, and every override as `{key, from, to, why}`;
  - the user's own policy numbers (`"fidelity": {"cell_scale", "tol_db", ...}` in the spec win over the level's);
  - `actual`: dims, cells, dt, minimum step per axis, the refine boxes, excitations, steps, seconds, stop
    reasons, objective worst and stability.

  `"mesh"` and `"objective"` are recorded for every run, a plain one included.

### Why the stop rule looks at the objective's spread

A 150 000-step full-mesh run of the MTH PIFA (`base_long`) shows how the radiated power settles. Accepted
power converges early. The radiated-power DFT first creeps up, then beats against a slowly decaying
resonance near 2.67 GHz (late Q ≈ 75, likely a cavity mode between the planes). On the worst in-band value,
the beat's half amplitude is 0.034 dB, period about 4.4 ns. A tolerance near 0.01 dB is therefore unreachable
in reasonable time: `validated` with a 0.01 dB tail estimate ran into `max_steps` at 150 k steps. A spread
over a window, by contrast, can be reached. The standard full-run stop (ports −60 dB and flux within 1 % over
5 checks) fires at 55 040 steps, 0.005 dB from the beat's centre.

Stop rule replayed on the `base_long` history (`analyze.py --replay`), with the error against the centre of
the late beat:

| window | tol 0.1 dB | tol 0.2 dB | tol 0.3 dB |
|---|---|---|---|
| 2 periods | 1.02 × standard, −0.005 dB | 0.89 ×, −0.048 dB | 0.78 ×, −0.128 dB |
| 3 periods | 1.13 ×, +0.009 dB | 1.04 ×, +0.003 dB | 0.95 ×, −0.020 dB |
| **4 periods** | 1.20 ×, +0.016 dB | **1.11 ×, +0.008 dB** (screening) | 1.04 ×, +0.003 dB |
| 6 periods | **1.38 ×, +0.026 dB** (refined) | 1.29 ×, +0.023 dB | 1.22 ×, +0.022 dB |
| 8 periods | **1.96 ×, −0.017 dB** (validated) | 1.46 ×, +0.024 dB | 1.38 ×, +0.026 dB |

So on this board no stop rule saves steps against the standard one without losing accuracy. All of
screening's saving comes from the mesh, and the stop rule's job is to be robust (it does not stop on a quiet
moment of the port signal). The longer windows of `refined` and `validated` cost 1.4–2 × the steps and buy
nothing here, since the beat dominates.

### Calibration on twelve MTH geometries

Twelve candidates of the MTH R7A PIFA (antenna alone at the C20 feed, whole board, three band frequencies)
come from the as-built Gerbers through the new spec key `"edits"`
(`[{"layer", "op": "add" | "clear", "rect", "net"}]`):

- the radiator's open end trimmed by 1.1 and 2.0 mm, or extended by 1.0 and 1.5 mm;
- the feed moved 0.5 or 1.0 mm towards the open end, or 0.3 mm away from it;
- the ground clearance under the radiator widened by 0.5 or 1.0 mm, or narrowed by 0.5 mm;
- a trim combined with a feed move.

Every candidate ran at `spec` level (the reference: 12.85 M cells, dt 88.1 fs, standard stop) and at
`screening` (4.45 M cells, dt 100.5 fs); six also ran at `refined` (6.60 M cells, dt 97.2 fs). The worst
in-band radiated power per watt incident, in dB:

| candidate | full | screening | error | refined | error | screening cells × steps saving |
|---|---|---|---|---|---|---|
| base | −0.914 | −1.321 | −0.41 | −1.375 | −0.46 | 3.23 × |
| trim 1.1 mm | −1.685 | −1.404 | +0.28 | | | 3.21 × |
| trim 2.0 mm | −3.220 | −3.002 | +0.22 | −2.921 | +0.30 | 2.89 × |
| extend 1.0 mm | −3.013 | −3.826 | −0.81 | | | 2.94 × |
| extend 1.5 mm | −4.495 | −5.559 | −1.06 | −5.825 | −1.33 | 2.44 × |
| feed −0.5 mm | −0.909 | −0.780 | +0.13 | | | 3.10 × |
| feed −1.0 mm | −1.657 | −1.386 | +0.27 | −1.359 | +0.30 | 3.62 × |
| feed +0.3 mm | −1.422 | −2.156 | −0.73 | | | 3.29 × |
| clearance +0.5 mm | −1.080 | −0.880 | +0.20 | | | 3.40 × |
| clearance +1.0 mm | −1.727 | −1.500 | +0.23 | −1.499 | +0.23 | 2.83 × |
| clearance −0.5 mm | −4.269 | −5.407 | −1.14 | −5.507 | −1.24 | 2.63 × |
| trim 1.1 + feed −0.5 | −2.458 | −2.195 | +0.26 | | | 3.04 × |

| | screening (n = 12) | refined (n = 6) |
|---|---|---|
| error mean / sd | −0.21 / 0.57 dB | −0.37 / 0.77 dB |
| largest error (worst in band / any band frequency) | 1.14 / 1.17 dB | 1.33 / 1.33 dB |
| Spearman rank correlation with full | 0.944 | 0.943 |
| candidate pairs in the same order | 61 / 66 | 14 / 15 |
| pairs more than 2 sd (1.15 dB) apart in the same order | 36 / 36 | 8 / 8 |
| same best candidate | yes | no: feed −1.0 mm by 0.015 dB over base, 0.74 dB behind it in full |
| cells × steps saving, median (range) | 3.10 × (2.44–3.62) | 1.35 × (1.27–1.57) |
| cells / time step / steps | ÷ 2.89 / × 1.14 / 0.94 × | ÷ 1.95 / × 1.10 / 1.44 × |
| wall-time saving, median (contended GPU) | 2.7 × | 1.05 × |

What the numbers say:

- **The error is a resonance shift, not noise.** The coarse staircase moves the PIFA's resonance down by
  about 1 %. A candidate whose worst band edge is the upper one (base, the extensions, clearance −0.5 mm,
  feed +0.3 mm) therefore looks worse, and one limited at the lower edge looks better. The signed error is
  as large as 1.1 dB, but it is systematic, so screening ranks well: every pair more than 1.15 dB apart in
  the reference is in the right order, and the five swapped pairs are all within 0.3 dB of each other
  in the reference.
- **`refined` does not improve on `screening` on this board.** Its errors have the same sign and size
  (1.5 × cells still staircase the radiator edge), and its 6-period window costs 1.44 × the steps, so it
  saves only 1.35 × against the full run. Until the model corrections are in, the useful ladder is
  screening for search and the full mesh (`spec` or `validated`) for the candidates screening picks. When
  re-calibrating, the refined window should drop to 4 periods (the replay above).
- **A screening result is a ranking signal with a ±1.2 dB absolute uncertainty.** Candidates within about
  1.2 dB of each other in screening need a full run to separate them.

### 7-port extraction for component tuning

The 7-port MTH run (all ladder positions plus the chip pin, band and two harmonics) at `screening`, fed to
`tune` with `"compare_em"`, against the full extraction (`mth7_board`, 10.8 M cells):

| | full extraction | screening extraction |
|---|---|---|
| cells, steps (7 excitations) | 10.83 M, 458 752 | 5.53 M, 502 528 (cells × steps 1.79 × fewer) |
| as built, worst in band | −4.59 dB | −4.24 dB |
| best (harmonics no worse than as built) | −1.36 dB | −1.22 dB (C59 open, C60 0.9 pF, L4 0 Ω, L5 1.2 nH, R70 4.3 pF series) |
| that best scored on the full extraction | | −1.42 dB, rank 2 of the full search's 41 candidates |
| 41 best candidates re-scored on the other extraction | | mean −0.20 dB, max 0.35 dB, Spearman 0.83 |

Screening on the 7-port saves less: protecting seven ports keeps more cells fine (÷ 1.96), and the
objective stop over 7 ports × 5 frequencies with coupling terms runs 1.10 × the full run's steps. The
screening extraction's best ladder is within 0.06 dB of the full one's, so for component tuning a
screening extraction plus `tune`, then one full run of the chosen ladder, would have found the same answer.

### Wall time

Measured while another session's FDTD queue shared the GPU, so the wall-time savings above are contended:
median 2.7 × for screening, 1.6–4.7 × per candidate. The reference runs took 458–1128 s and the
screening runs 202–298 s, with the spread coming from the other jobs. The machine-independent figure is cells × steps
(median 3.1 ×; the step count already reflects the 1.14 × longer time step). With the Metal kernels' throughput roughly
constant in cells per second, that is the wall-time saving to expect on an idle GPU. No idle window came
up during this work to confirm it.

### Reproducing

`tools/fieldsolve_calib/`: `make_specs.py MTH_7PORT_SPEC GERBERS` writes the twelve candidate specs,
`run_all.sh` runs full, screening and refined runs one GPU job at a time (`FS=` picks the binary), and
`analyze.py [--replay RUN]` prints the tables above and writes `calib.json`. About 3.5 h of GPU time on the
M3 Max.

## Follow-ups

- **Re-calibrate with the model corrections** from `fieldsolve-adaptive-mesh` (edge singularity, conformal
  metal), switched on per level via `fidelity.set`. Then decide whether `refined` earns its place, and drop
  its window to 4 periods.
- **Stage 3, surrogate-assisted outer search:** screening runs as the cheap model in a geometry search, with
  full runs for the shortlist. The ±1.2 dB screening uncertainty sets the shortlist width.
- **Stage 4:** adjoint sensitivities, or other solver backends (MoM or FEM for small structures).
- **A single-feed screening option for multiport tuning extractions**, since all-port screening saves only
  1.8 × in cells × steps.
- **A screening-level bias correction:** subtracting a per-board offset measured from one screening/full pair
  (base: −0.41 dB) would centre the error. The sign depends on which band edge is worst, though, so a
  frequency-shift correction is the better model.
