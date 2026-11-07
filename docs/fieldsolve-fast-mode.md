# fieldsolve fast iteration: component tuning and fidelity policies

November 2026. This implements the
[fast-mode proposal](../reports/pcb-field-solver-review-2097498/fast-mode/proposal.md) (outside the repository,
in `reports/pcb-field-solver-review-2097498/fast-mode/`) in its four stages:

1. **Component changes reuse one EM extraction.** `fieldsolve tune` scores part combinations on a saved
   N-port result, with no new EM run.
2. **Geometry changes get cheaper EM runs with a measured accuracy budget.** `fieldsolve run SPEC --fast`
   resolves to a recorded fidelity policy (`screening`, `refined` or `validated`).
3. **Geometry search spends full solves only where they matter.** `fieldsolve optimise SPACE.json` searches a
   few geometry parameters on screening solves, maps screening onto full fidelity from a few promoted
   full solves, and returns the best full-fidelity candidate.
4. **The next solver investment, chosen from a profile.** Where a screening run's time goes, and what adjoint,
   frequency-domain or local time stepping could and could not save ([Stage 4](#stage-4-where-the-remaining-time-goes)).

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
> the resolved spec and recorded. On that branch the corrections are the `"edges"` spec key, and the
> recommended per-level setting is:
>
> ```json
> "edges": true,
> "fidelity": {"set": {"screening": {"edges": {"dt_check": false}}, "refined": {"edges": true}, "validated": {"edges": true}}}
> ```
>
> Screening skips the CPU stability check (about 100 s on large boards); it could also leave `"edges"` out.
> `make_specs.py --edges` writes the calibration specs this way.

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

- the radiator's open end (x = 141.9) trimmed by 1.0 and 1.94 mm, or extended by 1.0 and 1.5 mm;
- the feed moved 0.5 or 1.0 mm towards the open end, or 0.3 mm away from it;
- the ground clearance under the radiator widened by 0.5 or 1.0 mm, or narrowed by 0.5 mm;
- a trim combined with a feed move.

**The as-built model differs from the earlier MTH reference.** These specs widen the fine-mesh ROI from
x = 145 to 141.5, so the radiator's open end (x 141.9–144.35) is meshed at 0.05 mm. In the reference run of
[mth-antenna-findings.md](mth-antenna-findings.md) it fell in graded cells of about 0.35–0.5 mm (about 5 nodes
along the 2.45 mm arm), and that coarse open end lowers the resonance. With the reference ROI, the same
one-port spec reproduces the reference's antenna-alone impedance to about 1 Ω. The widened ROI is the
better-resolved as-built model; its own mesh convergence is not demonstrated.

| C20 into 50 Ω, as built | 2.40 GHz | 2.44 GHz | 2.48 GHz | worst in band |
|---|---|---|---|---|
| reference 7-port, other ports open | 37.5+j5.1 Ω | 28.5+j16.1 Ω | 21.8+j26.3 Ω | |
| one-port, reference ROI (x ≥ 145) | 38.3+j5.6 Ω | 28.3+j17.2 Ω | 21.9+j27.7 Ω | −2.09 dB |
| one-port, widened ROI (x ≥ 141.5) | 59.8−j8.8 Ω | 44.3+j2.3 Ω | 33.4+j13.2 Ω | −0.91 dB |

Every candidate ran at `spec` level (the reference: 12.85 M cells, dt 88.1 fs, standard stop) and at
`screening` (4.45 M cells, dt 100.5 fs); six also ran at `refined` (6.60 M cells, dt 97.2 fs). The worst
in-band radiated power per watt incident, in dB:

| candidate | full | screening | error | refined | error | screening cells × steps saving |
|---|---|---|---|---|---|---|
| base | −0.914 | −1.321 | −0.41 | −1.375 | −0.46 | 3.23 × |
| trim 1.0 mm | −1.685 | −1.404 | +0.28 | | | 3.21 × |
| trim 2.0 mm | −3.220 | −3.002 | +0.22 | −2.921 | +0.30 | 2.89 × |
| extend 1.0 mm | −3.013 | −3.826 | −0.81 | | | 2.94 × |
| extend 1.5 mm | −4.495 | −5.559 | −1.06 | −5.825 | −1.33 | 2.44 × |
| feed −0.5 mm | −0.909 | −0.780 | +0.13 | | | 3.10 × |
| feed −1.0 mm | −1.657 | −1.386 | +0.27 | −1.359 | +0.30 | 3.62 × |
| feed +0.3 mm | −1.422 | −2.156 | −0.73 | | | 3.29 × |
| clearance +0.5 mm | −1.080 | −0.880 | +0.20 | | | 3.40 × |
| clearance +1.0 mm | −1.727 | −1.500 | +0.23 | −1.499 | +0.23 | 2.83 × |
| clearance −0.5 mm | −4.269 | −5.407 | −1.14 | −5.507 | −1.24 | 2.63 × |
| trim 1.0 + feed −0.5 | −2.458 | −2.195 | +0.26 | | | 3.04 × |

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
constant in cells per second, that is the wall-time saving to expect on an idle GPU.

An idle re-time of the as-built pair confirms it. No other FDTD job ran during it, apart from two 20-step dry runs:

| | cells × steps | stepping | total wall time | stepping rate |
|---|---|---|---|---|
| full (`spec`) | 12.85 M × 55 040 | 253.7 s | 258 s | 2788 Mcell/s |
| screening | 4.45 M × 49 152 | 75.3 s | 77 s | 2908 Mcell/s |
| ratio | 3.23 × | 3.37 × | 3.35 × | |

On an idle GPU, the wall-time saving matches the cells × steps saving: 3.35 × against 3.23 ×. The small
excess comes from the far-field transform, which shrinks with the box.

### Reproducing

`tools/fieldsolve_calib/`: `make_specs.py MTH_7PORT_SPEC GERBERS` writes the twelve candidate specs,
`run_all.sh` runs full, screening and refined runs one GPU job at a time (`FS=` picks the binary), and
`analyze.py [--replay RUN]` prints the tables above and writes `calib.json`. About 3.5 h of GPU time on the
M3 Max.

## Stage 3: geometry search with adaptive promotion

### What it does

```sh
target/release/fieldsolve optimise space.json      # [--method surrogate|pattern] [--budget N] [--out DIR]
target/release/fieldsolve optimise-bench           # [--seeds N] [--promote-every N] [--final N] [--json OUT]
```

`space.json` names a spec and a few parameters. The spec is parameterised by expressions: any string
starting with `=` is evaluated with the parameters (numbers, names, `+ - * /`, parentheses, `min`, `max`,
`abs`, `clamp`), and an edit whose rectangle comes out empty is skipped, so one spec covers both trimming
and extending. The MTH space, three parameters on a 0.1 mm grid:

```json
{ "spec": "mth_param.json", "out": "surrogate", "cache": "cache", "budget": 16,
  "parameters": [{"name": "arm",  "min": -1.9, "max": 1.5, "step": 0.1, "start": 0},
                 {"name": "feed", "min": -1.0, "max": 0.5, "step": 0.1, "start": 0},
                 {"name": "clr",  "min": -0.5, "max": 1.0, "step": 0.1, "start": 0}] }
```

with edits such as `{"layer": "F.Cu", "op": "add", "rect": ["=141.9 - max(arm, 0)", 74.64, "=141.9 + min(max(arm, 0), 0.05)", 74.94]}`
for the radiator's open end. `arm` lengthens (positive) or trims the radiator. `feed` moves the feed strip
towards the open end (negative) or towards the short. `clr` widens (positive) or narrows the ground
clearance under the radiator. All values are mm.

- **Screening solves explore.** The default method, `surrogate`, does the following:
  - starts from the start point and a Latin hypercube (2d + 1 points);
  - then screens the point of largest expected improvement of a Gaussian-process model, inside a trust
    region around the incumbent. The trust region doubles after an improvement and halves after repeated
    failures.

  `pattern` is a compass search from the start point. It is the simple bounded search the proposal asks to
  compare against, with the same budget and promotion rules.
- **Full solves are promotions:**
  - every `promote_every` (6) screening solves, the best predicted full-fidelity point gets a full solve;
  - at the end, up to `final` (3) more, while the prediction is within `tie_db` (1.2 dB, the calibration's
    2 sd) of the best full value.

  The answer is the best full-fidelity point.
- **Space mapping carries screening onto full fidelity.** The screening model is a Gaussian process per
  band frequency. Full fidelity is modelled as that model at shifted parameters plus an offset:
  `full_f(u) ≈ screen_f(u + s) + b`. The shift `s` and offset `b` are fitted by ridge-regularised least
  squares to every full solve's per-frequency values, so one full solve already constrains the shift.
  Expected improvement uses the mapped prediction, and the solve is made at the point's screening image
  `u + s`. `mapping: false` promotes the best screening candidates as they are.
- **Every solve is cached.** `cache/LEVEL/KEY/` holds `result.json` and the resolved `spec.json`. A
  repeated point is never solved twice, two searches can share a cache, and an interrupted search resumes.
- **The record** `optimise.json` holds:
  - every evaluation (level, point, per-frequency values, cells × steps, seconds, why it was chosen);
  - the trajectory of the best full value;
  - screening and full solve counts and costs;
  - the fitted mapping;
  - the full-minus-screening differences.

### Why space mapping

The calibration shows a systematic error: screening shifts the resonance down by about 1 %. On a sharp
worst-in-band objective, the screening optimum is then a different geometry, and in the synthetic test
below its full-fidelity value is about 1.3–1.5 dB below the true optimum. Any search that trusts
screening values as they are inherits that gap, however many solves it spends. A constant offset model
cannot remove it, because the error depends on where the band edge sits. A parameter shift is the right
model, and the per-frequency responses identify it from very few full solves.

### Synthetic benchmark

`optimise-bench` runs the methods on a two-fidelity test problem built to match the measured behaviour:
- three parameters, each on a 41-point grid;
- a resonance moved by the parameters (strongly by the first, like a radiator length), with a match factor
  on the others;
- worst of three band frequencies;
- the screening level shifts the resonance down 1 % and adds 0.03 dB of repeatable noise;
- cost is counted in full-solve equivalents, a screening solve being 1 / 3.1 (the measured cells × steps
  saving).

30 seeds; regret is the true full-fidelity optimum minus the best full value found:

| method | screening solves | full solves | cost (full solves) | mean regret | median | within 0.1 dB |
|---|---|---|---|---|---|---|
| surrogate + mapping | 12 | 4.0 | 7.9 | 0.29 dB | 0.22 dB | 33 % |
| surrogate + mapping | 16 | 5.0 | 10.2 | 0.16 dB | 0.09 dB | 53 % |
| surrogate + mapping | 24 | 6.0 | 13.7 | 0.13 dB | 0.08 dB | 60 % |
| surrogate + mapping | 39 | 8.9 | 21.6 | 0.04 dB | 0.03 dB | 90 % |
| pattern + mapping | 16 | 5.0 | 10.2 | 0.13 dB | 0.10 dB | 53 % |
| surrogate, no mapping | 16 | 4.9 | 10.1 | 1.26 dB | 1.43 dB | 0 % |
| pattern (no mapping) | 16 | 5.0 | 10.1 | 1.34 dB | 1.29 dB | 0 % |
| pattern (no mapping) | 36 | 8.4 | 19.9 | 1.53 dB | 1.55 dB | 0 % |
| surrogate, full fidelity only | 0 | 12 | 12 | 0.32 dB | 0.30 dB | 23 % |
| surrogate, full fidelity only | 0 | 16 | 16 | 0.12 dB | 0.07 dB | 67 % |
| pattern, full fidelity only | 0 | 24 | 24 | 0.26 dB | 0.21 dB | 20 % |

- **Mapping is what makes screening useful.** Without it, both search methods stay 1.3–1.5 dB from the
  optimum, and more screening solves do not help.
- **With mapping, screening plus promotion beats searching at full fidelity, by a modest margin.** At about
  10 full-solve equivalents it reaches 0.16 dB mean regret, which a full-fidelity surrogate needs about
  14–16 full solves to match. That is roughly a third less EM. It is not an order of magnitude: one
  screening solve costs a third of a full one, and a few full solves are needed to fit the mapping.
- **The surrogate beats the compass search at full fidelity** (0.12 dB at 16 solves against 0.26 dB at 24).
  With mapping on, the two are close at a budget of 16, because the mapped promotions do much of the work.
- **Promotion settings:** final 3 and promote every 6 are the defaults. In a sweep of promote every
  {4, 6} × final {2, 3}, the variants differed by about 0.05 dB in mean regret, against a cost difference of
  one or two full solves.

### MTH PIFA search

A surrogate search with the space above (budget 16), on the widened-ROI MTH model of the calibration:

| | screening solves | full solves | cells × steps | full-run equivalents | best full value |
|---|---|---|---|---|---|
| surrogate + mapping | 16 | 5 | 7.10e12 (3.56 screening + 3.54 full) | 10.0 | −0.467 dB |
| pattern (no mapping), the simple bounded search | 16 | 5 | 7.19e12 (3.55 screening + 3.64 full) | 10.2 | −0.709 dB |

At the same budget and nearly the same cost, the compass search ended 0.24 dB worse, at (−0.6, +0.4, 0).
It followed raw screening values, which were wrong by −0.27 to +0.34 dB on the points it promoted. Its
first poll, `arm` +0.8 alone, screened at −3.68 dB, so it turned towards shorter arms. The surrogate's
optimum lies the other way: a longer arm only pays off together with a wider clearance. Its best screening
point, (−0.8, +0.4, 0), scored −0.75 dB at screening and −0.85 dB at full fidelity. This is one run of each
method on one problem, so the synthetic benchmark above is the basis for comparing them. The MTH run agrees
with it.

One full run of the as-built model is 7.07e11 cells × steps. At the idle rates of Stage 2 (77 s per screening
run, 258 s per full run), a 16 + 5 search takes about 42 min of GPU time, against about 69 min for 16 full
solves. Trajectory of full solves (the as-built model
scores −0.914 dB):

| after | full solve at (arm, feed, clr) mm | predicted | full | why |
|---|---|---|---|---|
| 6 screening | (0, 0, 0), as built | −1.33 | −0.914 | periodic promotion |
| 12 screening | (+0.3, +0.3, +0.5) | −0.31 | −0.579 | periodic promotion |
| 16 screening | (+1.5, −0.3, +1.0) | −0.43 | −0.532 | final |
| | (+1.5, −0.2, +1.0) | −0.45 | −0.487 | final |
| | **(+1.4, −0.3, +0.9)** | −0.46 | **−0.467** | final |

At the end, the fitted mapping was a shift of +0.11 mm in `arm` and +0.19 mm in `clr`, with a 0.02 dB offset
and 0.035 dB rms residual. Its last prediction was within 0.004 dB of the full solve, and the shift is in the
direction the calibration predicts: the coarse mesh's lower resonance has to be compensated by a shorter
radiator in screening.

C20 impedance into 50 Ω and the objective's parts (full fidelity):

| | 2.40 GHz | 2.44 GHz | 2.48 GHz | radiated / accepted |
|---|---|---|---|---|
| as built | 59.8−j8.8 Ω | 44.3+j2.3 Ω | 33.4+j13.2 Ω | −0.43 / −0.48 / −0.58 dB |
| best (+1.4, −0.3, +0.9) | 79.5−j6.2 Ω | 62.6+j1.1 Ω | 49.3+j10.3 Ω | −0.25 / −0.26 / −0.37 dB |

About half the gain is radiation: the longer arm and wider clearance raise the share of the accepted power
that is radiated by 0.2 dB. The other half is the match, which is flatter across the band.

Caveats:
- **Antenna-alone values only.** These are free-space values for the antenna alone at the C20 feed into
  50 Ω, with the ladder open. A geometry change is evaluated properly together with a retuned ladder:
  `tune` on a 7-port extraction of the chosen geometry.
- **Two parameters are at the edge of their range.** The best point has `arm` and `clr` one step inside
  their upper bounds, so a wider range may do better.
- **Mesh convergence of the widened-ROI model is not demonstrated** (calibration note above).
- **Fabrication tolerance.** Each move is a fraction of a millimetre, so the result should be checked at
  ±0.1 mm on each parameter before the last tenth of a dB is trusted.

## Stage 4: where the remaining time goes

The proposal's fourth stage is a decision: adjoint sensitivities, another solver, or local time stepping,
"only after profiling identifies the remaining dominant cost". The profile, on the MTH screening run of the
as-built board (4.45 M cells, 49 152 steps), on an idle GPU:

| phase | time | share |
|---|---|---|
| set-up and mesh build | 0.4 s | 0.5 % |
| time stepping (Metal) | 75.3 s | 97.8 % |
| far-field transform (CPU, 1.0 M face samples, 3 frequencies) | 1.0 s | 1.3 % |

The time stepping is the field-update kernel itself, not the bookkeeping around it. Each run below is 8000
steps on the same screening mesh, on an idle GPU:

| configuration | rate |
|---|---|
| `fieldsolve bench` (update, PML, source, conductor loss, port probes; no DFTs, no checks) | 3080 Mcell/s |
| screening spec, no objective, no far-field box | 3043 Mcell/s |
| + far-field box (gathered DFT every 414 steps) | 3033 Mcell/s |
| + objective, evaluated every 1000 steps (port spectra and box flux on the CPU) | 3020 Mcell/s |

The accumulators and the objective checks cost about 2 % together. The whole 49 152-step run reaches
2908 Mcell/s, 6 % below the bench. Nothing outside the update is worth optimising first.

(An earlier profile of the same run, taken while another session's queue shared the GPU, showed 843 Mcell/s
and a 5 % far-field share. Both were contention.)

Inside the stepping, the step count is the product of two things:
- **Physical time:** 4.9 ns, set by the objective's settling against the slowly decaying 2.67 GHz
  resonance (Stage 2).
- **The time step:** 100.5 fs, fixed by the protected 0.047 mm cells around the feed.

Without the protected regions, the same screening mesh would have 3.44 M cells and a 165.4 fs step. The
protection therefore costs a factor of 1.29 × 1.65 = 2.1 in cells × steps. On a tensor-product grid, its
fine lines run through the whole domain, and its smallest cell sets the step everywhere.

What each option could and could not buy:

- **Adjoint sensitivities** give a gradient for about two solves, whatever the number of parameters.
  - With three to five parameters they do not pay: the MTH search found its optimum in about 10 full-run
    equivalents, and a finite-difference gradient costs d + 1 solves.
  - They also need geometry that is differentiable through the rasterisation. Today copper is binary per
    cell edge, so an adjoint would first need the conformal treatment of the adaptive-mesh branch.
  - They become the right tool for free-form copper (tens to thousands of pixel or boundary variables),
    which this work does not attempt.
- **Frequency-domain solves.** The objective needs three frequencies. A screening run's 49 000 time steps,
  shared between three 13 M-unknown systems with one matrix-vector product costing roughly one time step,
  are a budget of about 16 000 Krylov iterations per frequency.
  - That is reachable only with a good preconditioner for an open, PML-terminated, lossy problem, which is a
    research project in itself.
  - It would remove the high-Q ringing that sets the physical time.
  - Not recommended as the next step.
- **Local time stepping or subgridding of the protected regions** has a measured upper bound of about 2.1×
  on screening runs.
  - It is the largest single factor left, and it compounds with stages 2–3.
  - It needs a stable interface between grids (the proposal's caveat), so it is substantial work.
  - The full-fidelity run has 0.05 mm cells throughout its ROI, so subgridding only the feed would not help
    there.
- **Late-time extrapolation** of the objective, fitting the decaying resonances (as Harminv does) to stop
  before the beat settles, is cheaper to try. The 150 000-step history of Stage 2 is the test case: the
  standard stop fires at 55 k steps, 0.005 dB from the beat centre. An extrapolation that is right to
  0.05 dB at 25–30 k steps would halve every run.

**Recommendation:**
1. **Late-time extrapolation** first: an offline experiment on the saved objective histories, then a stop
   rule. It is small, and it halves the step count of every run if it holds.
2. **Then subgridding** of protected features for screening (up to 2.1×).
3. **Adjoint** only if free-form geometry becomes a goal.
4. **No frequency-domain backend** for now.

None of these is implemented here.

## Follow-ups

- **Re-calibrate with the model corrections** from `fieldsolve-adaptive-mesh` (`make_specs.py --edges`),
  then decide whether `refined` earns its place and drop its window to 4 periods.
- **Mesh convergence of the widened-ROI MTH model** (a 0.035 mm run, or the corrections on), before
  geometry results are compared with measurements.
- **MTH geometry:**
  - widen the `arm` and `clr` ranges;
  - check fabrication tolerance (±0.1 mm) at the best point;
  - take a 7-port extraction of it for `tune`, so the ladder is retuned for the new geometry.
- **Late-time extrapolation** of the objective (Stage 4), offline on the saved histories first.
- **Subgridding or local time stepping** of protected features (up to 2.1× on screening).
- **A single-feed screening option for multiport tuning extractions**, since all-port screening saves only
  1.8× in cells × steps.
