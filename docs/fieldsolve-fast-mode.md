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
| `tune`, default search (61 starts, 8 rounds, all 1- and 2-part changes) | 44 210 | 0.03-0.07 s search, 0.27 s total with outputs | 6e5-1.6e6 | −1.15 dB |
| `tune`, 2001 starts, 1 thread | 1.1 M | 4.2-6.7 s | 1.6-2.8e5 | |
| `tune`, 2001 starts, 4 / 8 / 16 threads | 1.1 M | 1.4-2.1 s | 5.5-8.3e5 | |

STAGE1_SCALING

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


STAGE2_PLACEHOLDER
