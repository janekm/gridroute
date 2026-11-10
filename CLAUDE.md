# gridroute

## fieldsolve (`src/bin/fieldsolve`, docs in `docs/fieldsolve*.md`)

- **Read `docs/fieldsolve-settings.md` before writing or reviewing a `fieldsolve run` spec.** It lists the
  keys to use for each job and the mistakes that have already cost results.
- **Antenna and far-field runs:** `"preset": "antenna"`. That gives surface-impedance copper
  (`"conductor"`), wideband FR4 loss, sub-cell edges, mesh refinement and the whole board in the domain.
  Never use `"loss": {"copper": true}`: it is a plane-pair stand-in that leaves the copper lossless.
- **`"roi"` must cover every piece of copper carrying RF current,** the whole antenna included. Do a 20-step
  dry run first (`"max_steps": 20`), read the `model:` line and every `WARNING:` line, and fix them or
  acknowledge them (`"coarse_copper_ok"`) before starting a long run.
- **Before acting on a number, do a two-mesh check:** the run again at about 0.7 × `cell`.
- **Tests:** `RAYON_NUM_THREADS=4 cargo test --release --bin fieldsolve -- --test-threads=1`. Parallel test
  threads contend for the GPU and fail spuriously. So does sharing the GPU with another job, for timing
  numbers.
