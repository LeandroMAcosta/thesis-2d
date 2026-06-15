# CLAUDE.md — thesis-2d

This file provides guidance to Claude Code (claude.ai/code) when
working in this repository.

## What this repo is

CPU + OpenMP implementation of the GCAS gas simulation in **2D**,
using **event-driven analytical integration**. Spinoff of the
[`thesis`](https://github.com/LeandroMAcosta/thesis) parent repo (1D).

## Quick facts

- Single backend: CPU + OpenMP. **No CUDA** (event-driven scales worse
  on GPU due to per-particle variable workload — see parent repo's
  `event-driven/docs/ALGORITHM.md` §9).
- State per particle: $(x, y, p_x, p_y)$ — 4 doubles.
- Box: $[-0.5, 0.5]^2$, 4 walls.
- Algorithm: per-particle ballistic motion + analytical computation of
  next wall crossing time + stochastic bounce.
- **Second engine**: `SIM_ENGINE=super` selects the super-event
  closed-form aggregation (O(1) per particle per batch): **96× over
  the exact engine** at 1M-step batches, statistically equivalent,
  exact in the conservative limit (the bit-perfect cycle test passes
  with it too). See [`docs/SUPER_EVENT.md`](docs/SUPER_EVENT.md);
  implementation in `src/physics_super.c`. Default engine remains the
  exact event loop.

## Conventions

- `SIM_SEED` env var for deterministic runs.
- `SIM_ENGINE=super` opts into the super-event engine (default: exact).
- `OMP_NUM_THREADS` controls threads (default 8 for benches).
- Config TOML schema same as the 1D parent.
- Output `X<evolution>.dat` extends the 1D format with extra header
  fields (`chi2y`, `chi2px`, `chi2py`, `isotropy`, `corr`) and a 2D
  joint histogram block at the end (`HXY_BINS × HXY_BINS`).
- Checkpoint binary `graba.dmp` layout:
  `[uint evolution] [N doubles x] [N doubles y] [N doubles px] [N doubles py]`.

## Build / Run

```bash
make             # main-event + symlink
./main config.toml
make test-all    # full test suite
```

## Tests

| Target | What |
|---|---|
| `make test` | Smoke regression (~0.1s) |
| `make verify` | Postconditions relax + periodic (~1.5s) |
| `make periodicity` | Strict 1M-step alfa=0 energy conservation |
| `make isotropy` | $\langle p_x^2 \rangle \approx \langle p_y^2 \rangle$ |
| `make independence` | $\langle p_x p_y \rangle \approx 0$ |
| `make test-all` | All of the above |

## Tools (in `tools/`)

(In progress as of initial commit; see `docs/ALGORITHM.md` for the
full plan.)

| Tool | What |
|---|---|
| `plot.py` | Static plots: marginals, joint heatmaps, radial / angular |
| `animate.py` | Animated GIFs of evolution |
| `tui.py` | Live multi-panel terminal monitor (plotext) |
| `bench.py` | SQLite-recorded benchmarks |

## Things specifically NOT supported here

- CUDA backend. Event-driven matches CPU model better.
- Inter-particle interactions. The model assumes independent particles.
- Bit-equivalence with the parent 1D project. The 2D model is a
  different physics (4 walls, 2 momentum components); only the
  algorithmic technique is shared.

## Things NOT to do

- Don't add CUDA support — it underperforms CPU here for the same
  reasons as in the 1D event-driven implementation.
- Don't change the per-component physics model without updating
  `docs/ALGORITHM.md` §2 and re-running all postconditions.
- Don't break compatibility with the `tools/plot.py` reader without
  updating both ends in the same commit.
