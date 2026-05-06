# Performance — gas 2D event-driven

Benchmarks y comparación con el repo padre 1D.

**Última actualización**: 2026-05-07.

## Resumen

A $N = 2^{21}$ partículas, 500 k pasos equivalentes (≈ 0.22 s simulados):

| Implementación | Wall (s) | Mst/s | vs OMP-1D-stepping |
|----------------|---------:|------:|-------------------:|
| 1D OMP stepping (paper) | 125.26 | 0.004 | 1.00× |
| 1D CUDA tuned (df-trig + round 2) | 33.48 | 0.015 | 3.74× |
| **1D event-driven CPU OMP** | 1.01 | 0.494 | **124×** |
| **2D event-driven CPU OMP** | 2.85 | 0.176 | **44× equivalente** |

El 2D es ~2.8× más lento que el 1D event-driven, **lo esperado**: dobla la
memoria por partícula (4 doubles vs 2), aproximadamente dobla el número de
bounces (porque la partícula se mueve por dos dimensiones), y mantiene
aproximadamente el mismo costo por bounce. Aun así supera al CUDA-1D
optimizado del repo padre por **un factor 12×**, en CPU pura.

## Escalabilidad por N

500 k pasos equivalentes, 8 threads OMP, $\sigma_L = \alpha = 10^{-4}$:

| $N$ | Wall (s) | Wall por partícula (ns) |
|-----|---------:|------------------------:|
| $2^{16}$ = 65 536 | 0.11 | 1.68 |
| $2^{18}$ = 262 144 | 0.38 | 1.45 |
| $2^{20}$ = 1 048 576 | 1.49 | 1.42 |
| $2^{21}$ = 2 097 152 | 2.85 | 1.36 |
| $2^{22}$ = 4 194 304 | 5.98 | 1.43 |

**Escalabilidad lineal en $N$** (esperable para algoritmo embarazosamente
paralelo). Los ~1.4 ns/partícula son constantes; la pequeña variación a
$N$ chico viene de overhead fijo (init, alloc, finalize) y del
underutilization de OMP cuando hay menos partículas que threads × chunk
size.

## Setup experimental

- **Hardware**: AMD Ryzen 7 7800X3D (8 cores Zen 4 con AVX2). RAM DDR5.
- **Software**: Linux 6.6 (WSL2 Ubuntu sobre Windows 11). gcc 13 con
  `-O3 -march=native -ffast-math -fopenmp`.
- **Configuración canónica**: `tests/config.bench.big.toml`
  ($N = 2^{21}$, 5 batches × 100k steps = 500k steps, $\alpha = 10^{-4}$).
- **Mediciones**: 1 corrida por punto reportado. Variancia entre corridas
  es < 5 % en este host (cf. el repo padre 1D).
- **Determinismo**: `SIM_SEED = 42`.

## Cómo reproducir

```bash
make
./tools/setup.sh             # one-time venv setup
./tools/bench run tests/config.bench.big.toml --notes "..."
./tools/bench compare
./tools/bench list --limit 10
```

## Por qué no CUDA

Mismo argumento que en el repo padre (`event-driven/docs/ALGORITHM.md` §9):

- Cada partícula tarda un tiempo variable (depende del número de bounces,
  que es estocástico).
- En GPU, las partículas viven en warps de 32 hilos lockstep; la varianza
  intra-warp en el número de eventos desperdicia el warp.
- En CPU OMP, schedule `dynamic, 1024` balancea automáticamente.

Estimación: una versión CUDA event-driven tardaría 5-10 s en lugar de 3 s
en este hardware, perdiendo el speedup sobre 1D-stepping-CUDA-tuned. La
GPU es ineficiente para cargas con alta varianza intra-warp.

## Trabajo futuro

- Bench en cluster (Serafín / Nabu) con CPUs más grandes y sin WSL overhead.
- Variante CUDA con "bucket scheduling" (agrupar partículas por velocidad
  similar para reducir varianza intra-warp).
- Comparación con `make verify` runtime entre 1D y 2D para tener una
  segunda métrica.
- Profiling con `perf` para ver si el cuello sigue siendo el RNG o ya las
  trig de Box-Muller.
