# thesis-2d — gas 2D event-driven

Implementación CPU + OpenMP de la simulación GCAS extendida a **2 dimensiones**,
usando **integración por eventos analítica** (la misma técnica que el repo
[thesis](https://github.com/LeandroMAcosta/thesis) en su variante
`event-driven/`).

## Resumen

- $N$ partículas en una caja $[-0.5, 0.5]^2$ rebotando contra 4 paredes.
- Cada partícula tiene estado $(x, y, p_x, p_y)$.
- Entre rebotes, movimiento balístico en línea recta.
- En cada rebote, perturbación estocástica de la posición tangencial y de la
  magnitud del momento perpendicular (mismo modelo del paper original 1D
  extendido a 2D de la manera natural).
- Integración **event-driven**: cada partícula salta de un rebote al siguiente
  con cómputo $O(1)$ por evento.

## Build / Run

```bash
make            # main-event + symlink ./main → main-event
./main tests/config.regression.toml
```

## Tests

```bash
make test            # smoke (4096 particles, 300 steps, ~0.1s)
make verify          # postconditions relax + periodic (~1.5s)
make periodicity     # alfa=0, energy bit-conservation
make isotropy        # <p_x²> ≈ <p_y²>
make independence    # <p_x·p_y> ≈ 0
make test-all        # todo en orden
```

Todas las postcondiciones documentadas en [`docs/ALGORITHM.md` §7](docs/ALGORITHM.md).

## Resultados de validación (a tag inicial)

A $N = 2^{18}$, 1M pasos equivalentes (alfa=1e-4):

| Postcondición | Cota | Observado |
|---|---|---|
| $\chi^2_x$ | < 5 | **1.047** |
| $\chi^2_y$ | < 5 | **0.990** |
| $\chi^2_{p_x}$ | < 5 | **1.034** |
| $\chi^2_{p_y}$ | < 5 | **0.925** |
| Isotropía $\langle p_x^2 \rangle / \langle p_y^2 \rangle$ | $\in [0.9, 1.1]$ | **0.9944** |
| Correlación $\rho(p_x, p_y)$ | $|\cdot| < 0.05$ | **+2.4e-3** |
| Deriva de energía | < 100 ppm | **+0.02 ppm** |
| Energía conservada bit-for-bit (alfa=0) | exacto | $E = 1.087 \times 10^{-15}$ J en los 11 snapshots ✓ |

## Performance

A $N = 2^{21}$, 500 k pasos equivalentes — pendiente bench formal (ver
`docs/PERFORMANCE.md` cuando se complete).

## Estructura

```
thesis-2d/
├── README.md              # este archivo
├── CLAUDE.md              # instrucciones para LLMs
├── docs/
│   ├── ALGORITHM.md       # documentación técnica completa
│   ├── PERFORMANCE.md     # bench y comparación con el repo 1D
│   └── figures/           # plots y animaciones del informe
├── Makefile
├── include/               # headers C
├── src/                   # implementación
├── third_party/tomlc99    # parser TOML vendored
├── tests/                 # configs, scripts, checkers Python
└── tools/                 # plot, animate, tui, bench
```

## Documentación

- [`docs/ALGORITHM.md`](docs/ALGORITHM.md) — derivación matemática, equivalencia
  con stepping, casos borde, validación.
- [`CLAUDE.md`](CLAUDE.md) — contexto para sesiones futuras de Claude Code.

## Repo padre

- [`thesis`](https://github.com/LeandroMAcosta/thesis) — implementación 1D
  original con backends OMP-stepping, CUDA (FP64, df64, df-trig) y
  event-driven CPU (la base de esta extensión 2D).

## Licencia / autoría

- Tesis de licenciatura en Ciencias de la Computación, FAMAF, UNC.
- Director: Nicolás Wolovick. Co-director: Gustavo Castellano.
- Autor: Leandro Acosta, 2026.
