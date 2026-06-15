# Super-event en 2D: agregación O(1) por batch

La derivación completa, la lección del RNG y el harness de validación
viven en el repo 1D padre:
[`thesis/event-driven/docs/SUPER_EVENT.md`](../../thesis/event-driven/docs/SUPER_EVENT.md).
Acá: lo específico de 2D y los números de este repo.

**Selección**: `SIM_ENGINE=super` (default: engine event-driven exacto).

```bash
SIM_SEED=42 SIM_ENGINE=super ./main-event config.toml
```

## Qué agrega 2D sobre la versión 1D

Cada eje es un problema 1D independiente salvo los kicks tangenciales:
una pared vertical kickea `y`, una horizontal kickea `x` (la coordenada
que rebota se snapea exacto, sin kick propio). A nivel agregado:

- varianza de kicks del eje x = `n_y · σ_L²` (y viceversa);
- el paseo de p² de cada eje usa sólo sus propios rebotes.

**Híbrido descompuesto por eje** (`src/physics_super.c`): cada eje
decide por separado — super cerrado si pasa las guardas TCL (n ≥ 16,
sin riesgo de clamp de p², D≈const), mini-loop exacto 1D si no (1
uniforme por rebote, sin Box-Muller) con el kick agregado aplicado al
final. Sólo si ambos ejes son lentos (~0,1% a batches de 1M pasos) cae
al loop exacto acoplado, con el stream de RNG intacto.

## Validación

- `SIM_ENGINE=super make test-all` ✓ — verify (relax + periodic),
  periodicity con **cycle bit-perfect** (0/1004 bins en los 4
  marginales, 0/4096 celdas en h(x,y), energía bit-idéntica),
  isotropía, independencia.
- A/B contra el engine exacto (5×1M pasos, N=2²¹, seed 42): todos los
  chi² ≈ 1 en ambos, isotropía idéntica a 4 decimales, drift de
  energía −0,08 ppm (super) vs −3,4 ppm (exacto) — ambos dentro del
  ruido.

En el límite conservativo (alfa=0, sigma_l=0) el camino super es
exacto (fold cerrado, |p| copiado bit a bit) — por eso el cycle test
pasa bit-perfect.

## Benchmarks (Ryzen 7 7800X3D, 8 cores, N=2²¹, seed 42)

| Config | event (exacto) | super | speedup |
|---|---:|---:|---:|
| 5×100k pasos | 2.66 s | 0.35 s | 7.6× |
| 5×1M pasos | 24.93 s | 0.26 s | **96×** |

El costo del super es O(1) por partícula por batch: el speedup crece
linealmente con los pasos por batch, y a escala de producción el wall
time queda dominado por init + histogramas + I/O.
