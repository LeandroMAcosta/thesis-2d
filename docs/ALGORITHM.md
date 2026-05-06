# Algoritmo: gas 2D event-driven

Documentación técnica completa de la implementación. Adapta el
documento del repo padre (`event-driven/docs/ALGORITHM.md` 1D) al
caso 2D, manteniendo la estructura para facilitar comparación.

**Última actualización**: 2026-05-07.

---

## Tabla de contenidos

1. [Resumen](#1-resumen)
2. [Modelo físico 2D](#2-modelo-físico-2d)
3. [Dinámica entre eventos](#3-dinámica-entre-eventos)
4. [Cálculo del próximo evento](#4-cálculo-del-próximo-evento)
5. [El bounce: dos casos](#5-el-bounce-dos-casos)
6. [Implementación](#6-implementación)
7. [Postcondiciones físicas](#7-postcondiciones-físicas)
8. [Validación experimental](#8-validación-experimental)
9. [Performance](#9-performance)
10. [Diferencias con el modelo 1D](#10-diferencias-con-el-modelo-1d)
11. [Referencias](#11-referencias)

---

## 1. Resumen

Extensión natural del modelo 1D del repo padre a 2 dimensiones. Las
partículas viven en una caja cuadrada $[-0.5, 0.5]^2$, con 4 paredes
(2 verticales en $x = \pm 0.5$, 2 horizontales en $y = \pm 0.5$).
Entre rebotes el movimiento es balístico bidimensional (línea recta
en el plano). En cada rebote, el modelo aplica una de dos
transformaciones según qué pared se tocó.

El algoritmo es **event-driven**: por cada partícula se calcula
analíticamente el tiempo hasta el próximo rebote ($O(1)$ por evento)
y se salta directo, sin iterar.

---

## 2. Modelo físico 2D

### 2.1 Variables de estado

Cada partícula $i \in \{1, \ldots, N\}$ tiene 4 variables:

| Variable | Descripción | Unidades |
|---|---|---|
| $x_i$ | posición horizontal | m, $\in [-0.5, 0.5]$ |
| $y_i$ | posición vertical | m, $\in [-0.5, 0.5]$ |
| $p_{x,i}$ | momento horizontal | kg·m/s |
| $p_{y,i}$ | momento vertical | kg·m/s |

Velocidad: $\vec{v} = \vec{p}/m = (p_x/m, p_y/m)$.

Energía cinética total del sistema:
$$E_{\text{tot}} = \frac{1}{2m} \sum_{i=1}^N \left(p_{x,i}^2 + p_{y,i}^2\right)$$

### 2.2 La caja y sus paredes

Cuatro paredes, indexadas:

| Wall | Plano | Normal $\hat{n}$ | Componente perpendicular |
|---|---|---|---|
| izquierda | $x = -0.5$ | $(+1, 0)$ | $p_x$ |
| derecha | $x = +0.5$ | $(-1, 0)$ | $p_x$ |
| inferior | $y = -0.5$ | $(0, +1)$ | $p_y$ |
| superior | $y = +0.5$ | $(0, -1)$ | $p_y$ |

### 2.3 Distribuciones esperadas en equilibrio

Tras suficiente tiempo de relajación con $\alpha > 0, \sigma_L > 0$:

- **Posiciones**: uniforme en la caja:
  $f(x, y) = 1$ (densidad por unidad de área).
  Marginales: $f(x) = 1$, $f(y) = 1$ separadamente.

- **Momentos**: Maxwell-Boltzmann 2D, separable:
  $$f(p_x, p_y) = \frac{1}{2\pi \sigma^2} \exp\!\left(-\frac{p_x^2 + p_y^2}{2\sigma^2}\right)$$
  Marginales: $f(p_x), f(p_y)$ son gaussianas $N(0, \sigma^2)$ con
  $\sigma = $ `P_SIGMA`.

- **Coordenadas polares de $\vec{p}$**:
  - Magnitud $|\vec{p}| = \sqrt{p_x^2 + p_y^2}$ sigue una **distribución
    de Rayleigh**:
    $$f(p) = \frac{p}{\sigma^2} \exp(-p^2 / 2\sigma^2)$$
  - Ángulo $\theta = \arctan(p_y / p_x)$ es **uniforme en $[0, 2\pi)$**.

- **Independencia**: $\langle p_x p_y \rangle = 0$, los componentes
  son estadísticamente independientes.

- **Isotropía**: $\langle p_x^2 \rangle = \langle p_y^2 \rangle = \sigma^2$.

Estas son las predicciones que se chequean en las postcondiciones (§7).

---

## 3. Dinámica entre eventos

Mientras todas las partículas estén dentro de la caja **sin tocar
paredes**, el sistema es trivial: cada una se mueve en línea recta
con momento (y energía) constante.

$$\vec{x}_i(t + \Delta t) = \vec{x}_i(t) + \frac{\vec{p}_i}{m} \Delta t$$

Esta es la propiedad clave que habilita event-driven: el estado
futuro entre eventos es **una función cerrada** del estado presente
y el tiempo, no requiere integración numérica.

---

## 4. Cálculo del próximo evento

Dada una partícula con estado $(x, y, p_x, p_y)$ en el instante
$t$, el próximo evento ocurre cuando alcance una de las 4 paredes.
Se calcula el tiempo a cada candidato y se toma el mínimo positivo.

### 4.1 Tiempo a cada pared

**Pared vertical** (cualquier $x = x_w \in \{-0.5, +0.5\}$):

$$t_x = \frac{x_w - x}{v_x} \quad \text{con} \quad v_x = \frac{p_x}{m}$$

Sólo es positivo (físicamente alcanzable en el futuro) si:
- $v_x > 0$ y $x_w = +0.5$, o
- $v_x < 0$ y $x_w = -0.5$.

Análogamente para $y$:

$$t_y = \frac{y_w - y}{v_y}$$

con $v_y = p_y / m$.

### 4.2 Selección del próximo evento

Entre las 4 candidatas, sólo dos son finitas y positivas (las que
respetan el signo de la velocidad). El próximo evento es:

$$t_{\text{event}} = \min(t_x, t_y)$$

Si $t_x < t_y$, la partícula toca **primero** una pared vertical
(side = 0). Si $t_y < t_x$, toca primero una horizontal (side = 1).
La probabilidad de empate exacto es 0 (medida cero); en la práctica
los empates por error de redondeo se rompen consistentemente por
el orden del `<`.

### 4.3 Edge cases

- **$v_x = 0$ o $v_y = 0$**: la pared correspondiente nunca se
  alcanza; tratado como $t = +\infty$.
- **Partícula completamente detenida** ($p_x = p_y = 0$): no se
  mueve, terminamos su evolución.
- **$t < 0$ por drift numérico**: forzar $t \to 0$ para que el
  siguiente loop iter procese el bounce inmediato.

---

## 5. El bounce: dos casos

Sea $\vec{v}^{\text{pre}}$ la velocidad inmediatamente antes del bounce.

### 5.1 Pared vertical (side = 0)

La pared está en $x = \pm 0.5$. La transformación es:

**a) Inversión del momento perpendicular:**
$$p_x \mapsto -p_x$$

**b) Perturbación tangencial de la posición** (rugosidad de la pared):
$$y \mapsto y + \delta y, \quad \delta y = \xi_1 \cos(\xi_2) \sigma_L$$

con $\xi_1, \xi_2$ del par de Box-Muller a partir de dos uniformes.

**c) Perturbación de $|p_x|$** (intercambio energético con el baño térmico):
$$|p_x|^2 \mapsto |p_x|^2 + \alpha (|p_x| - p_{\min})(p_{\max} - |p_x|)(r_3 - 0.5)$$

donde $r_3$ es una uniforme adicional. Después se toma raíz para
recuperar el nuevo $|p_x|$.

**d) Reflexión geométrica** (clamp para $|y| > 0.502$):
$$y \mapsto \begin{cases} 1.004 \operatorname{sgn}(y) - y & \text{si } |y| > 0.502 \\ y & \text{en caso contrario}\end{cases}$$

### 5.2 Pared horizontal (side = 1)

Análoga con los roles de $x \leftrightarrow y$ y $p_x \leftrightarrow p_y$
intercambiados.

### 5.3 Conservación de energía aproximada

En cada bounce, $|p_x|$ (o $|p_y|$) cambia por:

$$\Delta(|p|) = \sqrt{|p|^2 + \alpha(|p| - p_{\min})(p_{\max} - |p|)(r-0.5)} - |p|$$

Para $\alpha$ chico, $\Delta E$ por bounce es pequeño y de signo
aleatorio. Sumado sobre todos los bounces y todas las partículas, la
energía media se conserva con deriva $< 1$ ppm a lo largo de la
corrida (cf. postcondiciones).

Para $\alpha = 0$ (régimen periódico), $|p|$ no cambia → conservación
**bit-for-bit** verificable empíricamente.

---

## 6. Implementación

### 6.1 Estructura del kernel principal

```c
void evolve_one_particle(x, y, px, py, seed, T_target, ...) {
    t = 0
    while (true) {
        if (px == 0 && py == 0) break;       // partícula detenida

        vx = px / m;  vy = py / m;
        compute t_x (= ∞ if vx == 0)
        compute t_y (= ∞ if vy == 0)
        clamp negative drift to 0

        t_event = min(t_x, t_y)
        side = (t_x < t_y) ? 0 : 1

        if (t + t_event >= T_target) {
            advance ballistically to T_target
            break
        }

        t += t_event
        x += vx * t_event
        y += vy * t_event
        snap perpendicular coord exactly to ±0.5

        apply_bounce(side, sign of perpendicular velocity)

        if (event_count > EVENTS_MAX) break  // panic guard
    }
}
```

Loop principal sobre partículas, paralelizado con OpenMP:

```c
#pragma omp parallel for schedule(dynamic, 1024)
for (int i = 0; i < N; i++) {
    seed = rng_seed_mix(base, i)
    evolve_one_particle(s->x[i], s->y[i], s->px[i], s->py[i], &seed, ...)
}
```

### 6.2 OpenMP scheduling

`schedule(dynamic, 1024)` por la misma razón que el repo padre 1D:
las partículas tardan tiempos variables (una rápida hace muchos
bounces, una lenta hace pocos), y dynamic balancea automáticamente.

Empíricamente, el paso de static a dynamic mejora el wall en
~30-50% para corridas medianas con 8 cores.

### 6.3 RNG y reproducibilidad

Cada partícula usa una secuencia xorshift32 sembrada con
`rng_seed_mix(base, particle_index)`. `base` se deriva del `SIM_SEED`
y del índice del batch. Esto garantiza:

- **Reproducibilidad bit-for-bit** entre corridas con el mismo
  `SIM_SEED`.
- **Independencia entre threads**: como cada partícula tiene su
  propia secuencia, el orden en que los threads procesen las
  partículas no afecta el resultado.

### 6.4 Compatibilidad de I/O con visualización

El formato `.dat` extiende el del repo padre 1D con campos extra:

- Header: `chi2x, chi2y, chi2px, chi2py, isotropy, corr, Et, N, BINS, HXY_BINS`.
- Bloque marginales: 2*BINS + 4 filas con columnas `x_val, hx, y_val, hy, p_val, gpx, gpy`.
- Bloque joint: `HXY_BINS × HXY_BINS` matriz de counts del histograma 2D
  de posiciones (resolución coarse para mantener tamaño manejable).

El binario `graba.dmp` también se extiende:
`[uint evolution] [N x] [N y] [N px] [N py]` (4 arrays double en lugar de 2).

---

## 7. Postcondiciones físicas

### 7.1 Postcondiciones generales (todos los regímenes)

| # | Condición | Cota |
|---|---|---|
| 1 | Conservación de partículas: $\sum h(x_i) = N$, idem para los 4 marginales | exacto |
| 2 | Todos los counts no-negativos | exacto |

### 7.2 Postcondiciones del régimen de relajación ($\alpha > 0, \sigma_L > 0$)

| # | Condición | Cota |
|---|---|---|
| 3 | $\chi^2_x < 5$ | (vs uniforme) |
| 4 | $\chi^2_y < 5$ | (vs uniforme) |
| 5 | $\chi^2_{p_x} < 5$ | (vs Maxwell-Boltzmann) |
| 6 | $\chi^2_{p_y} < 5$ | (vs MB) |
| 7 | Isotropía: $\langle p_x^2 \rangle / \langle p_y^2 \rangle \in [0.9, 1.1]$ | |
| 8 | Independencia: $\rho(p_x, p_y) = \langle p_x p_y \rangle / \sqrt{\langle p_x^2 \rangle \langle p_y^2 \rangle} \in [-0.05, +0.05]$ | |
| 9 | Deriva de energía: $|\Delta E_{\text{tot}}/E_{\text{tot}}| < 100$ ppm | |

### 7.3 Postcondiciones del régimen periódico ($\alpha = 0, \sigma_L = 0$)

| # | Condición | Cota |
|---|---|---|
| 10 | Energía conservada bit-for-bit en todos los snapshots | exacto (o spread $< 10^{-12}$) |

Todas estas postcondiciones están implementadas en
`tests/check_postconditions.py` y pasan en el commit inicial (ver §8).

---

## 8. Validación experimental

Resultados a fecha del primer commit:

### Régimen de relajación ($N = 2^{18}$, 1M pasos, $\alpha = \sigma_L = 10^{-4}$)

| Postcondición | Cota | Observado | Status |
|---|---|---|---|
| $\chi^2_x$ | < 5 | **1.047** | ✅ |
| $\chi^2_y$ | < 5 | **0.990** | ✅ |
| $\chi^2_{p_x}$ | < 5 | **1.034** | ✅ |
| $\chi^2_{p_y}$ | < 5 | **0.925** | ✅ |
| Isotropía | $\in [0.9, 1.1]$ | **0.9944** | ✅ |
| Correlación | $\in [-0.05, +0.05]$ | **+2.36e-3** | ✅ |
| Deriva de energía | < 100 ppm | **+0.02 ppm** | ✅ |
| Conservación de $N$ | exacto | $N = 262\,144$ en 11 snapshots | ✅ |

### Régimen periódico ($\alpha = 0, \sigma_L = 0$)

| Postcondición | Status |
|---|---|
| Conservación de partículas | ✅ |
| Energía bit-idéntica en 11 snapshots ($E_t = 1.0867 \times 10^{-15}$ J) | ✅ |

---

## 9. Performance

A medir formalmente. Preview empírico:

- **Verify completo (relax + periodic, $N = 2^{18}$, 2 × 1M pasos)**: ~1.5 s wall.
- **Smoke regression** ($N = 4096$, 300 pasos): ~0.05 s wall.
- **Periodicity strict** ($N = 2^{18}$, 1M pasos): ~0.6 s wall.

Comparación cuantitativa con el repo 1D y escalabilidad por $N$:
ver `docs/PERFORMANCE.md` (TBD).

---

## 10. Diferencias con el modelo 1D

| Aspecto | 1D (repo padre) | 2D (este repo) |
|---|---|---|
| State per particle | $(x, p)$ | $(x, y, p_x, p_y)$ |
| Walls | 2 ($x = \pm 0.5$) | 4 (vertical + horizontal) |
| Operating cost per event | 3 RNG, 1 sqrt + log + cos, 1 alfa-loop | mismo (+ comparación de 2 candidatos para el min) |
| Distribución de equilibrio | uniforme en $x$, gaussiana en $p$ | uniforme en $\vec{x}$, MB 2D en $\vec{p}$ |
| Postcondiciones extra | — | isotropía, independencia, distribución radial / angular |
| Output formato | hx, gpx | hx, hy, gpx, gpy + h_xy joint |

---

## 11. Referencias

Comparte la base bibliográfica del documento padre. Para citar este
trabajo específicamente:

- Acosta, L. *Tesis de licenciatura, FAMAF, UNC* (2026), capítulo
  sobre extensión 2D de event-driven.
- Repo 1D: `event-driven/docs/ALGORITHM.md` en
  https://github.com/LeandroMAcosta/thesis.

Para Maxwell-Boltzmann 2D y sus marginales (distribución de Rayleigh
para $|\vec{p}|$, uniforme para $\theta$):

- Reichl. *A Modern Course in Statistical Physics*. Wiley, 2016. §1.2-1.4.
- Pathria, Beale. *Statistical Mechanics*, Academic Press, 4ta ed., 2021.
