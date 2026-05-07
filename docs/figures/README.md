# Galería de figuras

Salida del toolkit `tools/plot` y `tools/animate`. Tres regímenes
ilustrados, en orden de "fuerza física" del retorno a la condición
inicial.

---

## Los regímenes y sus parámetros

| Régimen | $\alpha$ | $\sigma_L$ | `init_uniform_p` | Vuelve a la esquina? |
|---|---|---|---|---|
| **Cycle** (§2) | 0 | 0 | **sí**, calibrado | **sí, exactamente** (L1 distance = 0) |
| **Periodic** (§4) | 0 | 0 | no (gaussian) | sólo en parte (Poincaré recurrence astronómica) |
| **Relax** (§3) | 1e-4 | 1e-4 | no | no — se equilibra a Maxwell-Boltzmann |

Las dos demos sin retorno usan `--palindrome` en `animate.py` para
forzar visualmente el ida-y-vuelta. La demo cycle muestra el retorno
**real**.

> 📐 **El truco del cycle**: si todas las partículas tienen $|p_x| = |p_y| = p_0$
> con $p_0$ elegido tal que $T_{\text{periodo}} = 2\,m\,L / p_0$ divide
> exactamente $N_{\text{steps}} \cdot dt$, todas vuelven a su estado
> inicial al mismo tiempo. Concretamente, con $p_0 \approx 5.2500\times 10^{-24}$
> kg·m/s y los parámetros estándar, encajan exactamente **175 períodos** en
> 1 000 000 pasos. Este es el **único caso donde la simulación es
> verdaderamente cíclica** sin trampa visual.

---

## 1. Estado inicial — todas en una esquina

65 536 partículas distribuidas uniformemente en $[0, 0.5]^2$ (cuadrante
superior-derecho). Joint $h(x, y)$:

![Joint initial corner](joint_initial_corner.png)

Marginales iniciales — $h(x), h(y)$ concentrados en $[0, 0.5]$;
$g(p_x), g(p_y)$ ya gaussianas:

![Marginales iniciales](marginals_initial.png)

Scatter del estado inicial:

![Scatter corner](scatter_corner_initial.png)

```bash
./tools/plot joint X0000000.dat -o joint_initial_corner.png
./tools/plot marginals X0000000.dat -o marginals_initial.png
./tools/plot scatter --snapshot X0000000.dat --dump graba.dmp.0000150 \
    --n-sub 10000 -o scatter_corner_initial.png
```

---

## 2. **Ciclo periódico real** — todas vuelven exactamente

**Configuración**: $N = 65\,536$, $\alpha = 0$, $\sigma_L = 0$,
**`init_uniform_p_magnitude = true`** con $p_0 = 5.25\times10^{-24}$
calibrado para 175 períodos exactos en 1 000 000 pasos.

**Cadencia exponencial simétrica** (snapshots concentrados en los
extremos para ver bien la salida y el retorno):

```python
# evolutions: 0, 5, 11, 26, 59, 134, 305, 695, 1581, 3598, 8189, 18638,
#             42417, 96535, 219699, 500000, 780301, 903465, 957583,
#             981362, 991811, 996402, 998419, 999305, 999695, 999866,
#             999941, 999974, 999989, 999995, 1000000
steps = [5, 6, 15, 33, 75, 171, 390, 886, 2017, 4591, 10449, 23779,
         54118, 123164, 280301, 280301, 123164, 54118, 23779, 10449,
         4591, 2017, 886, 390, 171, 75, 33, 15, 6, 5]
```

> ✅ **Verificación numérica**: L1 distance per particle entre
> snapshot 0 y snapshot 1M = **0.000000** exacto. Las 65 536 partículas
> están en las **mismas celdas** del histograma 64×64 al inicio y al
> final. El joint $h(x,y)$ es **bit-idéntico**.

### Animación: la difusión y el retorno

| | |
|---|---|
| ![Cycle joint](anim_cycle_joint.gif) | ![Cycle scatter](anim_cycle_scatter.gif) |

Vas a ver: salida lenta y suave de la esquina superior-derecha, llegada
al equilibrio aparente en el medio (todas las posiciones ocupadas), y
**retorno smooth a la esquina** al final. Sin trampa: es el sistema real
ejecutándose.

### Marginales — vuelven exactas

![Cycle marginales](anim_cycle_marginals.gif)

$h(x), h(y)$ comienzan concentrados en $[0, 0.5]$, se difunden, vuelven.
$g(p_x), g(p_y)$ son delta-functions discretas (sólo dos valores: $\pm p_0$),
permanecen invariantes — sólo cambian de signo en cada bounce.

### Momentos en el ciclo

![Cycle pscatter](anim_cycle_pscatter.gif)

Los momentos toman exactamente 4 valores: $(\pm p_0, \pm p_0)$. La nube
es **4 puntos** que rotan entre los cuatro cuadrantes según los bounces.

### Por qué funciona

- Cada componente $p_x, p_y$ tiene el mismo $|p|$ → mismo período de bounce.
- $T_{\text{period}} = 2 \cdot L \cdot m / p_0$.
- Eligiendo $p_0$ tal que $N_{\text{steps}} \cdot dt = k \cdot T_{\text{period}}$
  para $k$ entero, todas las partículas, sin importar su posición inicial
  dentro del cuadrante, vuelven a su estado inicial al mismo tiempo.
- Conservación bit-exacta de la energía (verificable: $E_t$ idéntico en los 31
  snapshots).
- La dinámica es **completamente reversible** y **periódica con período T_total**.

---

## 3. Régimen de relajación — difusión hasta equilibrio (con `--palindrome`)

**Configuración**: $N = 65\,536$, $\alpha = 10^{-4}$, $\sigma_L = 10^{-4}$,
30 batches con cadencia campana ($\cos^2$).

> 🔁 **Animaciones con `--palindrome`** en este régimen: las partículas
> NO vuelven realmente a la esquina (la disipación rompe la
> reversibilidad y, además, las trayectorias gaussianas son
> incommensurables). El GIF reproduce los frames de ida y los mismos
> frames al revés, como **ilustración** del concepto. Para retorno
> real, ver §2.

| Joint | Scatter |
|---|---|
| ![Relax joint](anim_joint_relax.gif) | ![Relax scatter](anim_scatter_relax.gif) |

Marginales convergiendo a uniforme (espacial) y MB (momento):

![Relax marginales](anim_marginals_relax.gif)

Momentos relajando — la nube se reorganiza estadísticamente:

![Relax pscatter](anim_pscatter_relax.gif)

### Estado al final del relax

![Dashboard relax](dashboard_relax.png)

Validación: $\chi^2_x = 1.05$, $\chi^2_y = 1.04$, isotropy = 0.99,
correlación ≈ 0, deriva de energía pocos ppm.

### Distribuciones $|\vec p|$ y $\theta$ al final

| Radial vs Rayleigh | Angular uniforme |
|---|---|
| ![Radial relax](radial_relax_final.png) | ![Angular relax](angular_relax_final.png) |

### Deriva de energía relax

![Energy drift relax](energy_drift_relax.png)

---

## 4. Régimen periódico — gaussian momenta, $\alpha=\sigma_L=0$ (con `--palindrome`)

**Configuración**: misma cadencia campana, pero **momentos gaussianos**
(no uniformes). Sin disipación, energía conservada bit-for-bit, pero
los momentos tienen 1000 magnitudes distintas (cuantizadas a bin
centers) → **Poincaré recurrence astronómica**, no vuelve a la
esquina en 1M pasos.

| Joint | Scatter |
|---|---|
| ![Periodic joint](anim_joint_periodic.gif) | ![Periodic scatter](anim_scatter_periodic.gif) |

Marginales periódicas — $h(x), h(y)$ phase-mixean, $g(p_x), g(p_y)$
**permanecen idénticas**:

![Periodic marginales](anim_marginals_periodic.gif)

Momentos en régimen periódico — la nube **NO se redistribuye**, cada
partícula conserva su $|p|$:

![Periodic pscatter](anim_pscatter_periodic.gif)

Energía conservada bit-for-bit (línea perfectamente plana):

![Energy drift periodic](energy_drift_periodic.png)

Dashboard final:

![Dashboard periodic](dashboard_periodic.png)

---

## 5. Comparación de los tres regímenes

| | Cycle (§2) | Relax (§3) | Periodic (§4) |
|---|---|---|---|
| Init de momentos | uniforme $\pm p_0$ | gaussian | gaussian |
| Disipación / ruido pared | no | sí | no |
| Posiciones | vuelven exacto | Maxwell-Boltzmann | phase-mixing |
| Momentos | 4 valores discretos | redistribuyen | invariantes |
| Energía | bit-exacta | drift ppm | bit-exacta |
| Retorno a la esquina | **sí, real** | ilustrativo (palíndromo) | ilustrativo (palíndromo) |

La distinción visual más clara está entre las animaciones `pscatter`:

- **Cycle**: 4 puntos, rotación discreta.
- **Periodic**: nube fija (cada partícula conserva $|p|$).
- **Relax**: nube reorganizándose (estocástico).

---

## Comandos de reproducción

```bash
# Demo cycle (retorno exacto)
make
mkdir -p /tmp/cycle && cp <config con uniform_p_magnitude> /tmp/cycle/config.toml
SIM_SEED=42 OMP_NUM_THREADS=8 ./main-event config.toml

# Plots
./tools/plot dashboard X1000000.dat --dump graba.dmp.1000000 -o dashboard.png

# Anims sin palíndromo (cycle: el retorno es real)
./tools/animate --duration 0.18 joint     /tmp/cycle -o anim_cycle_joint.gif
./tools/animate --duration 0.18 marginals /tmp/cycle -o anim_cycle_marginals.gif
./tools/animate --duration 0.18 scatter   /tmp/cycle -o anim_cycle_scatter.gif
./tools/animate --duration 0.18 pscatter  /tmp/cycle -o anim_cycle_pscatter.gif

# Anims con palíndromo (regímenes que no vuelven realmente)
./tools/animate --duration 0.12 --palindrome joint /tmp/demo_relax -o anim_joint_relax.gif
# ... (idem para los otros tres modes en los dos regímenes)
```
