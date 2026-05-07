# Galería de figuras

Salida del toolkit `tools/plot` y `tools/animate`. Tres regímenes
ilustrados.

Todas las animaciones tienen **evolución monótona de 0 a 1 000 000
pasos** con **cadencia tipo sigmoide**: pocos pasos al principio,
muchos en el medio, pocos al final. Concretamente:

```
evol = 0 → 5 → 11 → 26 → 59 → 134 → 305 → 695 → 1581 → 3598 →
       8189 → 18638 → 42417 → 96535 → 219699 → 500000 →
       780301 → 903465 → 957583 → 981362 → 991811 → 996402 →
       998419 → 999305 → 999695 → 999866 → 999941 → 999974 →
       999989 → 999995 → 1000000
```

(31 snapshots, 30 batches; los pasos `steps[i]` son simétricos en
torno al medio, formando una "campana" que produce densidad de
snapshots concentrada en los dos extremos.)

---

## Los regímenes

| Régimen | $\alpha$ | $\sigma_L$ | `init_uniform_p` | En evol = 1M |
|---|---|---|---|---|
| **Cycle** (§2) | 0 | 0 | **sí**, calibrado | **vuelve EXACTO a la esquina** (L1 distance = 0) |
| **Periodic** (§4) | 0 | 0 | no (gaussian) | sigue disperso (Poincaré recurrence astronómica) |
| **Relax** (§3) | 1e-4 | 1e-4 | no | en equilibrio (Maxwell-Boltzmann + uniforme espacial) |

> 📐 **El truco del cycle**: si todas las partículas tienen $|p_x| = |p_y| = p_0$
> con $p_0$ elegido tal que $T_{\text{periodo}} = 2\,m\,L / p_0$ divide
> exactamente $N_{\text{steps}} \cdot dt$, todas vuelven a su estado
> inicial al mismo tiempo. Concretamente, con $p_0 \approx 5.2500\times 10^{-24}$
> kg·m/s y los parámetros estándar, encajan exactamente **175 períodos** en
> 1 000 000 pasos. Este es el **único caso donde la simulación es
> verdaderamente cíclica**.

---

## 1. Estado inicial — todas en una esquina

65 536 partículas distribuidas uniformemente en $[0, 0.5]^2$ (cuadrante
superior-derecho).

![Joint initial corner](joint_initial_corner.png)

```bash
./tools/plot joint X0000000.dat -o joint_initial_corner.png
```

---

## 2. **Ciclo periódico exacto** — todas vuelven a la esquina en evol = 1M

**Configuración**: $N = 65\,536$, $\alpha = 0$, $\sigma_L = 0$,
**`init_uniform_p_magnitude = true`** con $p_0 = 5.25\times10^{-24}$
calibrado para 175 períodos exactos en 1 000 000 pasos.

> ✅ **Verificación numérica**: L1 distance per particle entre
> snapshot 0 y snapshot 1M = **0.000000** exacto. Las 65 536 partículas
> están en las **mismas celdas** del histograma 64×64 al inicio y al
> final.

### Heatmap $h(x, y)$ — sale de la esquina y vuelve

![Cycle joint](anim_cycle_joint.gif)

### Scatter de posiciones — la difusión y el retorno

![Cycle scatter](anim_cycle_scatter.gif)

### Marginales — vuelven exactas a su forma inicial

![Cycle marginales](anim_cycle_marginals.gif)

### Momentos — 4 puntos discretos rotando entre los cuadrantes

![Cycle pscatter](anim_cycle_pscatter.gif)

Los momentos toman exactamente 4 valores: $(\pm p_0, \pm p_0)$, así que
"la nube" son 4 puntos discretos. Cada bounce los rota entre los
cuatro cuadrantes del plano $(p_x, p_y)$.

---

## 3. Régimen de relajación — difusión hasta equilibrio

**Configuración**: $N = 65\,536$, $\alpha = 10^{-4}$, $\sigma_L = 10^{-4}$,
mismos pasos sigmoidales. **Las partículas no vuelven a la esquina**:
la disipación rompe la reversibilidad y el sistema se equilibra.

| Joint | Scatter |
|---|---|
| ![Relax joint](anim_joint_relax.gif) | ![Relax scatter](anim_scatter_relax.gif) |

Marginales convergiendo:

![Relax marginales](anim_marginals_relax.gif)

Momentos relajando — la nube se reorganiza estocásticamente:

![Relax pscatter](anim_pscatter_relax.gif)

### Estado al final del relax

![Dashboard relax](dashboard_relax.png)

Validación: $\chi^2_x \approx 1$, $\chi^2_y \approx 1$, isotropy ≈ 1,
correlación ≈ 0, deriva de energía pocos ppm.

### Distribuciones $|\vec p|$ y $\theta$ al final

| Radial vs Rayleigh | Angular uniforme |
|---|---|
| ![Radial relax](radial_relax_final.png) | ![Angular relax](angular_relax_final.png) |

### Deriva de energía relax

![Energy drift relax](energy_drift_relax.png)

---

## 4. Régimen periódico (gaussian momenta) — sin disipación pero sin retorno

**Configuración**: misma cadencia sigmoidal, $\alpha = 0$, $\sigma_L = 0$,
**momentos gaussianos** (no uniformes). Sin disipación → energía
conservada bit-for-bit. Pero los momentos tienen ~1000 magnitudes
distintas (cuantizadas a bin centers); la **Poincaré recurrence** es
astronómica → no vuelve a la esquina en 1M pasos, las partículas
quedan repartidas.

| Joint | Scatter |
|---|---|
| ![Periodic joint](anim_joint_periodic.gif) | ![Periodic scatter](anim_scatter_periodic.gif) |

Marginales — $h(x), h(y)$ phase-mixean; $g(p_x), g(p_y)$
**permanecen idénticas snapshot tras snapshot**:

![Periodic marginales](anim_marginals_periodic.gif)

Momentos — la nube **NO se redistribuye**, cada partícula conserva
$|p|$ (sólo cambia de signo en los bounces):

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
| Posiciones en evol = 1M | **vuelven exacto** | uniforme (equilibrio) | uniforme (phase-mixing) |
| Momentos | 4 valores discretos | redistribuyen | invariantes |
| Energía | bit-exacta | drift ppm | bit-exacta |

La distinción visual más clara está entre las animaciones `pscatter`:

- **Cycle**: 4 puntos, rotación discreta.
- **Periodic**: nube fija (cada partícula conserva $|p|$).
- **Relax**: nube reorganizándose.

---

## Comandos de reproducción

```bash
make
mkdir -p /tmp/cycle && cp <config con uniform_p_magnitude> /tmp/cycle/config.toml
SIM_SEED=42 OMP_NUM_THREADS=8 ./main-event config.toml

# Plots estáticos
./tools/plot dashboard X1000000.dat --dump graba.dmp.1000000 -o dashboard.png

# Animaciones — cadencia sigmoidal monótona 0→1M (no usa --palindrome)
./tools/animate --duration 0.18 joint     /tmp/cycle -o anim_cycle_joint.gif
./tools/animate --duration 0.18 marginals /tmp/cycle -o anim_cycle_marginals.gif
./tools/animate --duration 0.18 scatter   /tmp/cycle -o anim_cycle_scatter.gif
./tools/animate --duration 0.18 pscatter  /tmp/cycle -o anim_cycle_pscatter.gif
```

> El flag `--palindrome` en `animate.py` existe pero **no se usa en
> esta galería**: las animaciones son monótonamente crecientes en
> evol. El palíndromo es opcional para casos donde quieras visualizar
> la reversibilidad temporal del sistema; las demos cycle ya muestran
> el retorno físico real, sin necesidad de truco.
