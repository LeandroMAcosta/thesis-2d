# Galería de figuras

Salida del toolkit `tools/plot` y `tools/animate`. Todas las figuras se
pueden regenerar con los comandos en cada sección.

## Régimen de relajación ($\alpha = 10^{-4}, \sigma_L = 10^{-4}$)

Corrida: $N = 2^{18} = 262\,144$, 1M pasos en 10 batches. Las
partículas arrancan distribuidas uniformemente en posición y con
momentos gaussianos de las condiciones iniciales, y evolucionan hasta
equilibrio térmico.

### Dashboard combinado

Vista global en un solo PNG: marginales de posición y momento,
heatmap del joint $h(x,y)$, joint $g(p_x, p_y)$ desde partículas,
distribución radial $|\vec{p}|$, distribución angular $\theta$,
deriva de energía a lo largo de los 11 snapshots.

![Dashboard relax](dashboard_relax.png)

Comando: `./tools/plot dashboard X1000000.dat --dump graba.dmp.1000000`

### Marginales: inicial vs final

Inicial — distribuciones todavía con sello de las condiciones iniciales:

![Marginales inicial](marginals_relax_initial.png)

Final — convergencia a uniforme en $x, y$ y a Maxwell-Boltzmann en $p_x, p_y$:

![Marginales final](marginals_relax_final.png)

### Joint $h(x,y)$ y $g(p_x, p_y)$

![Joint](joint_relax_final.png)

### Distribución radial $|\vec{p}|$ vs Rayleigh

![Radial](radial_relax_final.png)

### Distribución angular $\theta = \arctan(p_y / p_x)$

![Angular](angular_relax_final.png)

### Scatter: posiciones y momentos (subsampleados)

![Scatter](scatter_relax_final.png)

### Deriva de energía (relax)

![Energy drift relax](energy_drift_relax.png)

### Animaciones (relax)

| Nombre | Descripción |
|---|---|
| ![Joint anim](anim_joint.gif) | Heatmap $h(x,y)$ evolucionando. |
| ![Marginals anim](anim_marginals.gif) | Marginales convergiendo. |
| ![Scatter anim](anim_scatter.gif) | Posiciones de 5000 partículas evolucionando. |
| ![pscatter anim](anim_pscatter.gif) | Momentos $(p_x, p_y)$ relajando. |

## Régimen periódico ($\alpha = 0, \sigma_L = 0$)

Corrida: $N = 2^{16} = 65\,536$, 500 k pasos. Energía conservada
**bit-for-bit** en los 11 snapshots (cf. test de periodicidad).

### Dashboard

![Dashboard periodic](dashboard_periodic.png)

### Deriva de energía (periodic)

Idealmente cero (energía bit-conservada). Lo verifica visualmente:

![Energy drift periodic](energy_drift_periodic.png)

### Animaciones (periodic)

| Nombre | Descripción |
|---|---|
| ![Joint anim](anim_joint_periodic.gif) | Heatmap $h(x,y)$ — sin disipación, no relaja a uniforme exacto sino a phase-mixed. |
| ![Scatter anim](anim_scatter_periodic.gif) | Posiciones de 3000 partículas en régimen determinista. |
| ![pscatter anim](anim_pscatter_periodic.gif) | Momentos $(p_x, p_y)$ — los signos cambian pero las magnitudes no. |
