/* 2D histogram accumulation, chi² metrics, output writer.
 *
 * We bin four 1D marginals (hx, hy, gpx, gpy) and one coarse 2D joint
 * (h_xy). The .dat file format extends the 1D parent's: it has the
 * marginals as columns and the joint as a separate block.
 *
 * chi² metrics computed in the file header:
 *   chi2x  — h(x) marginal vs uniform reference
 *   chi2y  — h(y) marginal vs uniform reference
 *   chi2px — g(p_x) marginal vs Maxwell-Boltzmann reference
 *   chi2py — g(p_y) marginal vs Maxwell-Boltzmann reference
 *   isotropy_ratio — <p_x²>/<p_y²>, should be close to 1
 *   correlation    — <p_x p_y> / sqrt(<p_x²><p_y²>), should be close to 0 */

#include "histogram.h"
#include "constants.h"
#include "state.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#define X_BIN_SCALE 1.99999999999999
#define P_BIN_SCALE 0.999999999999994

void histogram_accumulate(SimState *s)
{
    int N = s->params.N_PART;
    int B = s->params.BINS;

    state_zero_histograms(s);

    /* Bin into per-thread local arrays first to avoid atomic contention,
     * then sum locally at the end. For simplicity we use #pragma omp atomic
     * — at N=2^20 with BINS=500 this is fast enough (a few ms). */

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        /* hx and hy use the same binning convention as the 1D parent. */
        int hx_idx = (int)floor((s->x[i] + 0.5) * (X_BIN_SCALE * B) + 2.0);
        int hy_idx = (int)floor((s->y[i] + 0.5) * (X_BIN_SCALE * B) + 2.0);
        if (hx_idx < 0)            hx_idx = 0;
        if (hx_idx > 2 * B + 3)    hx_idx = 2 * B + 3;
        if (hy_idx < 0)            hy_idx = 0;
        if (hy_idx > 2 * B + 3)    hy_idx = 2 * B + 3;

        int gpx_idx = (int)floor((s->px[i] / P_RANGE + 1.0) * (P_BIN_SCALE * B));
        int gpy_idx = (int)floor((s->py[i] / P_RANGE + 1.0) * (P_BIN_SCALE * B));
        if (gpx_idx < 0)           gpx_idx = 0;
        if (gpx_idx > 2 * B - 1)   gpx_idx = 2 * B - 1;
        if (gpy_idx < 0)           gpy_idx = 0;
        if (gpy_idx > 2 * B - 1)   gpy_idx = 2 * B - 1;

        /* Coarse 2D joint. Map x,y in [-0.5, 0.5] to [0, HXY_BINS). */
        int hxy_x = (int)floor((s->x[i] + 0.5) * HXY_BINS);
        int hxy_y = (int)floor((s->y[i] + 0.5) * HXY_BINS);
        if (hxy_x < 0)             hxy_x = 0;
        if (hxy_x >= HXY_BINS)     hxy_x = HXY_BINS - 1;
        if (hxy_y < 0)             hxy_y = 0;
        if (hxy_y >= HXY_BINS)     hxy_y = HXY_BINS - 1;

        #pragma omp atomic update
        s->hx[hx_idx]++;
        #pragma omp atomic update
        s->hy[hy_idx]++;
        #pragma omp atomic update
        s->gpx[gpx_idx]++;
        #pragma omp atomic update
        s->gpy[gpy_idx]++;
        #pragma omp atomic update
        s->h_xy[hxy_y * HXY_BINS + hxy_x]++;
    }
}

/* Compute <p_x²>, <p_y²>, <p_x p_y> on the particle arrays for the
 * isotropy and independence checks. */
static void moments(SimState *s, double *m_pxx, double *m_pyy, double *m_pxy)
{
    int N = s->params.N_PART;
    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    #pragma omp parallel for reduction(+:sxx,syy,sxy) schedule(static)
    for (int i = 0; i < N; i++) {
        sxx += s->px[i] * s->px[i];
        syy += s->py[i] * s->py[i];
        sxy += s->px[i] * s->py[i];
    }
    *m_pxx = sxx / (double)N;
    *m_pyy = syy / (double)N;
    *m_pxy = sxy / (double)N;
}

int histogram_write(SimState *s, const char *filename, double Et)
{
    int B = s->params.BINS;
    int *hx  = s->hx;
    int *hy  = s->hy;
    int *gpx = s->gpx;
    int *gpy = s->gpy;
    double *DxE = s->DxE;
    double *DpE = s->DpE;

    /* Two independent chi² for x and y marginals (each vs uniform). */
    double chi2x = 0.0, chi2y = 0.0;
    for (int i = 4; i < 2 * B; i++) {
        double dx = hx[i] - DxE[i];
        double dy = hy[i] - DxE[i];
        chi2x += (dx * dx) / DxE[i];
        chi2y += (dy * dy) / DxE[i];
    }
    chi2x /= (2.0 * B - 4);
    chi2y /= (2.0 * B - 4);

    /* chi² for momentum marginals (each vs MB reference). */
    double chi2px = 0.0, chi2py = 0.0;
    for (int i = 0; i < 2 * (B - BORDES); i++) {
        double dx = gpx[i + BORDES] - DpE[i + BORDES];
        double dy = gpy[i + BORDES] - DpE[i + BORDES];
        chi2px += (dx * dx) / DpE[i + BORDES];
        chi2py += (dy * dy) / DpE[i + BORDES];
    }
    chi2px /= (2.0 * (B - BORDES));
    chi2py /= (2.0 * (B - BORDES));

    /* Isotropy and correlation from particle arrays. */
    double m_pxx, m_pyy, m_pxy;
    moments(s, &m_pxx, &m_pyy, &m_pxy);
    double isotropy_ratio = (m_pyy > 0) ? (m_pxx / m_pyy) : 1.0;
    double correlation    = (m_pxx > 0 && m_pyy > 0)
                          ? (m_pxy / sqrt(m_pxx * m_pyy)) : 0.0;

    FILE *out = fopen(filename, "w");
    if (out == NULL) {
        fprintf(stderr, "histogram_write: cannot open %s\n", filename);
        return 1;
    }
    /* Header line — extends the 1D format with 2D-specific metrics. */
    fprintf(out,
            "# 2D-event-driven  chi2x =%9.6f  chi2y =%9.6f  "
            "chi2px =%9.6f  chi2py =%9.6f  "
            "isotropy =%9.6f  corr =%12.6e  Et=%12.9E  N=%d  BINS=%d  HXY_BINS=%d\n",
            chi2x, chi2y, chi2px, chi2py,
            isotropy_ratio, correlation, Et,
            s->params.N_PART, B, HXY_BINS);

    /* Marginals block: 2*B + 4 rows with x_val, hx, y_val, hy, p_val, gpx, gpy. */
    fprintf(out, "# MARGINALS\n");
    fprintf(out, "# x_val  hx  y_val  hy  p_val  gpx  gpy\n");
    fprintf(out, "%8.5f %6d %8.5f %6d %24.12E %6d %6d\n",
            -0.5015, hx[0], -0.5015, hy[0], -2.997e-23, gpx[0], gpy[0]);
    fprintf(out, "%8.5f %6d %8.5f %6d %24.12E %6d %6d\n",
            -0.5005, hx[1], -0.5005, hy[1], -2.997e-23, gpx[0], gpy[0]);
    for (int i = 0; i < 2 * B; i++) {
        fprintf(out, "%8.5f %6d %8.5f %6d %24.12E %6d %6d\n",
                (0.5 * i / B - 0.4995),
                hx[i + 2],
                (0.5 * i / B - 0.4995),
                hy[i + 2],
                (P_RANGE * i / B - 2.997e-23),
                gpx[i],
                gpy[i]);
    }
    fprintf(out, "%8.5f %6d %8.5f %6d %24.12E %6d %6d\n",
            0.5005, hx[2 * B + 2], 0.5005, hy[2 * B + 2],
            2.997e-23, gpx[2 * B - 1], gpy[2 * B - 1]);
    fprintf(out, "%8.5f %6d %8.5f %6d %24.12E %6d %6d\n",
            0.5015, hx[2 * B + 3], 0.5015, hy[2 * B + 3],
            2.997e-23, gpx[2 * B - 1], gpy[2 * B - 1]);

    /* Coarse 2D joint h_xy block: HXY_BINS rows of HXY_BINS counts. */
    fprintf(out, "# JOINT h(x,y) — HXY_BINS x HXY_BINS\n");
    for (int j = 0; j < HXY_BINS; j++) {
        for (int i = 0; i < HXY_BINS; i++) {
            fprintf(out, "%6d ", s->h_xy[j * HXY_BINS + i]);
        }
        fprintf(out, "\n");
    }

    fclose(out);
    return 0;
}
