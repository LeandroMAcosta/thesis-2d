/* Simulation state for the 2D model. Each particle has 4 scalar
 * variables (x, y, p_x, p_y) instead of the 1D model's 2 (x, p).
 *
 * Histograms:
 *   hx, hy      — 1D marginals of position over each axis (size 2*BINS+4)
 *   gpx, gpy    — 1D marginals of momentum over each component (size 2*BINS)
 *   h_xy        — 2D joint histogram of position (size HXY_BINS × HXY_BINS,
 *                 a coarser bineo since 500x500 = 250k cells is heavy)
 *
 * For chi² we use only the marginals (hx, hy, gpx, gpy) plus the
 * isotropy / independence checks computed directly from particle
 * arrays in postcondition scripts. The 2D joint h_xy is used only
 * for visualization and for the "h(x,y) phase-mixes uniformly"
 * test in the periodic regime. */

#ifndef STATE_H
#define STATE_H

#include <stdbool.h>
#include <stdint.h>

#define MAX_TANDAS    500
#define MAX_FILENAME  256
#define HXY_BINS      64    /* coarse 2D joint histogram, 64×64 = 4096 cells */

typedef struct {
    int N_PART;
    int BINS;
    int N_THREADS;
    double DT;
    double M;
    double sigmaL;
    double alfa;
    double pmin;
    double pmax;

    unsigned int Ntandas;
    int steps[MAX_TANDAS];

    char inputFilename[MAX_FILENAME];
    char saveFilename[MAX_FILENAME];
    bool should_resume;
    bool should_dump;
    bool keep_snapshots;
    bool discretize_momenta;

    uint64_t seed;
} SimParams;

typedef struct {
    SimParams params;

    double *x;     /* N_PART */
    double *y;     /* N_PART */
    double *px;    /* N_PART */
    double *py;    /* N_PART */

    int *hx;       /* 2*BINS + 4 */
    int *hy;       /* 2*BINS + 4 */
    int *gpx;      /* 2*BINS     */
    int *gpy;      /* 2*BINS     */

    double *DxE;   /* 2*BINS + 4 — reference uniform distribution    */
    double *DpE;   /* 2*BINS     — reference Maxwell-Boltzmann       */

    int *h_xy;     /* HXY_BINS × HXY_BINS — coarse 2D joint           */

    unsigned int evolution;
} SimState;

void state_alloc(SimState *s);
void state_free(SimState *s);
void state_zero_histograms(SimState *s);

#endif
