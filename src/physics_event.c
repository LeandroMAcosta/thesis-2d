/* 2D event-driven CPU + OpenMP backend.
 *
 * The 2D model: each particle (x, y, p_x, p_y) lives in a square box
 * [-0.5, 0.5]^2. Between wall hits it moves ballistically. The next
 * event is the first of four candidate wall crossings:
 *
 *   t_x+ = ( 0.5 - x) / v_x   if v_x > 0
 *   t_x- = (-0.5 - x) / v_x   if v_x < 0
 *   t_y+ = ( 0.5 - y) / v_y   if v_y > 0
 *   t_y- = (-0.5 - y) / v_y   if v_y < 0
 *
 * Of the four, exactly two are finite (positive) at any given time —
 * those that match the sign of v_x and v_y. We take the minimum.
 *
 * On a vertical-wall hit (x = ±0.5):
 *   - p_x flips sign, |p_x| is perturbed via the alfa-loop;
 *   - y receives a Gaussian kick (Box-Muller, scaled by sigmaL).
 * On a horizontal-wall hit (y = ±0.5):
 *   - p_y flips sign, |p_y| is perturbed;
 *   - x receives the kick.
 *
 * This is the natural 2D extension of the 1D model: the wall affects
 * the perpendicular momentum component (energy exchange) and induces
 * a tangential position kick (rugosidad). The other component is
 * untouched. */

#include "constants.h"
#include "physics.h"
#include "rng.h"
#include "state.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

/* ---------------------------------------------------------------------------
 * Initial reference distributions (for chi² in histogram_omp.c).
 * --------------------------------------------------------------------------- */

void physics_init_distributions(SimState *s)
{
    int N = s->params.N_PART;
    int B = s->params.BINS;

    /* Per-component Maxwell-Boltzmann. Normalised to N counts (each
     * particle contributes one count per histogram, since the marginal
     * over p_x sums all N particles). */
    double numerator   = P_DELTA * (double)N;
    double denominator = P_SIGMA * sqrt(2.0 * PI);
    for (int i = 0; i < 2 * B; i++) {
        double z = P_RANGE * ((double)i / B - 0.999) / P_SIGMA;
        s->DpE[i] = (numerator / denominator) * exp(-(z * z) / 2.0);
    }

    /* Per-axis uniform on [-0.5, 0.5]. The marginal over x has N
     * particles spread in 2*BINS interior bins → N / (2*BINS) per
     * bin times 1/L = 1 → DX_BIN_DENSITY = 1/(2*BINS).
     * To match parent project's normalisation we use N / (2*BINS) per
     * interior bin. */
    int n = 2 * B + 4;
    double bin_density = (double)N / (double)(2 * B);
    for (int i = 0; i < n; i++) {
        if (i == 0 || i == 1 || i == 2 * B + 2 || i == 2 * B + 3) {
            s->DxE[i] = 0.0;
        } else {
            s->DxE[i] = bin_density;
        }
    }
}

/* ---------------------------------------------------------------------------
 * Initial random state.
 * --------------------------------------------------------------------------- */

void physics_init_state_random(SimState *s)
{
    int N = s->params.N_PART;
    uint64_t base_xy = s->params.seed;
    uint64_t base_p  = s->params.seed ^ 0xDEADBEEFCAFEBABEull;

    /* Positions: uniform in the upper-right quadrant [0, X_INIT_HALF)^2.
     * Mimics the 1D parent project (which fills only [0, 0.5)) and gives
     * visually dramatic "diffusion from a corner" demos. The system
     * relaxes to uniform fill of the whole box at equilibrium. */
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        uint32_t seed = rng_seed_mix(base_xy, (uint32_t)i);
        double rx = rng_xorshift32(&seed);
        double ry = rng_xorshift32(&seed);
        s->x[i] = rx * X_INIT_HALF;
        s->y[i] = ry * X_INIT_HALF;
    }

    if (s->params.init_uniform_p_magnitude && s->params.uniform_p_magnitude > 0.0) {
        /* Cycle-demo init: all particles with |p_x| = |p_y| = p_0,
         * random ± signs. With p_0 chosen so that 2L/v_perp divides
         * evenly into the run length, the system returns EXACTLY to
         * its initial state at the end of the run (Poincaré recurrence
         * by construction). */
        double p0 = s->params.uniform_p_magnitude;
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) {
            uint32_t seed = rng_seed_mix(base_p, (uint32_t)i);
            uint32_t bits = (uint32_t)(rng_xorshift32(&seed) * 4.0);
            int sx = (bits & 1) ? +1 : -1;
            int sy = (bits & 2) ? +1 : -1;
            s->px[i] = sx * p0;
            s->py[i] = sy * p0;
        }
    } else {
        /* Default: per-component gaussian via Box-Muller on pairs. */
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) {
            uint32_t seed = rng_seed_mix(base_p, (uint32_t)i);
            double r1 = rng_xorshift32(&seed);
            double r2 = rng_xorshift32(&seed);
            double r3 = rng_xorshift32(&seed);
            double r4 = rng_xorshift32(&seed);
            double xi1 = sqrt(-2.0 * log(r1 + LOG_GUARD));
            double xi2 = 2.0 * PI * r2;
            double xi3 = sqrt(-2.0 * log(r3 + LOG_GUARD));
            double xi4 = 2.0 * PI * r4;
            s->px[i] = xi1 * cos(xi2) * P_SIGMA;
            s->py[i] = xi3 * cos(xi4) * P_SIGMA;
        }
    }
}

/* ---------------------------------------------------------------------------
 * Bounce: applied when a wall is hit.
 * Side: 0 = vertical wall (x = ±0.5), 1 = horizontal wall (y = ±0.5).
 * v_perp_sign: sign of the perpendicular velocity component pre-bounce.
 * --------------------------------------------------------------------------- */

static inline void apply_bounce(
    double *x, double *y, double *px, double *py,
    uint32_t *seed, double sigmaL, double alfa,
    double pmin, double pmax,
    int side, double v_perp_sign)
{
    /* Tangential kick (Gaussian) on the OTHER coordinate. */
    double r1 = rng_xorshift32(seed);
    double r2 = rng_xorshift32(seed);
    double xi1 = sqrt(-2.0 * log(r1 + LOG_GUARD));
    double xi2 = 2.0 * PI * r2;
    double delta_tan = xi1 * cos(xi2) * sigmaL;
    if (fabs(delta_tan) > 1.0) delta_tan = copysign(1.0, delta_tan);

    /* Energy perturbation on the PERPENDICULAR component magnitude. */
    double *p_perp = (side == 0) ? px : py;
    double *t_pos  = (side == 0) ? y  : x;
    double *p_perp_at_wall = (side == 0) ? x : y;  /* the one fixed at ±0.5 */

    /* Apply tangential kick on the tangential coordinate. */
    *t_pos += delta_tan;
    if (!isfinite(*t_pos) || fabs(*t_pos) > X_REFLECT_THRESH) {
        *t_pos = X_REFLECT_FACTOR * copysign(1.0, *t_pos) - *t_pos;
        if (!isfinite(*t_pos) || fabs(*t_pos) > 0.5) *t_pos = 0.0;
    }

    /* alfa-loop on |p_perp|, single iteration (k = 1 always). */
    double p_abs = fabs(*p_perp);
    if (alfa != 0.0) {
        double DeltaE = alfa * (p_abs - pmin) * (pmax - p_abs);
        double r3 = rng_xorshift32(seed);
        double sq = p_abs * p_abs + DeltaE * (r3 - 0.5);
        p_abs = sqrt(sq > 0.0 ? sq : 0.0);
    }
    /* Flip sign of the perpendicular component. */
    *p_perp = -v_perp_sign * p_abs;

    /* Snap the perpendicular coordinate back into the box just in case
     * (it was set to ±0.5 exactly at the bounce; this guard is for
     * safety against subsequent FP noise pushing it out by ULPs). */
    (void)p_perp_at_wall;  /* kept for documentation only */
}

/* ---------------------------------------------------------------------------
 * Per-particle event-driven loop.
 * --------------------------------------------------------------------------- */

static inline void evolve_one_particle(
    double *x_io, double *y_io, double *px_io, double *py_io,
    uint32_t *seed, double T_target, double M, double sigmaL,
    double alfa, double pmin, double pmax)
{
    double x  = *x_io;
    double y  = *y_io;
    double px = *px_io;
    double py = *py_io;
    double t  = 0.0;
    long event_count = 0;

    /* Hoist 1/M outside the loop: same micro-opt as the 3D variant. */
    const double invM = 1.0 / M;

    while (1) {
        /* If both momenta are zero the particle is frozen. */
        if (px == 0.0 && py == 0.0) break;

        double vx = px * invM;
        double vy = py * invM;

        /* Branchless wall targets: copysign(0.5, v) picks ±0.5 with the
         * sign of v. v=0 (impossible in practice) would still produce
         * a target, but the other axis wins the min anyway. */
        double x_target = copysign(0.5, vx);
        double y_target = copysign(0.5, vy);
        double t_x = (x_target - x) / vx;
        double t_y = (y_target - y) / vy;

        /* Edge: numerical drift can make the candidate negative. */
        if (t_x < 0.0) t_x = 0.0;
        if (t_y < 0.0) t_y = 0.0;

        double t_event = (t_x < t_y) ? t_x : t_y;
        int side = (t_x < t_y) ? 0 : 1;  /* 0 = vertical wall, 1 = horizontal */

        if (t + t_event >= T_target) {
            /* Final ballistic segment, no more bounces this batch. */
            double dt = T_target - t;
            x += vx * dt;
            y += vy * dt;
            break;
        }

        /* Advance to the wall. */
        t += t_event;
        x += vx * t_event;
        y += vy * t_event;

        /* Snap the perpendicular coordinate to the wall exactly to
         * avoid FP drift (the next iteration's t_event would be 0
         * for that side otherwise). */
        if (side == 0) x = x_target;
        else           y = y_target;

        double v_perp_sign = (side == 0) ? ((vx > 0) ? 1.0 : -1.0)
                                          : ((vy > 0) ? 1.0 : -1.0);

        apply_bounce(&x, &y, &px, &py, seed, sigmaL, alfa,
                     pmin, pmax, side, v_perp_sign);

        event_count++;
        if (event_count > EVENTS_MAX) break;
    }

    *x_io  = x;
    *y_io  = y;
    *px_io = px;
    *py_io = py;
}

/* ---------------------------------------------------------------------------
 * Batch evolution.
 * --------------------------------------------------------------------------- */

void physics_evolve_batch(SimState *s, int n_steps, unsigned int batch_index)
{
    int N = s->params.N_PART;
    double T_target = (double)n_steps * s->params.DT;
    uint64_t base = s->params.seed
                  ^ ((uint64_t)batch_index * 0x9E3779B97F4A7C15ull);

    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < N; i++) {
        uint32_t seed = rng_seed_mix(base, (uint32_t)i);
        double xi  = s->x[i];
        double yi  = s->y[i];
        double pxi = s->px[i];
        double pyi = s->py[i];
        evolve_one_particle(&xi, &yi, &pxi, &pyi, &seed, T_target,
                            s->params.M, s->params.sigmaL,
                            s->params.alfa, s->params.pmin, s->params.pmax);
        s->x[i]  = xi;
        s->y[i]  = yi;
        s->px[i] = pxi;
        s->py[i] = pyi;
    }
}

double physics_energy(const SimState *s)
{
    int N = s->params.N_PART;
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (int i = 0; i < N; i++) {
        sum += s->px[i] * s->px[i] + s->py[i] * s->py[i];
    }
    double Et = sum / (2.0 * s->params.M);
    printf("N° de pasos %6u\tEnergía total = %12.9E\n", s->evolution, Et);
    fflush(stdout);
    return Et;
}
