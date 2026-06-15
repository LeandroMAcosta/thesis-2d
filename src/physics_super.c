/* Super-event engine: O(1)-per-batch aggregation of the 2D event-driven
 * dynamics. See include/physics_super.h for the model summary; the
 * derivation and validation harness live in the 1D parent repo
 * (thesis/event-driven/experiments/superevent_proto.c). */

#include "constants.h"
#include "physics_super.h"
#include "rng.h"

#include <math.h>
#include <stdlib.h>

/* Minimum bounce count for the CLT aggregation. Sums of uniforms are
 * close to Gaussian already at n ≈ 12; 16 keeps a margin and particles
 * below it run the exact loop cheaply anyway. */
#define SUPER_N_MIN 16L

int physics_super_enabled(void)
{
    static int cached = -1;
    if (cached < 0) {
        const char *e = getenv("SIM_ENGINE");
        cached = (e != NULL && e[0] == 's') ? 1 : 0;
    }
    return cached;
}

/* One full Box-Muller pair (2 uniforms → 2 independent N(0,1)). */
static inline void gauss_pair(uint32_t *seed, double *g1, double *g2)
{
    double r1 = rng_xorshift32(seed);
    double r2 = rng_xorshift32(seed);
    double R = sqrt(-2.0 * log(r1 + LOG_GUARD));
    double th = 2.0 * PI * r2;
    *g1 = R * cos(th);
    *g2 = R * sin(th);
}

/* Per-axis pre-pass: bounce count + CLT guards. No RNG consumed.
 * Returns -1 when the axis must fall back to the exact loop. */
static inline long axis_count_and_guard(double x, double p, double T,
                                        double M, double alfa,
                                        double pmin, double pmax,
                                        double *u_det_out, double *D_out,
                                        double *p2_out)
{
    if (p == 0.0) return -1; /* frozen axis: exact loop handles it */

    double v0 = p / M;
    double w0 = x + 0.5;
    double u_det = w0 + v0 * T;
    double lo = fmin(w0, u_det), hi = fmax(w0, u_det);
    long n = (long)floor(hi) - (long)floor(lo);
    if (n < SUPER_N_MIN) return -1;

    double p_abs = fabs(p);
    double p2 = p_abs * p_abs;
    double D = alfa * (p_abs - pmin) * (pmax - p_abs);
    double sigS = sqrt((double)n / 12.0);
    if (6.0 * fabs(D) * sigS > 0.5 * p2) return -1;  /* p² clamp risk */
    if (1.5 * fabs(D) * sigS / p2 > 0.05) return -1; /* D drifts      */

    *u_det_out = u_det;
    *D_out = D;
    *p2_out = p2;
    return n;
}

/* Per-axis commit: draw the aggregated noise, fold, write back.
 * kick_var = (Σ other axes' bounce counts) · σ_L². */
static inline void axis_apply(double *x_io, double *p_io, uint32_t *seed,
                              long n, double u_det, double D, double p2,
                              double kick_var)
{
    const double L = 1.0;
    double v0_sign = copysign(1.0, *p_io);

    double u = u_det;
    double p_absf = fabs(*p_io); /* conservative limit: magnitude bit-exact */

    if (D != 0.0) {
        /* Correlated pair (S, A): endpoint and integral of the U-walk.
         *   S = σ√n Z₁
         *   A = σ n^{3/2} (Z₁/2 + Z₂/√12),  σ = 1/√12. */
        double Z1, Z2;
        gauss_pair(seed, &Z1, &Z2);
        const double SIG = 0.28867513459481287; /* 1/sqrt(12) */
        double sqn = sqrt((double)n);
        double S = SIG * sqn * Z1;
        double A = SIG * (double)n * sqn * (0.5 * Z1 + Z2 * SIG);

        double p2f = p2 + D * S;
        if (p2f < 0.0) p2f = 0.0; /* unreachable given the guards */
        p_absf = sqrt(p2f);

        /* Velocity-diffusion phase (correlated with S, along motion). */
        u += v0_sign * (L * D / (2.0 * p2)) * A;
    }

    if (kick_var != 0.0) {
        /* Cross-axis kick aggregate (independent). */
        double Zk, Zspare;
        gauss_pair(seed, &Zk, &Zspare);
        (void)Zspare;
        u += sqrt(kick_var) * Zk;
    }

    /* Fold back with mirror period 2L. */
    double m2 = fmod(u, 2.0);
    if (m2 < 0.0) m2 += 2.0;
    double wf;
    int mirrored;
    if (m2 <= 1.0) { wf = m2;       mirrored = 0; }
    else           { wf = 2.0 - m2; mirrored = 1; }

    *x_io = wf - 0.5;
    *p_io = copysign(p_absf, mirrored ? -v0_sign : v0_sign);
}

/* Exact per-axis mini-loop for axes that fail the CLT guards. In 2D a
 * wall hit never kicks its own coordinate (the bouncing coordinate is
 * snapped to the wall exactly), so a single axis evolves exactly with
 * one uniform draw per bounce — no Box-Muller in the loop. The kicks
 * it receives from the OTHER axis's bounces are aggregated into one
 * Gaussian applied at the end (kick_var), mirroring the position into
 * the box without touching p — same behaviour as the exact engine's
 * reflection guard, which reflects position only. */
static inline void axis_mini_loop(double *x_io, double *p_io, uint32_t *seed,
                                  double T, double M, double alfa,
                                  double pmin, double pmax, double kick_var)
{
    double x = *x_io, p = *p_io, t = 0.0;
    const double invM = 1.0 / M;

    while (p != 0.0) {
        double v = p * invM;
        double v_sign = (v > 0.0) ? 1.0 : -1.0;
        double target = copysign(0.5, v);
        double tw = (target - x) / v;
        if (tw < 0.0) tw = 0.0;
        if (t + tw >= T) { x += v * (T - t); break; }
        t += tw;
        x = target;

        double p_abs = fabs(p);
        if (alfa != 0.0) {
            double D = alfa * (p_abs - pmin) * (pmax - p_abs);
            double r = rng_xorshift32(seed);
            double sq = p_abs * p_abs + D * (r - 0.5);
            p_abs = sqrt(sq > 0.0 ? sq : 0.0);
        }
        p = -v_sign * p_abs;
    }

    if (kick_var != 0.0) {
        double Zk, Zspare;
        gauss_pair(seed, &Zk, &Zspare);
        (void)Zspare;
        double w = (x + 0.5) + sqrt(kick_var) * Zk;
        /* Mirror-fold the displaced position; momentum untouched. */
        double m2 = fmod(w, 2.0);
        if (m2 < 0.0) m2 += 2.0;
        x = ((m2 <= 1.0) ? m2 : 2.0 - m2) - 0.5;
    }

    *x_io = x;
    *p_io = p;
}

int physics_super_particle(double *x_io, double *y_io,
                           double *px_io, double *py_io, uint32_t *seed,
                           double T, double M, double sigmaL,
                           double alfa, double pmin, double pmax)
{
    /* Phase 1: per-axis bounce counts + guards. No RNG consumed, so a
     * full fallback re-runs the exact loop on the untouched stream. */
    double ux, Dx, p2x, uy, Dy, p2y;
    long nx = axis_count_and_guard(*x_io, *px_io, T, M, alfa, pmin, pmax,
                                   &ux, &Dx, &p2x);
    long ny = axis_count_and_guard(*y_io, *py_io, T, M, alfa, pmin, pmax,
                                   &uy, &Dy, &p2y);
    if (nx < 0 && ny < 0) return 0; /* both axes slow: coupled exact loop */

    /* Deterministic bounce-count estimates for the kick coupling (also
     * needed for guard-failing axes). */
    long nx_est = nx, ny_est = ny;
    if (nx_est < 0) {
        double v = *px_io / M, w0 = *x_io + 0.5, u = w0 + v * T;
        double lo = fmin(w0, u), hi = fmax(w0, u);
        nx_est = (long)floor(hi) - (long)floor(lo);
    }
    if (ny_est < 0) {
        double v = *py_io / M, w0 = *y_io + 0.5, u = w0 + v * T;
        double lo = fmin(w0, u), hi = fmax(w0, u);
        ny_est = (long)floor(hi) - (long)floor(lo);
    }

    /* The aggregated cross-kick is a Gaussian approximation of a sum of
     * n_other Gaussian kicks — exact in distribution for ANY count, but
     * the count itself is the deterministic estimate. Only when few
     * kicks carry non-negligible amplitude (exotic σ_L) does the
     * difference matter; production σ_L=1e-4 aggregates are ≤4e-4,
     * far below a histogram bin. */
    if (sigmaL != 0.0 &&
        ((ny_est < SUPER_N_MIN && sigmaL * sqrt((double)ny_est) > 0.01) ||
         (nx_est < SUPER_N_MIN && sigmaL * sqrt((double)nx_est) > 0.01)))
        return 0;

    /* Scramble the seed before drawing. Super consumes only the first
     * few values of each particle's stream, where the linear structure
     * of rng_seed_mix is still visible across particles — without
     * this, population aggregates pick up a systematic component
     * (measured −37 ppm energy drift in the 1D repo; see the harness
     * referenced in the header). A murmur3 finalizer decorrelates the
     * shallow draws. The exact engine consumes draws at hundreds of
     * varying depths and is unaffected. */
    uint32_t z = *seed;
    z ^= z >> 16; z *= 0x85ebca6bu;
    z ^= z >> 13; z *= 0xc2b2ae35u;
    z ^= z >> 16;
    *seed = z ? z : 1u;

    /* Phase 2: super for fast axes, exact mini-loop for slow ones.
     * Cross-axis kick coupling: x accumulates n_y kicks, y gets n_x. */
    double sl2 = sigmaL * sigmaL;
    if (nx >= 0)
        axis_apply(x_io, px_io, seed, nx, ux, Dx, p2x, (double)ny_est * sl2);
    else
        axis_mini_loop(x_io, px_io, seed, T, M, alfa, pmin, pmax,
                       (double)ny_est * sl2);
    if (ny >= 0)
        axis_apply(y_io, py_io, seed, ny, uy, Dy, p2y, (double)nx_est * sl2);
    else
        axis_mini_loop(y_io, py_io, seed, T, M, alfa, pmin, pmax,
                       (double)nx_est * sl2);
    return 1;
}
