#ifndef PHYSICS_SUPER_H
#define PHYSICS_SUPER_H

#include <stdint.h>

/* Super-event engine: O(1)-per-batch closed-form aggregation of the
 * event-driven dynamics, selected at runtime with SIM_ENGINE=super
 * (the default engine remains the exact per-bounce event loop).
 *
 * Per axis, over a batch with n ≫ 1 bounces:
 *   - p² gains D·S with S = Σ U(-½,½) → N(0, n/12) by CLT;
 *   - the unfolded coordinate gains (L·D/2p²)·A from velocity
 *     diffusion, A jointly Gaussian with S (Var n³/36, Cov n²/24);
 *   - tangential kicks from the OTHER axis's bounces aggregate to
 *     N(0, n_other·σ_L²)  (in 2D the wall kick lands on the other
 *     coordinate; the bouncing coordinate is snapped exactly);
 *   - the bounce count and final fold come from exact trajectory
 *     unfolding, so the conservative limit (alfa = sigma_l = 0) is
 *     exact, not approximate.
 *
 * Validated against the exact loop in the 1D parent repo
 * (thesis/event-driven/experiments/superevent_proto.c): chi² ≈ 1
 * between methods, energy drift within sampling noise. */

/* Reads SIM_ENGINE once. Returns 1 when SIM_ENGINE=super. */
int physics_super_enabled(void);

/* Evolve one particle for time T using the super-event closed form.
 * Returns 1 if handled; 0 when any axis violates the CLT guards (too
 * few bounces, |p| near the clamp region). On 0 the seed is untouched
 * and the caller must run the exact event loop. */
int physics_super_particle(double *x_io, double *y_io,
                           double *px_io, double *py_io, uint32_t *seed,
                           double T, double M, double sigmaL,
                           double alfa, double pmin, double pmax);

#endif
