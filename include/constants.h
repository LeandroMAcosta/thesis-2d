/* Physical and numerical constants for the 2D event-driven gas
 * simulation. The 2D model extends the 1D parent project
 * (../thesis/) keeping the same overall philosophy: ballistic motion
 * + stochastic wall reflections + analytical event-driven
 * integration. */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define PI 3.14159265358979323846

/* Physics scales — same as 1D parent (per-component, since we model
 * each component p_x, p_y as if it were the 1D momentum). */
#define P_SIGMA  5.24684e-24     /* gaussian width per component (kg·m/s)   */
#define P_RANGE  3.0e-23         /* histogram half-range per component       */
#define P_DELTA  6.0e-26         /* momentum bin scale factor                */

/* Initial state and reflection geometry (2D box [-0.5, 0.5]^2). */
#define X_INIT_HALF       0.5    /* fills x, y in [-X_INIT_HALF, X_INIT_HALF) */
#define X_REFLECT_THRESH  0.502
#define X_REFLECT_FACTOR  1.004

/* Histogram trimming for chi^2 over momentum components.
 * Same fraction as 1D parent (BORDES=237 over BINS=500 → ~47%). */
#define BORDES 237

/* Numerical guards. */
#define LOG_GUARD 1.0e-35

/* Cap per-particle event count (panic mode for runaway bounces). */
#define EVENTS_MAX 1048576L

#endif
