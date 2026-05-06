#ifndef PHYSICS_H
#define PHYSICS_H

#include "state.h"

void physics_init_distributions(SimState *s);
void physics_init_state_random(SimState *s);

/* Evolve every particle for `n_steps * DT` seconds of simulated time
 * via 2D event-driven integration. */
void physics_evolve_batch(SimState *s, int n_steps, unsigned int batch_index);

double physics_energy(const SimState *s);

#endif
