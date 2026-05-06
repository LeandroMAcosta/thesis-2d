#include "backend.h"

#include <omp.h>

void backend_init(SimState *s)      { (void)s; }
void backend_finalize(SimState *s)  { (void)s; }

void backend_set_threads(int n) { omp_set_num_threads(n); }

const char *backend_name(void)  { return "2d-omp-event"; }
