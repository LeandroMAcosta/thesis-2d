#ifndef BACKEND_H
#define BACKEND_H

#include "state.h"

void backend_init(SimState *s);
void backend_finalize(SimState *s);
void backend_set_threads(int n);
const char *backend_name(void);

#endif
