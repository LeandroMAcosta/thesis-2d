#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "state.h"

void histogram_accumulate(SimState *s);
int  histogram_write(SimState *s, const char *filename, double Et);

#endif
