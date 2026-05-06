#ifndef IO_H
#define IO_H

#include "state.h"

void io_load_parameters(const char *filename, SimParams *params);
void io_read_state(const char *filename, SimState *s);
void io_write_state(const char *filename, const SimState *s);

#endif
