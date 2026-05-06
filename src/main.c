/* Main loop. Same control flow as the 1D parent: load config, init
 * state, run batches, snapshot histogram + checkpoint after each. */

#include "backend.h"
#include "histogram.h"
#include "io.h"
#include "physics.h"
#include "state.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void make_filename(char *out, size_t cap, unsigned int evolution)
{
    if (evolution < 10000000u) {
        snprintf(out, cap, "X%07u.dat", evolution);
    } else {
        snprintf(out, cap, "X%1.3e.dat", (double)evolution);
        char *e = memchr(out, 'e', cap);
        if (e != NULL) memmove(e + 1, e + 3, strlen(e + 3) + 1);
    }
}

int main(int argc, char **argv)
{
    setvbuf(stdout, NULL, _IOLBF, 0);

    const char *config = (argc > 1) ? argv[1] : "config.toml";

    SimState s = (SimState){0};

    const char *seed_env = getenv("SIM_SEED");
    if (seed_env != NULL) {
        s.params.seed = strtoull(seed_env, NULL, 10);
    } else {
        s.params.seed = (uint64_t)time(NULL);
    }

    io_load_parameters(config, &s.params);
    backend_set_threads(s.params.N_THREADS);

    state_alloc(&s);
    backend_init(&s);
    physics_init_distributions(&s);

    if (s.params.should_resume) {
        io_read_state(s.params.inputFilename, &s);
    } else {
        physics_init_state_random(&s);
        histogram_accumulate(&s);
        double Et = physics_energy(&s);
        histogram_write(&s, "X0000000.dat", Et);
    }

    printf("pmin=%12.9E      alfa=%12.9E   backend=%s\n",
           s.params.pmin, s.params.alfa, backend_name());
    fflush(stdout);

    char filename[MAX_FILENAME];
    for (unsigned int j = 0; j < s.params.Ntandas; j++) {
        physics_evolve_batch(&s, s.params.steps[j], j + 1);
        histogram_accumulate(&s);

        s.evolution += s.params.steps[j];
        make_filename(filename, sizeof(filename), s.evolution);

        if (s.params.should_dump) {
            io_write_state(s.params.saveFilename, &s);
            if (s.params.keep_snapshots) {
                char snap_name[MAX_FILENAME + 16];
                snprintf(snap_name, sizeof(snap_name), "%s.%07u",
                         s.params.saveFilename, s.evolution);
                io_write_state(snap_name, &s);
            }
        }

        double Et = physics_energy(&s);
        histogram_write(&s, filename, Et);
    }

    printf("Completo evolution = %u\n", s.evolution);

    backend_finalize(&s);
    state_free(&s);
    return 0;
}
