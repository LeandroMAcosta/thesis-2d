/* TOML config loader + binary checkpoint I/O for 2D state.
 * Same TOML schema as the 1D parent; checkpoint binary stores
 * x, y, px, py instead of just x, p. */

#include "io.h"
#include "state.h"

#include "../third_party/tomlc99/toml.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void die_toml(const char *path, const char *msg)
{
    fprintf(stderr, "config %s: %s\n", path, msg);
    exit(1);
}

static int64_t require_int(toml_table_t *t, const char *key, const char *path,
                           int64_t min, int64_t max)
{
    toml_datum_t d = toml_int_in(t, key);
    if (!d.ok) {
        fprintf(stderr, "config %s: missing or non-integer key `%s`\n", path, key);
        exit(1);
    }
    if (d.u.i < min || d.u.i > max) {
        fprintf(stderr, "config %s: %s = %lld out of range [%lld, %lld]\n",
                path, key, (long long)d.u.i, (long long)min, (long long)max);
        exit(1);
    }
    return d.u.i;
}

static double require_double(toml_table_t *t, const char *key, const char *path,
                             double min, double max)
{
    toml_datum_t d = toml_double_in(t, key);
    if (!d.ok) {
        toml_datum_t di = toml_int_in(t, key);
        if (di.ok) d.u.d = (double)di.u.i;
        else {
            fprintf(stderr, "config %s: missing or non-numeric key `%s`\n", path, key);
            exit(1);
        }
    }
    if (d.u.d < min || d.u.d > max) {
        fprintf(stderr, "config %s: %s = %g out of range [%g, %g]\n",
                path, key, d.u.d, min, max);
        exit(1);
    }
    return d.u.d;
}

static bool require_bool(toml_table_t *t, const char *key, const char *path)
{
    toml_datum_t d = toml_bool_in(t, key);
    if (!d.ok) {
        fprintf(stderr, "config %s: missing or non-boolean key `%s`\n", path, key);
        exit(1);
    }
    return d.u.b != 0;
}

static void require_string(toml_table_t *t, const char *key, const char *path,
                           char *out, size_t cap)
{
    toml_datum_t d = toml_string_in(t, key);
    if (!d.ok) {
        fprintf(stderr, "config %s: missing or non-string key `%s`\n", path, key);
        exit(1);
    }
    strncpy(out, d.u.s, cap - 1);
    out[cap - 1] = '\0';
    free(d.u.s);
}

static void require_steps(toml_table_t *t, const char *path, SimParams *p)
{
    toml_array_t *arr = toml_array_in(t, "steps");
    if (arr == NULL) {
        fprintf(stderr, "config %s: missing array `steps`\n", path);
        exit(1);
    }
    int n = toml_array_nelem(arr);
    if (n <= 0 || n > MAX_TANDAS) {
        fprintf(stderr, "config %s: `steps` has %d entries (max %d)\n", path, n, MAX_TANDAS);
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        toml_datum_t d = toml_int_at(arr, i);
        if (!d.ok || d.u.i <= 0 || d.u.i > (int64_t)2147483647) {
            fprintf(stderr, "config %s: steps[%d] invalid\n", path, i);
            exit(1);
        }
        p->steps[i] = (int)d.u.i;
    }
    p->Ntandas = (unsigned int)n;
}

void io_load_parameters(const char *filename, SimParams *params)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "io_load_parameters: cannot open %s: %s\n",
                filename, strerror(errno));
        exit(1);
    }

    char errbuf[256];
    toml_table_t *root = toml_parse_file(f, errbuf, sizeof(errbuf));
    fclose(f);
    if (root == NULL) die_toml(filename, errbuf);

    params->N_PART    = (int)require_int(root, "n_part",    filename, 1, 1 << 30);
    params->BINS      = (int)require_int(root, "bins",      filename, 1, 1 << 20);
    params->DT        = require_double  (root, "dt",        filename, 1e-30, 1e10);
    params->M         = require_double  (root, "mass",      filename, 1e-40, 1.0);
    params->N_THREADS = (int)require_int(root, "n_threads", filename, 1, 4096);
    params->sigmaL    = require_double  (root, "sigma_l",   filename, 0.0, 1.0);

    require_steps(root, filename, params);

    params->should_resume = require_bool(root, "resume",      filename);
    require_string(root, "input_file",  filename, params->inputFilename, MAX_FILENAME);
    params->should_dump   = require_bool(root, "dump",        filename);
    require_string(root, "output_file", filename, params->saveFilename,  MAX_FILENAME);

    toml_datum_t snap = toml_bool_in(root, "keep_snapshots");
    params->keep_snapshots = (snap.ok && snap.u.b) ? true : false;

    toml_datum_t disc = toml_bool_in(root, "discretize_momenta");
    params->discretize_momenta = (disc.ok && disc.u.b) ? true : false;

    toml_datum_t up = toml_bool_in(root, "init_uniform_p_magnitude");
    params->init_uniform_p_magnitude = (up.ok && up.u.b) ? true : false;

    toml_datum_t um = toml_double_in(root, "uniform_p_magnitude");
    if (!um.ok) { toml_datum_t umi = toml_int_in(root, "uniform_p_magnitude");
                  if (umi.ok) { um.ok = 1; um.u.d = (double)umi.u.i; } }
    params->uniform_p_magnitude = um.ok ? um.u.d : 0.0;

    toml_datum_t a = toml_double_in(root, "alfa");
    if (!a.ok) { toml_datum_t ai = toml_int_in(root, "alfa");
                 if (ai.ok) { a.ok = 1; a.u.d = (double)ai.u.i; } }
    params->alfa = a.ok ? a.u.d : 1e-4;

    toml_datum_t pn = toml_double_in(root, "pmin");
    if (!pn.ok) { toml_datum_t pni = toml_int_in(root, "pmin");
                  if (pni.ok) { pn.ok = 1; pn.u.d = (double)pni.u.i; } }
    params->pmin = pn.ok ? pn.u.d : 2e-26;

    toml_datum_t px = toml_double_in(root, "pmax");
    if (!px.ok) { toml_datum_t pxi = toml_int_in(root, "pmax");
                  if (pxi.ok) { px.ok = 1; px.u.d = (double)pxi.u.i; } }
    params->pmax = px.ok ? px.u.d : 3.0e-23;

    toml_free(root);

    printf("%s lee %s  (N_Part:%d, 2D)\t%s escribe %s\tsigma(L) = %le\n",
           params->should_resume ? "sí" : "no", params->inputFilename, params->N_PART,
           params->should_dump   ? "sí" : "no", params->saveFilename,
           params->sigmaL);
}

void io_read_state(const char *filename, SimState *s)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "io_read_state: cannot open %s: %s\n", filename, strerror(errno));
        exit(1);
    }
    int N = s->params.N_PART;
    if (fread(&s->evolution, sizeof(s->evolution), 1, f) != 1 ||
        fread(s->x,  sizeof(double) * N, 1, f) != 1 ||
        fread(s->y,  sizeof(double) * N, 1, f) != 1 ||
        fread(s->px, sizeof(double) * N, 1, f) != 1 ||
        fread(s->py, sizeof(double) * N, 1, f) != 1) {
        fprintf(stderr, "io_read_state: short read from %s\n", filename);
        exit(1);
    }
    fclose(f);
}

void io_write_state(const char *filename, const SimState *s)
{
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        fprintf(stderr, "io_write_state: cannot open %s: %s\n", filename, strerror(errno));
        exit(1);
    }
    int N = s->params.N_PART;
    fwrite(&s->evolution, sizeof(s->evolution), 1, f);
    fwrite(s->x,  sizeof(double) * N, 1, f);
    fwrite(s->y,  sizeof(double) * N, 1, f);
    fwrite(s->px, sizeof(double) * N, 1, f);
    fwrite(s->py, sizeof(double) * N, 1, f);
    fclose(f);
}
