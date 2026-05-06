#include "state.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void state_alloc(SimState *s)
{
    int N = s->params.N_PART;
    int B = s->params.BINS;

    s->x   = malloc(sizeof(double) * N);
    s->y   = malloc(sizeof(double) * N);
    s->px  = malloc(sizeof(double) * N);
    s->py  = malloc(sizeof(double) * N);

    s->hx  = malloc(sizeof(int) * (2 * B + 4));
    s->hy  = malloc(sizeof(int) * (2 * B + 4));
    s->gpx = malloc(sizeof(int) * (2 * B));
    s->gpy = malloc(sizeof(int) * (2 * B));

    s->DxE = malloc(sizeof(double) * (2 * B + 4));
    s->DpE = malloc(sizeof(double) * (2 * B));

    s->h_xy = malloc(sizeof(int) * HXY_BINS * HXY_BINS);

    if (!s->x || !s->y || !s->px || !s->py ||
        !s->hx || !s->hy || !s->gpx || !s->gpy ||
        !s->DxE || !s->DpE || !s->h_xy) {
        fprintf(stderr, "state_alloc: malloc failed\n");
        exit(1);
    }

    state_zero_histograms(s);
    s->evolution = 0u;
}

void state_free(SimState *s)
{
    free(s->x);    s->x   = NULL;
    free(s->y);    s->y   = NULL;
    free(s->px);   s->px  = NULL;
    free(s->py);   s->py  = NULL;
    free(s->hx);   s->hx  = NULL;
    free(s->hy);   s->hy  = NULL;
    free(s->gpx);  s->gpx = NULL;
    free(s->gpy);  s->gpy = NULL;
    free(s->DxE);  s->DxE = NULL;
    free(s->DpE);  s->DpE = NULL;
    free(s->h_xy); s->h_xy = NULL;
}

void state_zero_histograms(SimState *s)
{
    int B = s->params.BINS;
    memset(s->hx,   0, sizeof(int) * (2 * B + 4));
    memset(s->hy,   0, sizeof(int) * (2 * B + 4));
    memset(s->gpx,  0, sizeof(int) * (2 * B));
    memset(s->gpy,  0, sizeof(int) * (2 * B));
    memset(s->h_xy, 0, sizeof(int) * HXY_BINS * HXY_BINS);
}
