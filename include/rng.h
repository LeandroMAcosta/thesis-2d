#ifndef RNG_H
#define RNG_H

#include <stdint.h>

#ifdef __CUDACC__
#define RNG_HD __host__ __device__
#else
#define RNG_HD
#endif

/* xorshift32 — period 2^32 - 1, 32 bits of state. Good enough for short
 * runs; for long runs (>~1 G samples per thread) consider upgrading to
 * xoshiro256** or PCG. Callable from both host and device. */
static inline RNG_HD double rng_xorshift32(uint32_t *state)
{
    uint32_t x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    return (double)x / (double)UINT32_MAX;
}

/* Mix a 64-bit base seed and a small integer (thread id, batch index, ...)
 * into a non-zero 32-bit xorshift seed. The multipliers are the Knuth and
 * Fibonacci hashing constants — they decorrelate adjacent inputs that
 * would otherwise produce nearly-identical xorshift streams. */
static inline RNG_HD uint32_t rng_seed_mix(uint64_t base, uint32_t aux)
{
    uint64_t mix = base * 2654435761ull + (uint64_t)aux * 0x9E3779B1ull;
    uint32_t s = (uint32_t)(mix ^ (mix >> 32));
    return s == 0 ? 1u : s;
}

#endif
