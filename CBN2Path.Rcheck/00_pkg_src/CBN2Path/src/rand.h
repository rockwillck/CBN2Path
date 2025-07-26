#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

// PCG state and increment
static uint64_t pcg_state = 0x853c49e6748fea9bULL; // default initial state
static uint64_t pcg_inc = 0xda3e39cb94b95bdbULL;   // default stream

// PCG random number generator (like rand)
int pcg_rand(void) {
  uint64_t oldstate = pcg_state;
  // Advance internal state
  pcg_state = oldstate * 6364136223846793005ULL + pcg_inc;
  // Calculate output function (XSH RR), uses old state
  uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
  uint32_t rot = oldstate >> 59u;
  uint32_t result = (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
  return (int)(result % (RAND_MAX));
}

// Seed PCG (like srand)
void pcg_srand(uint64_t seed) {
  pcg_state = 0U;
  pcg_inc = (seed << 1u) | 1u;  // must be odd
  pcg_rand();                   // advance state
  pcg_state += seed;
  pcg_rand();
}
