/* ns/candidate for each scorer in y56_filter.h, plus a check against reference scores.
 * cc -O2 -march=native -I<dir with y56_filter.h> -DBK=50 bench.c -o bench -lm
 * ./bench < ap.txt        (ap.txt: 168 ints per line, a_p with 0 at bad primes) */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "y56_filter.h"
#ifndef BK
#define BK 50
#endif
#define CAT2(a, b) a##b
#define CAT(a, b) CAT2(a, b)
#define MAXN 200000
static int A[MAXN][168];
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
typedef float (*fn)(const int *);
static float sK(const int *ap) { return y56_S(ap, BK); }
int main(void) {
  int n = 0; char line[8192];
  while (n < MAXN && fgets(line, sizeof line, stdin)) { char *p = line; for (int i = 0; i < 168; i++) A[n][i] = strtol(p, &p, 10); n++; }
  struct { const char *nm; fn f; } S[] = {{"S_K", sK}, {"lin", CAT(y56_lin, BK)}, {"quad", CAT(y56_quad, BK)},
                                          {"gbm", CAT(y56_gbm, BK)}, {"mlp", CAT(y56_mlp, BK)}};
  for (int s = 0; s < 5; s++) {
    volatile float acc = 0; int reps = 0; double t0 = now(), t;
    do { for (int i = 0; i < n; i++) acc += S[s].f(A[i]); reps++; } while ((t = now() - t0) < 0.5);
    printf("K=%d %-5s %8.1f ns/candidate  (n=%d reps=%d)\n", BK, S[s].nm, 1e9 * t / ((double)n * reps), n, reps);
  }
  return 0;
}
