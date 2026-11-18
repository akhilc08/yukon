// nagao: estimate the Q(T)-rank of quadratic twists E^{cD} of a K3 family Y^2 = x^3 + A x^2 + B x
// by Nagao sums  r(D,c) ~ (1/X) sum_{l<=X} -(log l / l) sum_{t in P^1(F_l)} chi(c D(t)) a_l(E_t).
// A positive r(D,c) flags a quadratic base change T = T(s) (parametrising c D(T) = z^2) that gains rank.
// D(p,q) = a p^2 + b p q + c q^2 (a=0: one ramification point at infinity).
// Build: cc -O3 -o nagao nagao.c -lm.  Usage: nagao FAMFILE X [thr] < list  (lines: a b c)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static __int128 Ai[5], Bi[9];
static __int128 p128(const char *w) { int neg = 0; __int128 v = 0; if (*w == '-') { neg = 1; w++; } while (*w >= '0' && *w <= '9') v = v * 10 + (*w++ - '0'); return neg ? -v : v; }
static long long md(__int128 c, long long l) { long long r = (long long)(c % l); return r < 0 ? r + l : r; }
static long long hom(const __int128 *c, int d, long long p, long long q, long long l) {
  long long s = 0, pi = 1, qp[9]; qp[0] = 1; for (int i = 1; i <= d; i++) qp[i] = qp[i - 1] * q % l;
  for (int i = 0; i <= d; i++) { s = (s + md(c[i], l) * pi % l * qp[d - i]) % l; pi = pi * p % l; }
  return s;
}
static long long pw(long long a, long long e, long long m) { long long r = 1; a %= m; if (a < 0) a += m; while (e) { if (e & 1) r = r * a % m; a = a * a % m; e >>= 1; } return r; }
int np, pr[5000]; int *AP[5000]; signed char *CH[5000];
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); char tag[8];
  while (fscanf(f, "%7s", tag) == 1) {
    if (!strcmp(tag, "A")) for (int i = 0; i < 5; i++) { char w[64]; fscanf(f, "%63s", w); Ai[i] = p128(w); }
    else if (!strcmp(tag, "B")) for (int i = 0; i < 9; i++) { char w[64]; fscanf(f, "%63s", w); Bi[i] = p128(w); }
    else break;
  }
  fclose(f);
  int X = atoi(argv[2]);
  for (int l = 5; l <= X; l++) { int ok = 1; for (int d = 2; d * d <= l; d++) if (l % d == 0) { ok = 0; break; } if (ok) pr[np++] = l; }
  for (int k = 0; k < np; k++) {
    int l = pr[k]; signed char *chi = malloc(l); for (int x = 0; x < l; x++) chi[x] = -1; chi[0] = 0; for (long long x = 1; x < l; x++) chi[x * x % l] = 1;
    CH[k] = chi; AP[k] = malloc(sizeof(int) * (l + 1));
    for (int r = 0; r <= l; r++) {
      long long a = r < l ? hom(Ai, 4, r, 1, l) : hom(Ai, 4, 1, 0, l), b = r < l ? hom(Bi, 8, r, 1, l) : hom(Bi, 8, 1, 0, l);
      long long s = 0; int g = (int)b, dg = (int)((1 + a) % l);
      for (int x = 0; x < l; x++) { s += chi[x] * chi[g]; g += dg; if (g >= l) g -= l; dg += 2; if (dg >= l) dg -= l; }
      AP[k][r] = (int)(-s);
    }
  }
  double SL = 0; for (int k = 0; k < np; k++) SL += log((double)pr[k]);
  // untwisted estimate
  { double s = 0; for (int k = 0; k < np; k++) { int l = pr[k]; long long t = 0; for (int r = 0; r <= l; r++) t += AP[k][r]; s += -(double)t / l * log((double)l); } printf("BASE %.3f\n", s / SL); }
  int cs[] = {1, -1, 2, -2, 3, -3, 5, -5, 6, -6, 7, -7, 10, -10, 11, -11, 13, -13, 14, -14, 15, -15, 17, -17, 19, -19, 21, -21, 22, -22, 23, -23, 26, -26, 29, -29, 30, -30, 31, -31, 33, -33, 35, -35, 37, -37, 39, -39, 41, -41, 42, -42, 43, -43, 47, -47, 51, -51, 53, -53, 59, -59, 61, -61, 67, -67, 71, -71, 73, -73, 79, -79, 83, -83, 89, -89, 97, -97};
  int ncs = sizeof cs / sizeof cs[0];
  long long qa, qb, qc;
  double thr = argc > 3 ? atof(argv[3]) : -1e9;
  double *v = malloc(sizeof(double) * np);
  // input lines: a b c  -> D(p,q) = a p^2 + b p q + c q^2
  while (scanf("%lld %lld %lld", &qa, &qb, &qc) == 3) {
    for (int k = 0; k < np; k++) {
      int l = pr[k]; long long t = 0, A2 = md(qa, l), B2 = md(qb, l), C2 = md(qc, l);
      for (int r = 0; r < l; r++) { long long Dv = ((A2 * r + B2) % l * r + C2) % l; t += CH[k][Dv] * AP[k][r]; }
      t += CH[k][A2] * AP[k][l];
      v[k] = -(double)t / l * log((double)l);
    }
    double best = -1e9; int bc = 0;
    for (int j = 0; j < ncs; j++) {
      double s = 0; for (int k = 0; k < np; k++) { int l = pr[k]; long long cm = ((cs[j] % l) + l) % l; s += CH[k][cm] * v[k]; }
      s /= SL; if (s > best) { best = s; bc = cs[j]; }
    }
    if (best >= thr) { printf("%lld %lld %lld best_c %d r %.3f\n", qa, qb, qc, bc, best); fflush(stdout); }
  }
  return 0;
}
