// z2ek copy of y56/z2r20/ksieve.c; only change: env SH2 sets the 2-adic height allowance (default log 2)
// ksieve: height-first Mestre-Nagao sieve over fibres of a Z/2 elliptic K3 family
//   Y^2 = x^3 + A(T) x^2 + B(T) x,  T = p/q,  A deg<=4, B deg<=8 (homogenised to degrees 4 and 8).
// Order per fibre: (1) cheap height bound of the minimal model (double c4,c6 of the integral model,
// minus exact odd-prime shrinkage at the family's finite set of possible shrink primes, minus an
// optimistic 2-adic allowance), min over the 2-isogenous pair; reject if >= X.
// (2) S(P1) via per-prime lookup tables a_l[(p:q) mod l]; (3) S(P2); (4) S(P3).  S = sum -a_l log l / l.
// Build: cc -O3 -mcpu=native -fopenmp -o ksieve ksieve.c -lm
// Usage: ksieve FAMFILE X qmin qmax pmax T1 T2 T3 [P1 P2 P3] > out
//   FAMFILE: line "A a0 a1 a2 a3 a4", "B b0 .. b8", "SP l1 l2 ..." (odd shrink primes; 2 handled optimistically)
// Output lines: p q Hlow S1 S2 S3
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

static double Ac[5], Bc[9];
static __int128 Ai[5], Bi[9];  // integer coefficients (|c| < 2^127)
static double SHMAX = 8.0;
static uint8_t *ZA[64];  // ZA[j][r]=1 iff A(r:1)==0 mod sp[j] (r=l: point at infinity)
static int nsp, sp[64], spk[64];
static unsigned long long spM[64];
static int np, pr[10000];
static float *W[10000];
static int P1 = 200, P2 = 1000, P3 = 10000, i1, i2, i3;

static long long pw(long long a, long long e, long long m) {
  long long r = 1; a %= m; if (a < 0) a += m;
  while (e) { if (e & 1) r = (__int128)r * a % m; a = (__int128)a * a % m; e >>= 1; }
  return r;
}
static long long hom_mod(const __int128 *c, int d, long long p, long long q, unsigned long long M) {
  // sum c_i p^i q^(d-i) mod M
  __int128 mm = M, pp = ((p % (long long)M) + (long long)M) % (long long)M, qq = ((q % (long long)M) + (long long)M) % (long long)M;
  __int128 s = 0, pi = 1;
  __int128 qp[9]; qp[0] = 1; for (int i = 1; i <= d; i++) qp[i] = qp[i - 1] * qq % mm;
  for (int i = 0; i <= d; i++) {
    __int128 ci = ((c[i] % mm) + mm) % mm;
    s = (s + ci * pi % mm * qp[d - i]) % mm;
    pi = pi * pp % mm;
  }
  return (long long)s;
}
static int vval(long long r, int l, int k) { if (r == 0) return k; int v = 0; while (r % l == 0) { r /= l; v++; } return v; }

static void build(void) {
  int maxp = P3;
  char *sv = calloc(maxp + 1, 1);
  np = 0;
  for (int i = 2; i <= maxp; i++) if (!sv[i]) { if (i > 2) pr[np++] = i; for (int j = 2 * i; j <= maxp; j += i) sv[j] = 1; }
  free(sv);
  i1 = i2 = i3 = 0;
  for (int k = 0; k < np; k++) { if (pr[k] <= P1) i1 = k + 1; if (pr[k] <= P2) i2 = k + 1; if (pr[k] <= P3) i3 = k + 1; }
#pragma omp parallel for schedule(dynamic, 1)
  for (int k = 0; k < np; k++) {
    int l = pr[k];
    int8_t *chi = malloc(l);
    for (int x = 0; x < l; x++) chi[x] = -1;
    chi[0] = 0;
    for (long long x = 1; x < l; x++) chi[x * x % l] = 1;
    float *w = malloc(sizeof(float) * (l + 1));
    double lw = log((double)l) / l;
    for (int r = 0; r <= l; r++) {
      long long a, b;
      if (r < l) { a = hom_mod(Ai, 4, r, 1, l); b = hom_mod(Bi, 8, r, 1, l); }
      else { a = hom_mod(Ai, 4, 1, 0, l); b = hom_mod(Bi, 8, 1, 0, l); }
      // chi(x(x^2+ax+b)) = chi(x) chi(g(x)), g(x+1) = g(x) + 2x + 1 + a  (incremental, no division)
      long long s = 0;
      int g = (int)b, dg = (int)((1 + a) % l);  // g(0) = b, g(1)-g(0) = 1 + a
      for (int x = 0; x < l; x++) {
        s += chi[x] * chi[g];
        g += dg; if (g >= l) g -= l;
        dg += 2; if (dg >= l) dg -= l;
      }
      w[r] = (float)(s * lw);  // -a_l log l / l, a_l = -s
    }
    W[k] = w;
    free(chi);
  }
}

static inline double hmodel(double A, double B, int iso) {
  double c4, c6;
  if (!iso) { c4 = 16 * (A * A - 3 * B); c6 = -64 * A * A * A + 288 * A * B; }
  else { c4 = 16 * (A * A + 12 * B); c6 = -64 * A * A * A + 2304 * A * B; }
  double m = fmax(3 * log(fabs(c4) + 1e-300), 2 * log(fabs(c6) + 1e-300));
  return m / 12;
}

static __int128 p128(const char *w) { int neg = 0; __int128 v = 0; if (*w == '-') { neg = 1; w++; } while (*w >= '0' && *w <= '9') v = v * 10 + (*w++ - '0'); return neg ? -v : v; }

static int NG; static double *GT, *GG, *PMN, *PMX;
static int cmpg(const void *a, const void *b) { double x = GG[*(const int *)a], y = GG[*(const int *)b]; return x < y ? -1 : x > y; }
static void build_grid(void) {
  NG = 4000001; GT = malloc(sizeof(double) * NG); GG = malloc(sizeof(double) * NG);
  for (int k = 0; k < NG; k++) {
    double u = -22.0 + 44.0 * k / (NG - 1), t = sinh(u) * 1e-3;
    double A = 0, B = 0, pi = 1;
    for (int i = 0; i <= 8; i++) { if (i <= 4) A += Ac[i] * pi; B += Bc[i] * pi; pi *= t; }
    GT[k] = t; GG[k] = fmin(hmodel(A, B, 0), hmodel(A, B, 1));
  }
  int *ix = malloc(sizeof(int) * NG); for (int k = 0; k < NG; k++) ix[k] = k;
  qsort(ix, NG, sizeof(int), cmpg);
  double *g2 = malloc(sizeof(double) * NG); PMN = malloc(sizeof(double) * NG); PMX = malloc(sizeof(double) * NG);
  double mn = 1e300, mx = -1e300;
  for (int k = 0; k < NG; k++) { double t = GT[ix[k]]; if (t < mn) mn = t; if (t > mx) mx = t; PMN[k] = mn; PMX[k] = mx; g2[k] = GG[ix[k]]; }
  free(GG); GG = g2; free(ix);
}
// hull of {T : g(T) < lam}; returns 0 if empty
static int hull(double lam, double *lo, double *hi) {
  int a = 0, b = NG; while (a < b) { int m = (a + b) / 2; if (GG[m] < lam) a = m + 1; else b = m; }
  if (a == 0) return 0;
  double t1 = PMN[a - 1], t2 = PMX[a - 1];
  // widen by grid spacing
  *lo = t1 - fabs(t1) * 0.002 - 1e-9; *hi = t2 + fabs(t2) * 0.002 + 1e-9; return 1;
}

int main(int argc, char **argv) {
  if (argc < 9) { fprintf(stderr, "usage\n"); return 1; }
  FILE *f = fopen(argv[1], "r");
  char tag[8];
  while (fscanf(f, "%7s", tag) == 1) {
    if (!strcmp(tag, "A")) for (int i = 0; i < 5; i++) { char w[64]; fscanf(f, "%63s", w); Ai[i] = p128(w); Ac[i] = (double)Ai[i]; }
    else if (!strcmp(tag, "B")) for (int i = 0; i < 9; i++) { char w[64]; fscanf(f, "%63s", w); Bi[i] = p128(w); Bc[i] = (double)Bi[i]; }
    else if (!strcmp(tag, "SP")) { int l; while (fscanf(f, "%d", &l) == 1) { if (l >= 2) sp[nsp++] = l; } break; }
  }
  fclose(f);
  if (getenv("SHMAX")) SHMAX = atof(getenv("SHMAX"));
  for (int j = 0; j < nsp; j++) { int l = sp[j]; ZA[j] = calloc(l + 1, 1); for (int r = 0; r <= l; r++) ZA[j][r] = (r < l ? hom_mod(Ai, 4, r, 1, l) : hom_mod(Ai, 4, 1, 0, l)) == 0; }
  for (int j = 0; j < nsp; j++) { unsigned long long M = 1; int k = 0; while (M < (1ULL << 62) / sp[j]) { M *= sp[j]; k++; } spM[j] = M; spk[j] = k; }
  double X = atof(argv[2]);
  long long qmin = atoll(argv[3]), qmax = atoll(argv[4]), pmax = atoll(argv[5]);
  double T1 = atof(argv[6]), T2 = atof(argv[7]), T3 = atof(argv[8]);
  if (argc > 11) { P1 = atoi(argv[9]); P2 = atoi(argv[10]); P3 = atoi(argv[11]); }
  double t0 = omp_get_wtime();
  {
    // table cache: TABLES=<file> loads if present, else builds and saves
    const char *tf = getenv("TABLES");
    FILE *tfp = tf ? fopen(tf, "rb") : NULL;
    if (tfp) {
      int hdr[2]; fread(hdr, sizeof(int), 2, tfp);
      if (hdr[0] != P3) { fprintf(stderr, "table P3 mismatch\n"); return 1; }
      np = hdr[1]; fread(pr, sizeof(int), np, tfp);
      for (int k = 0; k < np; k++) { W[k] = malloc(sizeof(float) * (pr[k] + 1)); fread(W[k], sizeof(float), pr[k] + 1, tfp); }
      fclose(tfp);
      i1 = i2 = i3 = 0;
      for (int k = 0; k < np; k++) { if (pr[k] <= P1) i1 = k + 1; if (pr[k] <= P2) i2 = k + 1; if (pr[k] <= P3) i3 = k + 1; }
    } else {
      build();
      if (tf) { FILE *o = fopen(tf, "wb"); int hdr[2] = {P3, np}; fwrite(hdr, sizeof(int), 2, o); fwrite(pr, sizeof(int), np, o);
        for (int k = 0; k < np; k++) fwrite(W[k], sizeof(float), pr[k] + 1, o); fclose(o); }
    }
  }
  build_grid();
  fprintf(stderr, "tables built: %d primes <= %d in %.1fs\n", np, P3, omp_get_wtime() - t0);
  long long cnt[5] = {0, 0, 0, 0, 0};
  double L2 = getenv("SH2") ? atof(getenv("SH2")) : log(2.0);  /* 2-adic allowance (H units); z2ek: 1.5 */
#pragma omp parallel
  {
    long long lc[5] = {0, 0, 0, 0, 0};
    int *iq = malloc(sizeof(int) * np);
    char buf[256];
#pragma omp for schedule(dynamic, 1)
    for (long long q = qmin; q <= qmax; q++) {
      for (int k = 0; k < i3; k++) { int l = pr[k]; iq[k] = (q % l) ? (int)pw(q, l - 2, l) : -1; }
      double qd = (double)q;
      long long iqs[64]; for (int j = 0; j < nsp; j++) { long long qm = q % sp[j]; iqs[j] = qm ? pw(qm, sp[j] - 2, sp[j]) : 0; }
      double qp[9]; qp[0] = 1; for (int i = 1; i <= 8; i++) qp[i] = qp[i - 1] * qd;
      double tlo, thi;
      if (!hull(X + SHMAX + L2 - 2 * log(qd), &tlo, &thi)) continue;
      long long plo = (long long)floor(tlo * qd) - 1, phi = (long long)ceil(thi * qd) + 1;
      if (plo < -pmax) plo = -pmax; if (phi > pmax) phi = pmax;
      for (long long p = plo; p <= phi; p++) {
        double pd = (double)p, A = 0, B = 0, pi = 1;
        for (int i = 0; i <= 8; i++) { if (i <= 4) A += Ac[i] * pi * qp[4 - i]; B += Bc[i] * pi * qp[8 - i]; pi *= pd; }
        if (B == 0) continue;
        double h0 = hmodel(A, B, 0), h1 = hmodel(A, B, 1);
        if (fmin(h0, h1) - L2 - SHMAX >= X) continue;
        long long g = p < 0 ? -p : p, h = q; while (h) { long long t = g % h; g = h; h = t; }
        if (g != 1) continue;
        lc[0]++;
        double sh0 = 0, sh1 = 0;
        for (int j = 0; j < nsp; j++) {
          int l = sp[j];
          long long qm = q % l, pm = ((p % l) + l) % l;
          int r = qm ? (int)(pm * iqs[j] % l) : l;
          if (!ZA[j][r]) continue;
          unsigned long long M = spM[j]; int k = spk[j];
          long long a = hom_mod(Ai, 4, p, q, M), b = hom_mod(Bi, 8, p, q, M);
          int va = vval(a, l, k), vb = vval(b, l, k);
          int e0 = va / 2 < vb / 4 ? va / 2 : vb / 4;
          long long bb = (long long)((((__int128)a * a - 4 * (__int128)b) % (__int128)M + M) % M);
          int vbb = vval(bb, l, k);
          int e1 = va / 2 < vbb / 4 ? va / 2 : vbb / 4;
          sh0 += e0 * log((double)l); sh1 += e1 * log((double)l);
        }
        double Hlow = fmin(h0 - sh0, h1 - sh1) - L2;
        if (Hlow >= X) continue;
        lc[1]++;
        long long pa = p < 0 ? -p : p;
        double s = 0;
        int k = 0;
        for (; k < i1; k++) { int l = pr[k]; int r = iq[k] < 0 ? l : (int)((((p % l) + l) % l) * iq[k] % l); s += W[k][r]; }
        if (s < T1) continue;
        double s1 = s;
        lc[2]++;
        for (; k < i2; k++) { int l = pr[k]; int r = iq[k] < 0 ? l : (int)((((p % l) + l) % l) * iq[k] % l); s += W[k][r]; }
        if (s < T2) continue;
        double s2 = s;
        lc[3]++;
        for (; k < i3; k++) { int l = pr[k]; int r = iq[k] < 0 ? l : (int)((((p % l) + l) % l) * iq[k] % l); s += W[k][r]; }
        if (s < T3) continue;
        lc[4]++;
        (void)pa;
        snprintf(buf, sizeof buf, "%lld %lld %.4f %.3f %.3f %.3f\n", p, q, Hlow, s1, s2, s);
#pragma omp critical
        { fputs(buf, stdout); fflush(stdout); }
      }
    }
#pragma omp critical
    for (int i = 0; i < 5; i++) cnt[i] += lc[i];
    free(iq);
  }
  double el = omp_get_wtime() - t0;
  fprintf(stderr, "pairs %lld  Hpass %lld  S1pass %lld  S2pass %lld  S3pass %lld  time %.1fs threads %d\n", cnt[0], cnt[1], cnt[2], cnt[3], cnt[4], el, omp_get_max_threads());
  return 0;
}
