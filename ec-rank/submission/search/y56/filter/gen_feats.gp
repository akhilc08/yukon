\\ y56/filter: dataset generator for the learned-filter experiment (PARI/GP 2.17.4).
\\ Samples family specialisations, computes sieve features and a rank label, one line per curve:
\\   fam band param a1 a2 a3 a4 a6 H logD nbad tors rootno S100 S316 S1000 S3162 S10000 r R tfeat_ms tlab_ms ap(168, comma-joined)
\\ S(X) = sum_{3<=p<=X, p good} -a_p log p / p   (same convention as calib/msum; larger = better).
\\ Label: torsion/rand families -> ellrank(E,0) lower bound r and upper bound R (alarm LABT s; -1 on timeout).
\\        mestre -> #independent points from hyperellratpoints(r, MESB) (mestre-own/mestre.gp) as r; R = -1.
\\ Usage (see README.txt): echo 'FAM="Z6";BAND=1;N=1000;SEED=1;OUT="x.txt";read("gen_feats.gp");go();' | gp -q -s 1G
read(concat(GENDIR, "/../../calib/families.gp"));
read(concat(GENDIR, "/../../mestre-own/mestre.gp"));
if (type(LABT) != "t_INT", LABT = 120);
if (type(MESB) != "t_INT", MESB = 5000);
PR = primes([2, 1000]);

\\ band bounds: [lo, hi] on max(|u|, v) (Tate), max(|A|,|B|) (Z2Z2), coefficient scale (rand), root range (mestre)
bandb(fam, band) = {
  if (fam == "Z2Z2", return(if(band == 1, [1, 500], [500, 5000])));
  if (fam == "rand", return(if(band == 1, [1, 1], [1, 10])));
  if (fam == "mestre", return(if(band == 1, [20, 30], [31, 50])));
  if(band == 1, [1, 600], [600, 6000]);
}

\\ returns [minimal ellinit, param string, extra] or 0
sample(fam, band) = {
  my(b = bandb(fam, band), lo = b[1], hi = b[2], E, u, v, A, B);
  if (fam == "Z2Z2",
    A = random(2*hi + 1) - hi; B = random(2*hi + 1) - hi;
    if (A == 0 || B == 0 || A == B || max(abs(A), abs(B)) < lo, return(0));
    iferr(E = ellinit(ellminimalmodel(ellinit([0, A + B, 0, A*B, 0]))), e, return(0));
    return([E, Str(A, "/", B), 0]));
  if (fam == "rand",
    my(a = [random(2), random(3) - 1, random(2), random(2001*hi) - 1000*hi, random(200001*hi) - 100000*hi]);
    if (band == 2 && abs(a[4]) <= 1000 && abs(a[5]) <= 100000, return(0));
    iferr(E = ellinit(a), e, return(0)); if (#E == 0 || E.disc == 0, return(0));
    E = ellinit(ellminimalmodel(E));
    return([E, strjoin(apply(z -> Str(z), a), "/"), 0]));
  if (fam == "mestre",
    my(s = Set()); while (#s < 10, s = setunion(s, [random(2*hi + 1) - hi]));
    my(a = Vec(s)); if (a[10] - a[1] < 2*lo, return(0));
    a = a - vector(10, i, a[1]);  \\ translate so that min root = 0 (translation gives an isomorphic curve)
    my(F); iferr(F = mesPts(a, 1), e, return(0));
    return([ellinit(F[1]), strjoin(apply(z -> Str(z), a), "/"), a]));
  v = random(hi) + 1; u = random(2*hi + 1) - hi;
  if (gcd(u, v) != 1 || max(abs(u), v) < lo, return(0));
  E = famcurve(fam, u/v); if (E === 0, return(0));
  [ellinit(E), Str(u, "/", v), 0];
}

feats(E) = {
  my(an = ellan(E, 10000), D = E.disc, S = vector(5), X = [100, 316, 1000, 3162, 10000], s = 0., j = 1, nb = 0);
  forprime (p = 3, 10000,
    while (p > X[j], S[j] = s; j++);
    if (D % p, s += -an[p] * log(p) / p));
  S[5] = s;
  foreach (PR, p, if (D % p == 0, nb++));
  [log(abs(D)), nb, S, vector(#PR, i, an[PR[i]])];
}

go() = {
  my(f = fileopen(OUT, "a"), seen = Map(), n = 0, tries = 0);
  setrand(SEED);
  while (n < N && tries < 50*N, tries++;
    my(smp = sample(FAM, BAND)); if (smp === 0, next);
    my(E = smp[1], key = [E.a1, E.a2, E.a3, E.a4, E.a6]);
    if (mapisdefined(seen, key), next); mapput(seen, key, 1);
    my(t0 = getabstime(), F = feats(E), tors = elltors(E)[1], rn = 0);
    iferr(rn = alarm(5, ellrootno(E)), e, rn = 0);
    my(tf = getabstime() - t0, r = -1, R = -1);
    t0 = getabstime();
    if (FAM == "mestre",
      iferr(my(P = mesPts(smp[3], MESB)); r = #indep(ellinit(P[1]), P[2]), e, r = -1),
      iferr(my(L = alarm(LABT, ellrank(E, 0))); r = L[1]; R = L[2], e, r = -1));
    my(tl = getabstime() - t0);
    my(H = log(max(abs(E.c4)^3, E.c6^2)) / 12);
    filewrite(f, Str(FAM, " ", BAND, " ", smp[2], " ", E.a1, " ", E.a2, " ", E.a3, " ", E.a4, " ", E.a6, " ",
      Strprintf("%.5f %.3f", H, F[1]), " ", F[2], " ", tors, " ", rn, " ",
      strjoin(apply(z -> Strprintf("%.4f", z), F[3]), " "), " ", r, " ", R, " ", tf, " ", tl, " ",
      strjoin(apply(z -> Str(z), F[4]), ",")));
    fileflush(f); n++);
  fileclose(f);
}
