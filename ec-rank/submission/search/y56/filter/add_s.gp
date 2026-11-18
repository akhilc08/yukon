\\ y56/filter: extra S(X) partial sums up to 1e5 for curves already in a gen_feats.gp data file.
\\ Output per line: fam a1 a2 a3 a4 a6 S(X) for X in XS  (msum convention: 3 <= p <= X, p good, -a_p log p / p)
\\ Usage: echo 'IN="x.txt";OUT="x.s5";NMAX=20000;read("add_s.gp");go();' | gp -q -s 256M
XS = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000];
sx(E) = {
  my(D = E.disc, S = vector(#XS), s = 0., j = 1);
  forprime (p = 3, XS[#XS],
    while (p > XS[j], S[j] = s; j++);
    if (D % p, s -= ellap(E, p) * log(p) / p));
  S[#XS] = s; S;
}
go() = {
  my(L = readstr(IN), f = fileopen(OUT, "w"), n = min(#L, NMAX));
  for (i = 1, n,
    my(t = strsplit(L[i], " "), a = apply(eval, t[4..8]), E = ellinit(a));
    filewrite(f, Str(t[1], " ", strjoin(t[4..8], " "), " ", strjoin(apply(z -> Strprintf("%.4f", z), sx(E)), " "))));
  fileclose(f);
}
