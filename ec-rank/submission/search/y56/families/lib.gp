\\ Shared test helpers for the y56 family library.
\\ Each family file defines fam(t[,u]) -> a-invariants [a1,a2,a3,a4,a6] (rational functions OK)
\\ and secx(t[,u]) -> vector of x-coordinates of generic sections (on the SAME model as fam).
\\ famtest(ainvs, xs) specialises nothing: it takes a concrete curve + x list.

naiveH(E) = log(max(abs(E.c4)^3, E.c6^2)) / 12;

\\ lift x to a point on E (rational y), or 0 if x is not the x of a rational point
xtopt(E, x) = {
  my(a1=E.a1,a2=E.a2,a3=E.a3,a4=E.a4,a6=E.a6, b=a1*x+a3, c=x^3+a2*x^2+a4*x+a6, D=b^2+4*c, s);
  if (!issquare(D, &s), return(0));
  [x, (-b+s)/2];
}

\\ test one concrete specialisation. returns [torsion structure, #pts on curve, rank of height matrix (pos def => = #pts), H(minimal model), smallest eigenvalue]
spectest(ai, xs) = {
  my(E, Em, P, Q, M, ev, ok);
  E = ellinit(ai); if (#E == 0 || E.disc == 0, return(["singular"]));
  P = vector(#xs, i, xtopt(E, xs[i]));
  ok = sum(i=1,#P, P[i]!=0);
  Q = select(p -> p != 0, P);
  Em = ellminimalmodel(E);
  if (#Q, M = ellheightmatrix(E, Q); ev = vecmin(real(mateigen(M,1)[1])), ev = 0);
  [elltors(E)[2], ok, #xs, if(#Q, matrank(M), 0), naiveH(Em), ev];
}

\\ height of the minimal model only
minH(ai) = my(E=ellinit(ai)); if(#E==0 || E.disc==0, oo, naiveH(ellminimalmodel(E)));

\\ smallest H over t=a/b with |a|,|b|<=N (b>0), gcd 1 ; returns [Hmin, median-ish average]
hscan(f, N) = {
  my(L=List(), h);
  for(b=1,N, for(a=-N,N, if(gcd(a,b)==1, iferr(h=minH(f(a/b)), e, h=oo); if(h<oo, listput(L,h)))));
  L = vecsort(Vec(L)); [L[1], L[(#L+1)\2], L[#L]];
}
