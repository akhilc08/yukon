\\ z2ek: Elkies-Klagsbrun Z/2 family E_u (arXiv:2003.00077 Sec 9 eqs 3-4, App A eqs 7-8) at fixed u.
read("z2_ek_shimura230.gp");
Hn(E)=log(max(abs(E.c4)^3,E.c6^2))/12.;
\\ integral homogeneous model for fixed u: Y^2 = x^3 + A(p,q) x^2 + B(p,q) x, deg 4 / 8 in T
integ(a,b)={my(da=denominator(content(a)),db=denominator(content(b)),d=1,P=factor(da*db)[,1]);
  for(i=1,#P,my(p=P[i],e=max(ceil(valuation(da,p)/2),ceil(valuation(db,p)/4)));d*=p^e);
  my(A=a*d^2,B=b*d^4);
  \\ remove l^2 | A, l^4 | B common content
  my(ca=content(A),cb=content(B),g=gcd(ca^2,cb),F=factor(g)[,1]);
  if(A!=0&&B!=0,for(i=1,#F,my(l=F[i]);while(denominator(content(A/l^2))==1&&denominator(content(B/l^4))==1,A/=l^2;B/=l^4;d/=l)));
  [A,B,d]};
famu(u)=my(F=fam(T,u));integ(F[2],F[4]);
\\ shrink primes: divisors of Res(A,B) incl. leading coeffs, plus 2
shrinkp(A,B)={my(r=polresultant(A,B,T)*pollead(A,T)*pollead(B,T)*polcoef(A,0,T)*polcoef(B,0,T));
  my(P=factor(abs(r))[,1]~);setunion([2],P)};
writefam(u,fn)={my(ab=famu(u),A=ab[1],B=ab[2]);
  write(fn,"A ",strjoin(apply(x->Str(x),Vecrev(A))," "));
  write(fn,"B ",strjoin(apply(x->Str(x),Vecrev(B))," "));
  write(fn,"SP ",strjoin(apply(x->Str(x),shrinkp(A,B))," "))};
\\ section x-coords (scaled to the integral model), rank-9 if m given else rank 8
secs(u,m)={my(d=famu(u)[3],xs=if(m===0,secx(T,u),secxm(T,m)));apply(x->x*d^2,xs)};
