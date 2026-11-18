read("../families/z2_ek_shimura230.gp");
\\ integral model of y^2 = x^3 + a x^2 + b x over Q[T]: x -> x/d^2
integ(a,b)={my(da=denominator(content(a)),db=denominator(content(b)),d=1,P=factor(da*db)[,1]);
  for(i=1,#P,my(p=P[i],e=max(ceil(valuation(da,p)/2),ceil(valuation(db,p)/4)));d*=p^e);
  my(A=a*d^2,B=b*d^4,c=content([A,B]));
  \\ remove common square-scaling if possible
  [A,B,d]};
Hn(E)=log(max(abs(E.c4)^3,E.c6^2))/12.;
