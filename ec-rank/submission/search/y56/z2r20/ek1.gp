read("ekfam.gp");
{
foreach([-4,-5,2/7,-3],m, my(u=um(m),F=fam(T,u),ab=integ(F[2],F[4]),d=ab[3],xs=secxm(T,m));
  print("m=",m," u=",u," d=",d," A ",ab[1]);
  my(hs=List()); foreach([1/3,-2/7,5/11,13/17,-31/41, 101/97],t, my(E=ellinit([0,subst(ab[1],T,t),0,subst(ab[2],T,t),0]),M=ellminimalmodel(E));
     listput(hs,[t,round(Hn(E)-2*log(max(abs(numerator(t)),denominator(t))),&e), round(Hn(M)*100)/100.]));
  print(Vec(hs));
  my(t0=3/7,a0=subst(ab[1],T,t0),b0=subst(ab[2],T,t0),E=ellinit([0,a0,0,b0,0]),P=List());
  foreach(xs,x, my(x0=subst(x,T,t0)*d^2,y0); if(issquare(x0^3+a0*x0^2+b0*x0,&y0), listput(P,[x0,y0])));
  print(" sections on fibre: ",#P," rank ",if(#P,matrank(ellheightmatrix(E,Vec(P))),0)));
}
