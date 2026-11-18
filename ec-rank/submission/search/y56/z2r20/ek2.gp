read("ekfam.gp");
{
foreach([-1,-2,2,-3,3/2,-3/2,4,1/4,-1/4,3/4,-4/3,5/2,-2/5,2/3],m, my(u=um(m),F=fam(T,u),ab=integ(F[2],F[4]),d=ab[3],xs=secxm(T,m));
  my(fb=factor(ab[2])); 
  my(hs=List()); foreach([1/3,-2/7,5/11,13/17,-31/41,101/97],t, my(E=ellinit([0,subst(ab[1],T,t),0,subst(ab[2],T,t),0])); if(E.disc==0, listput(hs,"sing"), my(M=ellminimalmodel(E)); listput(hs, round(Hn(M)*10)/10.)));
  my(t0=3/7,a0=subst(ab[1],T,t0),b0=subst(ab[2],T,t0),E=ellinit([0,a0,0,b0,0]),P=List(),rk=0);
  foreach(xs,x, my(x0=subst(x,T,t0)*d^2,y0); if(x0!=0 && issquare(x0^3+a0*x0^2+b0*x0,&y0), listput(P,[x0,y0])));
  if(#P, my(M=ellheightmatrix(E,Vec(P))); rk=matrank(M));
  print("m=",m," u=",u," #Bfac=",#fb~," H:",Vec(hs)," secrank ",rk));
}
