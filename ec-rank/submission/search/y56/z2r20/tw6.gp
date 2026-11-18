read("fam.gp");
a=famA(552); b=famB(552); D=54*T^2+18*T+1;
{
for(t=-4,6, my(a0=subst(a*D,T,t),b0=subst(b*D^2,T,t),E=ellinit([0,a0,0,b0,0]),r=ellrank(E,1));
  print("t=",t," rank ",r[1..2]," gens x: ",apply(P->P[1],r[4])));
}
