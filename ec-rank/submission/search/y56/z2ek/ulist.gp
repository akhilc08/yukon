\\ prints candidate u slices: "tag U M" ; rank-9 slices u=um(m) first, then rank-8 slices u=r/s
um(m) = 2*(m^2-m-1)/(m^2+1);
{
my(S=Map(),out=List());
for(s=1,8, for(r=-8,8, if(gcd(r,s)!=1,next); my(m=r/s,u=abs(um(m)));
  if(u==0||u==3||u==4||mapisdefined(S,u),next); mapput(S,u,m); listput(out,[u,m])));
for(s=1,8, for(r=1,16, if(gcd(r,s)!=1,next); my(u=r/s);
  if(u==3||u==4||mapisdefined(S,u),next); mapput(S,u,0); listput(out,[u,0])));
foreach(out,x, my(u=x[1],m=x[2],mm=if(um(m)==u||m==0,m,0));
  \\ if um(m)=-u use E_u ~ E_{-u}: keep signed u matching m
  my(uu=if(m!=0,um(m),u));
  print("u",numerator(uu),"_",denominator(uu),"m",numerator(m),"_",denominator(m)," ",uu," ",m));
}
quit;
