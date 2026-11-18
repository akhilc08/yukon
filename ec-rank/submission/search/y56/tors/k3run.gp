read("k3fam.gp");
Ek(k,b,g,t0)=ellinit([t0^2+b,0,g*prod(i=1,3,t0^2-k[i]^2),0,0]);
rk(k,b,g)={my(m=99); foreach([2/7,3/13,5/17],t0, my(E=Ek(k,b,g,t0)); if(E.disc==0,next); m=min(m,ellrank(ellminimalmodel(E))[1])); m}
{
for(k3=K3LO,KMAX, for(k2=2,k3-1, for(k1=1,k2-1, if(gcd([k1,k2,k3])>1,next); my(k=[k1,k2,k3],r=k3scan(k),seen=List());
  for(i=1,#r, my(x=r[i]); if(x[1]+2*x[2]<12, break); my(b=betaof(k,x[4],x[3])); if(setsearch(Set(seen),[b,x[3]]),next); listput(seen,[b,x[3]]);
     my(E=Ek(k,b,x[3],2/7)); if(E.disc==0, next); my(tt=elltors(ellminimalmodel(E))[2]); if(tt!=[3], write(OUT,k," ",x," tors ",tt); next);
     write(OUT,k," ",x," beta=",b," rank>= ",rk(k,b,x[3]))))));
}
