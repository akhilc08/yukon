read("lsec.gp");
{
for(kk=1,2,K=[552,808][kk]; a=famA(K); b=famB(K);
  Lf=concat(factor(b)[,1]~, select(f->poldegree(f,T)==1, factor(a^2-4*b)[,1]~));
  Ls=List(); foreach(Lf,f,listput(Ls,f)); for(p=1,6,for(q=-12,12,if(gcd(p,q)==1,listput(Ls,p*T+q)))); Ls=Vec(Ls);
  print(K," #L=",#Ls);
  my(t0=getabstime()); r1=lsec(a,b,Ls,2); print(K," E hits ",#r1," ",getabstime()-t0,"ms");
  r2=lsec(-2*a,a^2-4*b,Ls,2); print(K," E' hits ",#r2);
  write("lsec"K".out",[r1,r2]));
}
