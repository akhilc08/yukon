\\ env FAMK, PMAXL, QMAXL, MODEL (0: E, 1: E'), OUT
read("lsec.gp");
K=eval(getenv("FAMK")); PL=eval(getenv("PMAXL")); QL=eval(getenv("QMAXL")); MD=eval(getenv("MODEL")); OUT=getenv("OUT");
{
a=famA(K); b=famB(K); if(MD, [a,b]=[-2*a,a^2-4*b]);
Ls=List(); for(p=1,PL,for(q=-QL,QL,if(gcd(p,q)==1,listput(Ls,p*T+q)))); listput(Ls,1+0*T);
r=lsec(a,b,Vec(Ls),2);
write(OUT,r);
foreach(r,v,if(v[1]>0, write(OUT, "BISECTION ", v)));
}
