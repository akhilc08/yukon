read("twsec.gp"); K=eval(getenv("FAMK")); MD=eval(getenv("MODEL")); OUT=getenv("OUT");
a=famA(K); b=famB(K); if(MD, [a,b]=[-2*a,a^2-4*b]);
r=twrun(a,b); write(OUT,r); print("DONE ",#r);
