\\ twists E^{cD} with x = k*P*L^2 (L = pT+q, |p|<=PL, |q|<=QL), P | D^2 B squarefree-subset, deg P <= DPM
\\ env FAMK MODEL PL QL DPM OUT ; shards over D index: env NSH, ISH
read("twsec.gp"); read("dsec.gp");
K=eval(getenv("FAMK")); MD=eval(getenv("MODEL")); PL=eval(getenv("PL")); QL=eval(getenv("QL")); DPM=eval(getenv("DPM")); OUT=getenv("OUT");
NSH=eval(getenv("NSH")); ISH=eval(getenv("ISH"));
{
a=famA(K); b=famB(K); if(MD, [a,b]=[-2*a,a^2-4*b]);
my(L=concat(factor(b)[,1]~, factor(a^2-4*b)[,1]~),lin=select(f->poldegree(f,T)==1,L),Ds=List());
listput(Ds,1); foreach(lin,f,listput(Ds,f)); for(i=1,#lin,for(j=i+1,#lin,listput(Ds,lin[i]*lin[j])));
Ls=List(); for(p=0,PL,for(q=-QL,QL,if(gcd(p,q)==1 && (p>0||q==1), listput(Ls,p*T+q)))); Ls=Vec(Ls);
for(id=1,#Ds, if((id-1)%NSH!=ISH, next); my(D=Ds[id],a2=a*D,b2=b*D^2,F=factor(b2),fl=F[,1],n=#fl,found=0);
  for(mask=0,2^n-1, my(P=prod(i=1,n,if(bittest(mask,i-1),fl[i],1)),Q=b2/P); if(poldegree(P,T)>DPM,next);
    foreach(Ls,LL, my(h=k*(k^2*P*LL^4+k*a2*LL^2+Q), cs, g);
      if(poldegree(h,T)%2, next);
      cs=sqconds(h); if(cs===0,next); g=0; for(i=1,#cs, if(cs[i]!=0, g=gcd(g,cs[i]))); if(g==0,next);
      foreach(nfroots(,g),kk, if(kk==0,next); my(ff=subst(h,k,kk),lc=pollead(ff));
        if(issquare(ff/lc), my(c=core(numerator(lc)*denominator(lc))); if(conicok(c*D),
           if(D!=1 || c!=1, found++; write(OUT,"HIT D=",D," c=",c," x=",kk*P*LL^2)))))));
  write(OUT,"Ddone ",id," ",D," found ",found));
}
