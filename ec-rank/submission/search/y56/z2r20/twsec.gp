\\ Forced extra section via quadratic base change: for D(T) = product of two bad-fibre linear forms
\\ (or one, ramified at infinity), search sections of the twist E^{cD}: y^2 = x^3 + cDA x^2 + c^2D^2B x.
\\ Ansatz x = k*P (P | D^2 B, with multiplicity), f = x^3+DAx^2+D^2Bx must be (const c)*square in Q[T].
\\ Then on the base change T = T(s) parametrising z^2 = c D(T) (needs a rational point), E gets a new section
\\ x_E = x / (c D(T)) -> rank r(E)+r(E^cD) family over Q(s).
read("sect.gp");
conicok(s)={my(d=poldegree(s,T));if(d<=0,return(issquare(s)));if(d==1,return(1));
  my(a=polcoef(s,2,T),b=polcoef(s,1,T),c=polcoef(s,0,T),M=[a,b/2,0;b/2,c,0;0,0,-1]);M=M*denominator(M);type(qfsolve(M))=="t_COL"};
twfind(a,b,D)={my(a2=a*D,b2=b*D^2,F=factor(b2),fl=F[,1],ml=F[,2],n=#fl,res=List());
  forvec(e=vector(n,i,[0,ml[i]]), my(P=prod(i=1,n,fl[i]^e[i]),Q=b2/P,f,cs,g,ks);
    f=k*(k^2*P+k*a2+Q); if(poldegree(f,T)%2,next);
    cs=sqconds(f); if(cs===0,next);
    g=0; for(i=1,#cs, if(cs[i]!=0, g=gcd(g,cs[i])));
    if(g==0, next);
    ks=[r|r<-nfroots(,g), r!=0];
    foreach(ks,kk, my(ff=subst(f,k,kk),lc=pollead(ff));
      if(issquare(ff/lc), my(c=core(numerator(lc)*denominator(lc))); if(conicok(c*D), listput(res,[D,c,kk*P]))))); 
  Vec(res)};
twrun(a,b)={my(L=concat(factor(b)[,1]~, factor(a^2-4*b)[,1]~),lin=select(f->poldegree(f,T)==1,L),Ds=List(),out=List());
  foreach(lin,f,listput(Ds,f)); for(i=1,#lin,for(j=i+1,#lin,listput(Ds,lin[i]*lin[j])));
  foreach(Ds,D, my(r=twfind(a,b,D)); if(#r, print("D=",factor(D)," -> ",#r," sections, c=",Set(apply(v->v[2],r))); foreach(r,v,listput(out,v))));
  Vec(out)};
