\\ bisection search: x = k*P*L^2 (P | b, deg P <= 2, L = T - m or linear factor), k free.
\\ h_k(T) = k(k^2 P L^4 + k a L^2 + b/P); need h = s*g^2 with deg s <= 2 (and z^2=s solvable).
read("dsec.gp");
lsec(a,b,Ls,maxdP)={my(F=factor(b),fl=F[,1],n=#fl,res=List());
  for(mask=0,2^n-1, my(P=prod(i=1,n,if(bittest(mask,i-1),fl[i],1)),Q=b/P); if(poldegree(P,T)>maxdP,next);
    foreach(Ls,L, my(h=k*(k^2*P*L^4+k*a*L^2+Q),D=poldisc(h,T),ks=List());
      if(D==0,next);
      foreach(nfroots(,numerator(D)),r,listput(ks,r));
      my(d=poldegree(h,T)); for(j=0,2, my(c=polcoef(h,d-j,T)); if(c!=0 && poldegree(c,k)>0, foreach(nfroots(,c),r,listput(ks,r))));
      foreach(Set(Vec(ks)),kk, if(kk==0,next); my(ff=subst(h,k,kk)); if(ff==0||poldegree(ff,T)<1,next);
        my(cc=content(ff), g=ff/cc, s=sqfpart(g)*core(numerator(cc)*denominator(cc)), ds=poldegree(s,T));
        if(ds<=2 && conicok(s), listput(res,[ds,kk*P*L^2,s]); print("  HIT ds=",ds," x=",factor(kk*P*L^2)," s=",s)))));
  Vec(res)};
