\\ Forced-extra-point search: x = k*P(T) on Y^2=x^3+a x^2+b x, P | b (product of irreducible factors).
\\ f_k(T) = k(k^2 P + k a + b/P).  Point exists at T iff f_k(T) is a square.
\\ For each P: rational k with f_k = s(T) g(T)^2, deg s <= 2 (repeated roots forced via disc in k),
\\ and s(T)=z^2 solvable -> genus-0 curve C: base change T=T(w) with a new section.
\\ deg s in 3..4 -> genus 1 curve C (report; points give sporadic T).
read("fam.gp");
sqfpart(f)={my(F=factor(f),r=1);for(i=1,#F~,if(F[i,2]%2,r*=F[i,1]));r};
conicok(s)={my(d=poldegree(s,T));if(d<=0,return(issquare(s)));if(d==1,return(1));
  my(a=polcoef(s,2,T),b=polcoef(s,1,T),c=polcoef(s,0,T),M=[a,b/2,0;b/2,c,0;0,0,-1]);M=M*denominator(M);type(qfsolve(M))=="t_COL"};
dsec(a,b)={my(F=factor(b),fl=F[,1],n=#fl,res=List());
  for(mask=1,2^n-2, my(P=prod(i=1,n,if(bittest(mask,i-1),fl[i],1)),Q=b/P,f,ks=List(),D);
    f=k*(k^2*P+k*a+Q);
    \\ candidate k: rational roots of disc_T(f) and of the top coefficients (degree drops)
    D=poldisc(f,T); if(D!=0, foreach(nfroots(,numerator(D)),r,listput(ks,r)));
    my(d=poldegree(f,T)); for(j=0,1, my(c=polcoef(f,d-j,T)); if(c!=0 && poldegree(c,k)>0, foreach(nfroots(,c),r,listput(ks,r))));
    foreach(Set(Vec(ks)),kk, if(kk==0,next); my(ff=subst(f,k,kk)); if(ff==0||poldegree(ff,T)<1,next);
      my(cc=content(ff), g=ff/cc, s=sqfpart(g)*core(numerator(cc)*denominator(cc)), ds=poldegree(s,T));
      if(ds<=2 && conicok(s), listput(res,[0,kk*P,s]));
      if(ds>=3 && ds<=4, listput(res,[1,kk*P,s]))));
  Vec(res)};
run(K)={my(a=famA(K),b=famB(K),r1,r2);
  r1=dsec(a,b); r2=dsec(-2*a,a^2-4*b);
  print(K," E: genus0 ",#select(v->v[1]==0,r1)," genus1 ",#select(v->v[1]==1,r1),
        " | E': genus0 ",#select(v->v[1]==0,r2)," genus1 ",#select(v->v[1]==1,r2));
  [r1,r2]};
