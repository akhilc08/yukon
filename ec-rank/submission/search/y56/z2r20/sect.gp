\\ Search sections x = k*P(T) on Y^2 = x^3 + a x^2 + b x over Q(T), P | b squarefree subset of irreducible factors.
read("fam.gp");
\\ formal square-root remainder conditions of polynomial f in T with coeffs in Q(k)
sqconds(f)={my(d=poldegree(f,T),m,lc,g,r);if(d%2,return(0));m=d/2;lc=pollead(f,T);
  g=T^m; \\ monic sqrt of f/lc to order
  my(h=f/lc);for(i=1,m, my(c=polcoef(h-g^2,2*m-i,T)/2); g+=c*T^(m-i));
  r=h-g^2; vector(m,i,numerator(polcoef(r,i-1,T)))};
findsec(a,b)={my(F=factor(b),fl=F[,1],n=#fl,cont=content(b),res=List());
  for(mask=0,2^n-1, my(P=prod(i=1,n,if(bittest(mask,i-1),fl[i],1)),Q=b/P,f,cs,g,ks);
    f=k*(k^2*P+k*a+Q); f=simplify(f);
    if(poldegree(f,T)%2,next);
    cs=sqconds(f); if(cs===0,next);
    g=0; for(i=1,#cs, if(cs[i]!=0, g=gcd(g,cs[i])));
    if(g==0, next); \\ identically square? skip (unlikely)
    ks=[r|r<-nfroots(,g), r!=0];
    for(j=1,#ks, my(kk=ks[j],ff=subst(f,k,kk),xx=kk*P,lcf=pollead(ff,T));
      if(issquare(ff), listput(res,[xx, ff])))); 
  Vec(res)};
\\ isogeny E'->E for y^2=x^3+a'x^2+b'x -> y^2=x^3+ax^2+bx with a'=-2a, b'=a^2-4b: x = y'^2/(4x'^2)
sects(K)={my(a=famA(K),b=famB(K),S1,S2,L=List());
  S1=findsec(a,b); S2=findsec(-2*a,a^2-4*b);
  print(K,": on E ",#S1,"  on E' ",#S2);
  for(i=1,#S1,listput(L,[S1[i][1],S1[i][2]]));
  for(i=1,#S2, my(xp=S2[i][1],yp2=S2[i][2]*0+ (xp^3-2*a*xp^2+(a^2-4*b)*xp)); listput(L,[yp2/(4*xp^2), 0]));
  Vec(L)};
\\ rank of the found sections at a specialisation
chk(K,L,t)={my(a=subst(famA(K),T,t),b=subst(famB(K),T,t),E=ellinit([0,a,0,b,0]),P=List());
  for(i=1,#L, my(x0=subst(L[i][1],T,t),y2=x0^3+a*x0^2+b*x0); if(issquare(y2,&y0), listput(P,[x0,y0]), print("not sq ",i)));
  P=Vec(P); my(M=ellheightmatrix(E,P)); [#P, matrank(M), matdet(M)]};
\\ greedy independent subset of sections, evaluated at fibre t0
indep(K,L,t0)={my(a0=subst(famA(K),T,t0),b0=subst(famB(K),T,t0),E=ellinit([0,a0,0,b0,0]),S=List(),Pt=List());
  for(i=1,#L, my(x0=subst(L[i][1],T,t0),y0); if(!issquare(x0^3+a0*x0^2+b0*x0,&y0),next);
     my(Q=concat(Vec(Pt),[[x0,y0]]),M=ellheightmatrix(E,Q),ev=vecmin(real(qfjacobi(M)[1])));
     if(ev>1e-8, listput(Pt,[x0,y0]); listput(S,L[i])));
  [Vec(S),Vec(Pt),E]};
