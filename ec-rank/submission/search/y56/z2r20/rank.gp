\\ survivors stage: exact minimal-model height of both 2-isogenous models, S(10^5), ellrank seeded with sections.
\\ usage: gp -q rank.gp with K, IN, OUT, HB, MODE defined via -D? we use read of a params file "par.gp"
read("fam.gp");
SECX = if(K==808, read("sec808.gp"), read("sec552.gp"));
SECX = if(K==808, SEC808, SEC552);
iso(P,b)=if(P==[0],[0],my(x=P[1],y=P[2]);[y^2/x^2, y*(b-x^2)/x^2]);
models(t)={my(a=subst(famA(K),T,t),b=subst(famB(K),T,t),E0=ellinit([0,a,0,b,0]),E1=ellinit([0,-2*a,0,a^2-4*b,0]),P0=List(),v0,v1,M0,M1);
  foreach(SECX,xx, my(x0=subst(xx,T,t),y0); if(x0!=0 && issquare(x0^3+a*x0^2+b*x0,&y0), listput(P0,[x0,y0])));
  P0=Vec(P0); my(P1=apply(P->iso(P,b),P0));
  M0=ellminimalmodel(E0,&v0); M1=ellminimalmodel(E1,&v1);
  [[M0,ellchangepoint(P0,v0),Hnaive(M0)],[M1,ellchangepoint(P1,v1),Hnaive(M1)]]};
Sx(E,X)={my(s=0.);forprime(p=3,X,if(E.disc%p,s-=ellap(E,p)*log(p)/p));s};
stageA(t)={my(m=models(t),j=if(m[1][3]<=m[2][3],1,2));[m[j][3],j,Sx(m[j][1],100000)]};
stageB(t,j,effort)={my(m=models(t)[j],E=m[1],P=m[2],r,L);
  r=ellrank(E,effort,P); L=r[4];
  L=ellsaturation(E,L,100);
  [r[1],r[2],L,E,m[3]]};
