\\ rankfam.gp: lines [F,p,q,...] -> ellrank of minimal model. Families: 1 X400a, 2 X312 (Z/2xZ/2), 3 z3k3 (Z/3).
\\ Output: F p q H rlo rhi npts ms [z4 isogenous: [H, ainvs]...]   globals INF OUTF
fA(F,p,q)=if(F==1, q*(121*q-99*p)*(26*p-23*q)*(75*p-44*q), -q*(29*p+2*q)*(400*p+49*q)*(221*p+86*q));
fB(F,p,q)=if(F==1, p*(21*q-74*p)*(25*p-36*q)*(p+77*q), -p*(8*p+53*q)*(25*p+88*q)*(392*p+221*q));
Ef(F,p,q)={if(F<=2, my(A=fA(F,p,q),B=fB(F,p,q),g=core(gcd(A,B),1)[2]); return(ellminimalmodel(ellinit([0,(A+B)/g^2,0,A*B/g^4,0]))));
  if(F==3, my(a1=7*p^2+309*q^2,a3=-91*(p^2-q^2)*(p^2-25*q^2)*(p^2-121*q^2)); return(ellminimalmodel(ellinit([a1,0,a3,0,0]))));}
H(E)=log(max(abs(E.c4)^3,E.c6^2))/12;
z4(E)={my(M=ellisomat(E,2),L=M[1],R=List());for(i=1,#L,my(e=ellminimalmodel(ellinit(L[i][1])));if(elltors(e)[2]==[4],listput(R,[H(e),e[1..5]])));Vec(R)}
{
my(L=readvec(INF));
for(i=1,#L, my(v=L[i],E=Ef(v[1],v[2],v[3]),r,t=getabstime(),zz=[]);
  iferr(r=ellrank(E), err, r=[-1,-1,0,[]]);
  if(v[1]<=2 && r[1]>=8, zz=z4(E));
  write(OUTF, v[1]," ",v[2]," ",v[3]," ",H(E)," ",r[1]," ",r[2]," ",#r[4]," ",getabstime()-t," ",zz);
  if(r[1]>=9, write(Str(OUTF,".pts"), [v[1],v[2],v[3],H(E),E[1..5],r[4]])));
}
