\\ survivors stage for EK fibre slice u=U (M = m-parameter giving 9 sections, or 0 for 8).
\\ env: U, M, IN (lines "p q ..."), OUT, HB (height bound), SB (min S(1e5) for ellrank), HITS (dir), TAG
read("ekfam.gp");
U=eval(getenv("U")); M=eval(getenv("M")); IN=getenv("IN"); OUT=getenv("OUT"); HB=eval(getenv("HB")); SB=eval(getenv("SB")); HITS=getenv("HITS"); TAG=getenv("TAG");
AB=famu(U); FA=AB[1]; FB=AB[2]; SECX=secs(U,M);
iso(P,b)=if(P==[0],[0],my(x=P[1],y=P[2]);[y^2/x^2, y*(b-x^2)/x^2]);
models(t)={my(a=subst(FA,T,t),b=subst(FB,T,t),E0=ellinit([0,a,0,b,0]),E1=ellinit([0,-2*a,0,a^2-4*b,0]),P0=List(),v0,v1,M0,M1);
  foreach(SECX,xx, my(x0=subst(xx,T,t),y0); if(x0!=0 && issquare(x0^3+a*x0^2+b*x0,&y0), listput(P0,[x0,y0])));
  P0=Vec(P0); my(P1=apply(P->iso(P,b),P0));
  M0=ellminimalmodel(E0,&v0); M1=ellminimalmodel(E1,&v1);
  [[M0,ellchangepoint(P0,v0),Hn(M0)],[M1,ellchangepoint(P1,v1),Hn(M1)]]};
Sx(E,X)={my(s=0.);forprime(p=3,X,if(E.disc%p,s-=ellap(E,p)*log(p)/p));s};
\\ board/baseline bound per certified rank (Z/2)
bnd(r)=if(r>=20,27.547796, r==19,26.528271, r==18,21.507228, r==17,18.911354, r==16,17.507728, r==15,15.138777, r==14,13.755516, r==13,12.354414, 0);
{
my(L=readstr(IN));
for(i=1,#L, my(w=strsplit(L[i]," "),t,m,j,h,s,r,Lp);
  if(#w<2,next); t=eval(w[1])/eval(w[2]);
  iferr(m=models(t), e, write(OUT,w[1]," ",w[2]," ERR"); next);
  j=if(m[1][3]<=m[2][3],1,2); h=m[j][3];
  if(h>=HB, write(OUT,w[1]," ",w[2]," H ",h); next);
  s=Sx(m[j][1],100000);
  if(s<SB, write(OUT,w[1]," ",w[2]," H ",h," S5 ",s); next);
  iferr(r=ellrank(m[j][1],1,m[j][2]), e, write(OUT,w[1]," ",w[2]," H ",h," S5 ",s," RKERR"); next);
  write(OUT,w[1]," ",w[2]," H ",h," S5 ",s," R ",r[1]," ",r[2]," j ",j);
  if(r[1]>=13 && h < bnd(r[1]),
     Lp=ellsaturation(m[j][1],r[4],100);
     write(HITS"/"TAG"_"w[1]"_"w[2]".txt", [U,M,w[1],w[2],h,r[1],r[2],j], "\n", m[j][1][1..5], "\n", Lp)));
}
quit;
