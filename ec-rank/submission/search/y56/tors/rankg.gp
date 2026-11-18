\\ rankg.gp: ellrank on minimal models of library families at t=p/q, seeded with the family's generic sections.
\\ globals: LIB (list of files to read), FF (closure t->ainvs), SX (closure t->x-coords of sections), TOR (expected torsion),
\\          INF (lines [p,q,...]), OUTF. Output: p q H rlo rhi npts nsec ms ; ".pts" gets [p,q,H,ainvs,points] for rlo>=MINR.
H(E)=log(max(abs(E.c4)^3,E.c6^2))/12;
xtopt(E,x)={my(b=E.a1*x+E.a3,c=x^3+E.a2*x^2+E.a4*x+E.a6,D=b^2+4*c,s);if(!issquare(D,&s),return(0));[x,(-b+s)/2]}
if(type(MINR)=="t_POL", MINR=9);
read("isocert.gp");
iso2(E)=my(R=List());iferr(foreach(isoq(E),z,my(tt=elltors(z[1])[2]);if(tt==[2,2]||tt==[4]||tt==[3],listput(R,[H(z[1]),tt]))),err,);Vec(R);
if(type(EFF)=="t_POL", EFF=0);
if(type(HCUT)=="t_POL", HCUT=1e9);
{
my(L=readvec(INF),i0=1,PF=Str(OUTF,".prog"));
iferr(my(v=readstr(PF)); if(#v, i0=eval(v[#v])+1), e, i0=1);
for(i=i0,#L, my(v=L[i],t0=v[1]/v[2],E0,Em,ch,P,r,tm=getabstime(),ns=0);
  iferr(E0=ellinit(FF(t0)), err, next); if(#E0==0||E0.disc==0, next);
  Em=ellminimalmodel(E0,&ch); if(H(Em)>=HCUT, next);
  P=List(); iferr(foreach(SX(t0),x, my(pt=xtopt(E0,x)); if(pt!=0, listput(P,pt))), err, );
  P=ellchangepoint(Vec(P),ch); ns=#P;
  iferr(r=ellrank(Em,EFF,P), err, r=[-1,-1,0,[]]);
  write(OUTF, v[1]," ",v[2]," ",H(Em)," ",r[1]," ",r[2]," ",#r[4]," ",ns," ",getabstime()-tm," ",elltors(Em)[2]," ",if(r[1]>=9 && (elltors(Em)[2]==[4]||elltors(Em)[2]==[2,2]||elltors(Em)[2]==[3]), iso2(Em), []));
  if(r[1]>=MINR, write(Str(OUTF,".pts"), [v[1],v[2],H(Em),Em[1..5],r[4]])); write(PF,i));
}
