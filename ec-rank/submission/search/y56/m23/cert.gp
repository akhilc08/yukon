\\ m23 certify/write: from roots a (ten-root Mestre), gather points (quartic points bound B + ellrank),
\\ keep an independent set (verifier-style: positive-definite height matrix), optionally saturate, write .ec.
\\ usage: echo 'ROOT="..";a=[...];B=20000;FN="..";HDR="# ...";read("cert.gp")' | gp -q -D parisizemax=8G
read(Str(ROOT,"/submission/search/mestre-own/mestre.gp"));
read(Str(ROOT,"/submission/search/y56/m23/lib.gp"));
if(type(EXTRA)!="t_VEC", EXTRA=[]);
indep2(M,P)={my(S=[]); P=vecsort(P,z->ellheight(M,z)); for(i=1,#P, my(T=concat(S,[P[i]])); my(G=ellheightmatrix(M,T)); if(vecmin(qfjacobi(G)[1])>1e-20, S=T)); S};
{my(F=mesPts(a,B),M=F[1],P=concat(F[2],EXTRA)); P=indep2(M,P); if(type(TARGET)!="t_INT", TARGET=99); if(type(TMAXR)!="t_INT", TMAXR=300);
 my(R=0); if(#P<TARGET, R=iferr(alarm(TMAXR,ellrank(M,1,P)),e,0); if(type(R)!="t_VEC", R=0)); if(R!=0, P=indep2(M,concat(P,R[4])));
 printf("rank bounds %s, indep %d\n",if(R!=0,R[1..2],"ellrank skipped/failed"),#P);
 my(Ps=iferr(ellsaturation(M,P,50),e,P)); P=indep2(M,Ps);
 my(G=ellheightmatrix(M,P)); printf("H=%.4f r=%d minev=%.3g\n",Hn(M),#P,vecmin(qfjacobi(G)[1]));
 my(s=HDR); s=Str(s,"\nE ",M.a1," ",M.a2," ",M.a3," ",M.a4," ",M.a6);
 for(i=1,#P, s=Str(s,"\nP ",fmtq(P[i][1])," ",fmtq(P[i][2]))); s=Str(s,"\n#CLAIM rank=",#P,"\n");
 my(f=fileopen(FN,"w")); filewrite(f,s); fileclose(f);}
quit;
