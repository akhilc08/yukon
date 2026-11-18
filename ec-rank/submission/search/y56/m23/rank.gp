\\ m23 rank stage (PARI/GP 2.17.x). Input lines "a1..a10 H S199 S1000 S3000" (m23f) or "b1..b6 v u S199 S1000 S3000" (r12/m12f: roots v*b, v*b-u).
\\ usage: echo 'ROOT="..";IN="..";OUT="..";B=3000;EFF=1;S2MIN=118;read("rank.gp")' | gp -q -D parisizemax=8G
\\ stage 1: minimal model, H, S(20000) for every line; stage 2 (S20k>=S2MIN): quartic points
\\ (hyperellratpoints, bound B) mapped to E, then ellrank with those points; lb = #independent points found.
read(Str(ROOT,"/submission/search/mestre-own/mestre.gp"));
read(Str(ROOT,"/submission/search/y56/m23/lib.gp"));
if(type(B)!="t_INT", B=3000);
if(type(EFF)!="t_INT", EFF=1);
if(type(S2MIN)=="t_POL", S2MIN=118);
if(type(HMAX)=="t_POL", HMAX=17.88);
fmtv(v)=strjoin(apply(z->Str(z),v)," ");
\\ stage 1 for the whole batch; stage 2 (unchanged threshold S20k>=S2MIN) runs in order of
\\ sc2 = S20k - 0.103*log|disc| + 1.49*nbad(p<1000) (re-ordering only); both scores are logged.
{L=readstr(IN); my(J=List());
 foreach(L,l, my(v=eval(strsplit(l," ")),a,s3); if(#v==14, a=v[1..10]; s3=v[14], if(#v==11, my(b=v[7]*v[1..6]); a=concat(b,apply(z->z-v[8],b)); s3=v[11], next));
  my(E=iferr(mesE(a),e,0)); if(E==0, next); my(h=Hn(E)); if(h>HMAX, write(OUT, strprintf("%s H=%.4f skip",fmtv(a),h)); next);
  my(S2=Sx(E,20000), ld=log(abs(E.disc)), nb=0); forprime(p=2,1000, if(E.disc%p==0, nb++)); my(sc2=S2-0.103*ld+1.49*nb);
  if(S2>=S2MIN, listput(J,[sc2,a,h,s3,S2]),
    write(OUT, strprintf("%s H=%.4f S3k=%.2f S20k=%.2f sc2=%.2f lb=-1 ub=-1 ms=0",fmtv(a),h,s3,S2,sc2))));
 J=vecsort(Vec(J),1,4);
 foreach(J,z, my(sc2=z[1],a=z[2],h=z[3],s3=z[4],S2=z[5],lb=-1,ub=-1,t0=getabstime());
    my(F=iferr(mesPts(a,B),e,0)); if(F!=0, my(M=F[1], I=indep(M,F[2]), R);
      R=iferr(alarm(600,ellrank(M,EFF,I)),e,0); if(type(R)!="t_VEC", R=0); if(R!=0, I=indep(M,concat(I,R[4])); ub=R[2]); lb=#I;
      if(lb>=16, write(Str(OUT,".pts"), Str(fmtv(a),"|",M[1..5],"|",I))));
    write(OUT, strprintf("%s H=%.4f S3k=%.2f S20k=%.2f sc2=%.2f lb=%d ub=%d ms=%d",fmtv(a),h,s3,S2,sc2,lb,ub,getabstime()-t0))) };
quit;
