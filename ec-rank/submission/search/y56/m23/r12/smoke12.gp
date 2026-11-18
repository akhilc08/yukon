read(Str(ROOT,"/submission/search/y56/m23/r12/r12.gp"));
F=readstr(FIN); F=vector(min(#F,NF),i,eval(strsplit(F[i]," ")));
{foreach(F,b, my(best=[]); for(t=1,TMAX, my(a=concat(b,apply(z->z-t,b))); if(#Set(a)<12,next); my(E=iferr(mesE(a),e,0)); if(E==0,next);
   my(h=Hn(E)); if(h>HMAX,next); best=concat(best,[[Sx(E,3000),t,h]]));
  best=vecsort(best,1,4); for(i=1,min(3,#best), printf("%s t=%d H=%.3f S3k=%.2f\n",b,best[i][2],best[i][3],best[i][1])));}
