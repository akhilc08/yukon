\\ scan the twelve-root family p=f(x)f(x+t): H and S(3000) for t in 1..TMAX for each f in FIN
read(Str(ROOT,"/submission/search/y56/m23/r12/r12.gp"));
F=readstr(FIN); F=vector(#F,i,eval(strsplit(F[i]," ")));
{foreach(F,b, if(b[1]+b[6]<0, next); for(t=1,TMAX, my(a=concat(b,apply(z->z-t,b))); if(#Set(a)<12,next); my(E=iferr(mesE(a),e,0)); if(E==0,next);
   my(h=Hn(E)); if(h>HMAX,next); my(s=Sx(E,1000)); if(s<SMIN,next); write(OUT,strprintf("%s %d %.3f %.2f %.2f",strjoin(apply(z->Str(z),b)," "),t,h,s,Sx(E,3000)))));}
