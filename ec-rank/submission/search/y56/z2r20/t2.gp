read("sect.gp");
{
for(kk=1,2,K=[808,552][kk]; L=sects(K); print("count ",#L); for(i=1,#L,print("  x=",factor(L[i][1])));
 t0=3/17; a0=subst(famA(K),T,t0); b0=subst(famB(K),T,t0);
 Pt=vector(#L,i,my(x0=subst(L[i][1],T,t0),y0);if(!issquare(x0^3+a0*x0^2+b0*x0,&y0),error("nsq"));[x0,y0]);
 M=ellheightmatrix(ellinit([0,a0,0,b0,0]),Pt); print("rank of found: ",matrank(round(M*10^20)));
);
}
