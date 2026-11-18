read("fam.gp");
{
for(kk=1,2,K=[808,552][kk];A=famA(K);B=famB(K);
 print("FAM ",K);print("A ",Vecrev(Vec(A),5));print("B ",Vecrev(Vec(B),9));
 F=factor(B)[,1]; S=[2,3]; S=concat(S,factor(pollead(A),10^6)[,1]~);
 for(i=1,#F, S=concat(S, factor(abs(polresultant(A,F[i])),10^6)[,1]~));
 Bp=A^2-4*B; F=factor(Bp)[,1]; for(i=1,#F, S=concat(S, factor(abs(polresultant(A,F[i])),10^6)[,1]~));
 S=vecsort(S,,8);print("SHRINKPRIMES ",[s|s<-S,s<10^7]);
);
}
