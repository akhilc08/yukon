\\ validate ksieve output lines "p q Hlow S1 S2 S3" (P1=200,P2=1000,P3=30000) against PARI. env U, IN
read("ekfam.gp");
U=eval(getenv("U")); IN=getenv("IN");
AB=famu(U); FA=AB[1]; FB=AB[2];
{
my(L=readstr(IN), mx=[0.,0.,0.], hv=0, n=0);
for(i=1,#L, my(w=strsplit(L[i]," ")); if(#w<6,next);
  my(p=eval(w[1]),q=eval(w[2]),a=substvec(FA,[T],[p/q])*q^4,b=substvec(FB,[T],[p/q])*q^8,E=ellinit([0,a,0,b,0]),s=0.,ss=vector(3));
  forprime(l=3,30000, my(ap=ellap(E,l)); s-=ap*log(l)/l; if(l<=200,ss[1]=s); if(l<=1000,ss[2]=s)); ss[3]=s;
  my(E2=ellinit([0,-2*a,0,a^2-4*b,0]),h=min(Hn(ellminimalmodel(E)),Hn(ellminimalmodel(E2))));
  for(k=1,3, mx[k]=max(mx[k],abs(ss[k]-eval(w[3+k]))));
  if(eval(w[3])>h+1e-6, hv++; print("HLOW VIOLATION ",p," ",q," Hlow ",w[3]," H ",h));
  n++);
print("n=",n," maxdiff S1,S2,S3 = ",mx," Hlow violations ",hv);
}
quit;
