read("../families/z2_ek_shimura230.gp");
read("fam.gp");
gT(a,b,t)={my(A=subst(a,T,t),B=subst(b,T,t),c4=16*(A^2-3*B),c6=-64*A^3+288*A*B,c4p=16*(A^2+12*B),c6p=-64*A^3+2304*A*B);
  min(log(max(abs(c4)^3,c6^2)),log(max(abs(c4p)^3,c6p^2)))/12};
vol(a,b,X)={my(s=0.,h=0.004);forstep(v=-12,12,h, my(t=sinh(v)); s+=exp(X-gT(a,b,t))*cosh(v)*h); 3/Pi^2*s};
integ(a,b)={my(d=denominator(content([a,b])));[a*d^2,b*d^4]};
{
print("X552 vol ", vol(famA(552),famB(552),27.55));
print("X808 vol ", vol(famA(808),famB(808),27.55));
foreach([2/5,11/5],u, my(F=fam(T,u),ab=integ(F[2],F[4])); print("u=",u," vol ",vol(ab[1],ab[2],27.55)));
L=List();
for(r=-8,8, for(s=1,8, if(gcd(r,s)==1,
  my(m=r/s,u=um(m));
  if(u!=0 && abs(u)!=3 && abs(u)!=4,
    my(F=fam(T,u),ab=integ(F[2],F[4]));
    listput(L,[vol(ab[1],ab[2],27.55),m,u])))));
L=vecsort(Vec(L),1,4); for(i=1,12,print(L[i]));
}
