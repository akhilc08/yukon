read("fam.gp");
Hint(K,p,q)={my(A=subst(famA(K),T,p/q)*q^4,B=subst(famB(K),T,p/q)*q^8,c4=16*(A^2-3*B),c6=-64*A^3+288*A*B,c4p=16*(A^2+12*B),c6p=-64*A^3+2304*A*B);
  [log(max(abs(c4)^3,c6^2))/12., log(max(abs(c4p)^3,c6p^2))/12.]};
{
foreach([552,808],K, my(c=vector(12),n=0,mx=0);
  for(i=1,3000, my(q=random(5000)+1,p=random(10000)-5000); if(gcd(p,q)!=1,next);
    my(A=subst(famA(K),T,p/q),B=subst(famB(K),T,p/q),E=ellinit([0,A,0,B,0]),E2=ellinit([0,-2*A,0,A^2-4*B,0]));
    if(E.disc==0,next);
    my(h=Hint(K,p,q), s=min(h[1]-Hnaive(ellminimalmodel(E)), h[2]-Hnaive(ellminimalmodel(E2))));
    n++; c[min(12,floor(s/log(2)+1.5))]++; if(s>mx,mx=s));
  print(K," n=",n," shrink hist in units log2 (bin 1 = 0): ",c," max ",mx));
}
