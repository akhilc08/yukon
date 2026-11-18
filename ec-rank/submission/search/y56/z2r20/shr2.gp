read("fam.gp");
{
foreach([10^5, 10^6, 2*10^6, 10^7, 2^20, 5^9, 2^10*5^5, 2^5*5^8, 2^12*5^6],q,
  my(best=[100,0]); forstep(p=1,q\20,max(1,q\20000), if(gcd(p,q)==1, my(h=Hmin(808,p/q)); if(h<best[1],best=[h,p]))); 
  print(q," ",factor(q)," best H ",best));
}
