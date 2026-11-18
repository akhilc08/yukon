\\ more rank-8 u slices: u=r/s, scale = max(log max|A|/2, log max|B|/4) of the integral model; also rank-9 u=um(m) for small m
read("ekfam.gp");
sc(u)=my(ab=famu(u),A=ab[1],B=ab[2]); if(poldegree(A,T)!=4||poldegree(B,T)!=8,return(99)); max(log(vecmax(abs(Vec(A))))/2, log(vecmax(abs(Vec(B))))/4);
{
my(L=List());
for(s=1,16,for(r=1,48, if(gcd(r,s)!=1,next); my(u=r/s); if(u==3||u==4,next); my(c=sc(u)); if(c<=10.6, listput(L,[c,u,0]))));
for(s=1,12,for(r=-24,24, if(gcd(r,s)!=1,next); my(m=r/s,u=um(m)); if(u==0||abs(u)==3||abs(u)==4,next); my(c=sc(u)); if(c<=11.5, listput(L,[c,u,m]))));
L=vecsort(Vec(L),1);
foreach(L,x, my(u=x[2],m=x[3]); print("u",numerator(u),"_",denominator(u),"m",numerator(m),"_",denominator(m)," ",u," ",m," ",x[1]));
}
quit;
