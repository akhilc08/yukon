\\ env U (rational), M (m-parameter or 0), FN (output .fam). Writes A/B coefficient lines and SP shrink primes.
read("ekfam.gp");
U=eval(getenv("U")); M=eval(getenv("M")); FN=getenv("FN");
{
my(ab=famu(U),A=ab[1],B=ab[2]);
if(poldegree(A,T)!=4||poldegree(B,T)!=8||issquarefree(B)==0, print("DEGENERATE ",U); quit);
my(r=abs(polresultant(A,B,T)*pollead(A,T)*pollead(B,T)), F, P=List([2]));
F=iferr(alarm(120,factor(r)),e,factor(r,2^24));
for(i=1,#F~, my(l=F[i,1]); if(l>2 && l<2^31 && isprime(l), listput(P,l)));
write(FN,"A ",strjoin(apply(x->Str(x),Vecrev(A))," "));
write(FN,"B ",strjoin(apply(x->Str(x),Vecrev(B))," "));
write(FN,"SP ",strjoin(apply(x->Str(x),Vec(P))," "));
print("OK ",U," d=",ab[3]," nsp=",#P," logmaxA=",round(log(vecmax(abs(Vec(A)))))," logmaxB=",round(log(vecmax(abs(Vec(B))))));
}
quit;
