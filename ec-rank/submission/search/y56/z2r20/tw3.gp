read("twsec.gp");
print("sanity D=1: ",#twfind(famA(552),famB(552),1));
t0=getabstime(); a=famA(552); b=famB(552); r=twrun(-2*a,a^2-4*b); print("552E' ",#r," ",getabstime()-t0,"ms"); write("tw552Ep.out",r);
t0=getabstime(); r=twrun(famA(808),famB(808)); print("808E ",#r," ",getabstime()-t0,"ms"); write("tw808E.out",r);
a=famA(808); b=famB(808); r=twrun(-2*a,a^2-4*b); print("808E' ",#r," ",getabstime()-t0,"ms"); write("tw808Ep.out",r);
