read("twsec.gp");
t0=getabstime(); r=twrun(famA(552),famB(552)); print("552E ",#r," ",getabstime()-t0,"ms");
write("tw552E.out",r);
