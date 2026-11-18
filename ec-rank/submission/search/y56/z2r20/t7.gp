read("sect.gp"); read("sec552.gp"); read("sec808.gp");
R=read("lsec552.out");
{
K=552; a=famA(K); b=famB(K);
L=apply(x->[x,0],SEC552);
foreach(R[1],v, L=concat(L,[[v[2],0]]));
foreach(R[2],v, my(xp=v[2]); L=concat(L,[[ (xp^3-2*a*xp^2+(a^2-4*b)*xp)/(4*xp^2),0]]));
J=indep(K,L,3/17); print("552 indep: ",#J[1]); J2=indep(K,J[1],-41/113); print(" check ",#J2[1]);
write("sec552.gp","SEC552=",apply(s->s[1],J[1]),";");
}
