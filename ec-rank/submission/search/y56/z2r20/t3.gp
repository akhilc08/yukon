read("sect.gp");
{
for(kk=1,2,K=[808,552][kk]; L=sects(K); 
  R=indep(K,L,3/17); print(K," indep at 3/17: ",#R[1]); R2=indep(K,R[1],-41/113); print("  at -41/113: ",#R2[1]);
  write("sec"K".gp", "SEC"K"=",apply(s->s[1],R[1]),";");
);
}
