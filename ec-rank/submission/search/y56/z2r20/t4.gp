read("dsec.gp");
{
for(kk=1,2,K=[808,552][kk]; R=run(K);
  for(m=1,2, foreach(R[m],v, if(v[1]==0, print("  ",["E","E'"][m]," g0 x=",factor(v[2])," s=",v[3]))));
  write("dsec"K".out",R));
}
