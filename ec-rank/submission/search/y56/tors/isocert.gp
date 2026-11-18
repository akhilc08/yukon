\\ isoq(E): quotients E/<T> for rational torsion points T of order 2 or 3: list of [Em, iso, ch] (minimal model + map data).
\\ isomap(rec, want): replace record [p,q,H,ainvs,pts] by a quotient curve with torsion structure want ([3],[4],[2,2]).
{
isoq(E)=my(T=elltors(E),K=List(),R=List());
  foreach(T[3],g, for(k=1,ellorder(E,g)-1, my(pt=ellmul(E,g,k),o=ellorder(E,pt)); if(o==2||o==3, listput(K,pt))));
  K=Set(Vec(K));
  foreach(K,k0, my(iso=ellisogeny(E,k0),ec=ellinit(iso[1]),c2,em); em=ellminimalmodel(ec,&c2); listput(R,[em,iso,c2]));
  Vec(R);
}
{
isomap(rec,want)=my(E=ellinit(rec[4]));
  foreach(isoq(E),z, my(em=z[1]); if(elltors(em)[2]==want, my(Q=apply(pt->ellchangepoint(ellisogenyapply(z[2][2],pt),z[3]), rec[5])); return([rec[1],rec[2],H(em),em[1..5],Q])));
  error("no isogenous curve with torsion ",want);
}
