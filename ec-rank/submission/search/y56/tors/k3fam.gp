\\ k3fam.gp: generalised Z/3 K3 families y^2 + a1 x y + a3 y = x^3, a1 = t^2 + beta, a3 = g*prod (t^2-k_i^2),
\\ (z3k3 on the board is k=(1,5,11), g=-13/49 after scaling a1 to be monic).
\\ A 3-factor section x=c*prod(t-e_i k_i) exists iff beta = e3/s - 4 g s^2 (s=sum e_i k_i, e3=prod e_i k_i), c=-2 g s.
\\ 4-factor sections x=c*P4 exist iff G=(a1 c+g a3/(g P4))^2+4c^3 P4 is a square: two polynomial conditions in (c,g).
lin(k)=[t-k[1],t+k[1],t-k[2],t+k[2],t-k[3],t+k[3]];
betaof(k,e,g)={my(v=vector(3,i,e[i]*k[i]),s=vecsum(v),e3=prod(i=1,3,v[i])); if(s==0,return(oo)); e3/s-4*g*s^2}
\\ rational g making 4-subset S give a section, with beta from pattern e
g4cands(k,e,S)={my(L=lin(k),P4=prod(i=1,4,L[S[i]]),A3,b,G,Gi,E1,E2,R,out=List());
  b=betaof(k,e,'g); A3='g*prod(i=1,3,t^2-k[i]^2);
  G=((t^2+b)*'c+A3/P4)^2+4*'c^3*P4; G=subst(G,t,'t);
  Gi=vector(5,i,polcoef(G,i-1,'t));
  E1=8*Gi[5]^2*Gi[2]-(4*Gi[4]*Gi[5]*Gi[3]-Gi[4]^3); E2=64*Gi[5]^3*Gi[1]-(4*Gi[5]*Gi[3]-Gi[4]^2)^2;
  E1=numerator(E1); E2=numerator(E2);
  R=polresultant(E1,E2,'c); if(R==0,return([]));
  foreach(nfroots(,R),r, if(r!=0, listput(out,r))); Vec(out)}
\\ check a g: count all subsets (size 1..5) giving sections via small-c search is expensive; instead count 4-subsets solved
\\ fixed: strip common factor before resultant
g4c(k,e,S)={my(L=lin(k),P4=prod(i=1,4,L[S[i]]),A3,b,G,Gi,E1,E2,D,R,out=List());
  b=betaof(k,e,'g); if(b==oo,return([])); A3='g*prod(i=1,3,t^2-k[i]^2);
  G=((t^2+b)*'c+A3/P4)^2+4*'c^3*P4; Gi=vector(5,i,polcoef(G,i-1,'t));
  E1=numerator(8*Gi[5]^2*Gi[2]-(4*Gi[4]*Gi[5]*Gi[3]-Gi[4]^3)); E2=numerator(64*Gi[5]^3*Gi[1]-(4*Gi[5]*Gi[3]-Gi[4]^2)^2);
  if(E1==0||E2==0,return([]));
  D=gcd(E1,E2); R=polresultant(E1/D,E2/D,'c); if(R==0,return([]));
  foreach(nfroots(,R),r, if(r!=0, listput(out,r))); Vec(out)}
PATS=[[1,1,1],[1,1,-1],[1,-1,1],[-1,1,1]];
\\ tally gammas for k triple: returns list of [count4, count3, g, e]
k3scan(k)={my(res=List());
  foreach(PATS,e, my(M=Map(),gs=List());
    forsubset([6,4],S, foreach(g4c(k,e,Vec(S)),g, my(c=0); if(mapisdefined(M,g,&c), mapput(M,g,c+1), mapput(M,g,1); listput(gs,g))));
    foreach(gs,g, my(c4=mapget(M,g), b=betaof(k,e,g), c3=sum(j=1,4,betaof(k,PATS[j],g)==b)); listput(res,[c4,c3,g,e])));
  vecsort(Vec(res),(x,y)->sign((y[1]+2*y[2])-(x[1]+2*x[2])))}
