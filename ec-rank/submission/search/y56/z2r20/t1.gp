read("fam.gp");
for(k=1,2,K=[808,552][k]; print(K," A=",factor(famA(K))," degB=",poldegree(famB(K))); print(" disc fac:",factor(poldisc(x^2+famA(K)*x+famB(K),x)*famB(K))));
\\ check vs board curves
E1=ellinit([0,-1,0,-479337962456030241119664973800901105,127023999849553558599285116689332611801173040541538897]);
F=fibre(808,30/769);F2=fibre2(808,30/769);print(Hnaive(E1)," ",Hnaive(F)," ",Hnaive(F2)," j eq:",E1.j==F.j,E1.j==F2.j);
E2=ellinit([1,1,1,-252979814033828149616633606835402096361712413,-462321266560434395510360229904992327482519890115183022109083267469]);
G=fibre(552,18371/22532);G2=fibre2(552,18371/22532);print(Hnaive(E2)," ",Hnaive(G)," ",Hnaive(G2)," j:",E2.j==G.j,E2.j==G2.j);
