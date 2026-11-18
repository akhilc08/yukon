\\ Family equations (from public ICARM commentary by wgxli: X808 under #4078, X552 under #4061)
\\ y^2 + c xy = x(x+u)(x+v)  ->  Y^2 = x^3 + A x^2 + B x, A = c^2/4+u+v, B = u v
famc(k)=if(k==808, 84*(35397*T^2-221*T+20), 238*(19*T-6)*(5*T-3));
famu(k)=if(k==808, (1139*T+115)*(1269*T-59)*(3776*T+34)*(539*T-5), (12*T-29)*(51*T+22)*(79*T+22)*(482*T+29));
famv(k)=if(k==808, 8*(638*T-23)*(250*T-1)*(7811*T+499), 256*(T-22)*(83*T-7)*(21*T-2)*(7*T+17));
famA(k)=famc(k)^2/4+famu(k)+famv(k);
famB(k)=famu(k)*famv(k);
\\ fibre at T=t, as minimal model
fibre(k,t)=my(A=subst(famA(k),T,t),B=subst(famB(k),T,t));ellminimalmodel(ellinit([0,A,0,B,0]));
Hnaive(E)=my(c4=E.c4,c6=E.c6);log(max(abs(c4)^3,c6^2))/12.;
\\ 2-isogenous: y^2=x^3-2Ax^2+(A^2-4B)x
fibre2(k,t)=my(A=subst(famA(k),T,t),B=subst(famB(k),T,t));ellminimalmodel(ellinit([0,-2*A,0,A^2-4*B,0]));
Hmin(k,t)=min(Hnaive(fibre(k,t)),Hnaive(fibre2(k,t)));
