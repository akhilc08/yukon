\\ Torsion Z/2 x Z/2, generic rank 4 over Q(a) (n = 7/3 slice of the Diophantine-triple surface),
\\ and rank 5 over Q(w) (quadratic-section base change a = 2(10w-1)(10w+1)/(9(10w^2-1))).
\\ Source: A. Dujella, J. C. Peral, "Elliptic curves induced by Diophantine triples",
\\   RACSAM 113 (2019) 791-806, arXiv:1712.02082, Sec. 2.3 eq (4), Thm 2 (rank 4);
\\   Sec. 2.4 eq (6) (rank 5).  Rank-6 subfamilies over Q(v): Sec. 2.5-2.6.
\\ DP found six rank-11 curves in the rank-4 family (Sec. 2.7).
\\ Model: Y^2 = X^3 + A X^2 + B X.
famA4(a) = -2*(-51200+109440*a+38880*a^2+55404*a^3+6561*a^4);
famB4(a) = 243*a^2*(20+3*a)*(-4+9*a)*(16+9*a)*(80+9*a)*(320+81*a^2);
fam(a) = [0, famA4(a), 0, famB4(a), 0];
secx(a) = { [81*a^2*(-4+9*a)*(80+9*a), 27*a*(20+3*a)*(-4+9*a)*(80+9*a),
  (-4+9*a)*(80+9*a)*(160+171*a)^2/441, 3*(20+3*a)*(-4+9*a)*(320+81*a^2)]; }
famA5(w) = -2*(-169-12020*w^2+678000*w^4-12680000*w^6+80000000*w^8);
famB5(w) = (-1+10*w)^2*(1+10*w)^2*(-1+20*w^2)*(1+80*w^2)*(-31+400*w^2)*(-41+500*w^2)*(9-200*w^2+2000*w^4);
fam5(w) = [0, famA5(w), 0, famB5(w), 0];
secx5(w) = { [(-1+10*w)^2*(1+10*w)^2*(1+80*w^2)*(-41+500*w^2)/9,
  (-1+10*w)*(1+10*w)*(1+80*w^2)*(-31+400*w^2)*(-41+500*w^2)/9,   \\ PDF prints 9*(..); 1/9 is the correct factor
  (1+80*w^2)*(-11+300*w^2)^2*(-41+500*w^2)/49,
  (1+80*w^2)*(-31+400*w^2)*(9-200*w^2+2000*w^4),
  9*(-1+10*w)*(1+10*w)*(-1+20*w^2)*(-41+500*w^2)]; }
