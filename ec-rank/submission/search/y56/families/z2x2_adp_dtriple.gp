\\ Torsion Z/2 x Z/2, generic rank 5 over Q(m) (curves induced by rational Diophantine triples).
\\ Source: J. Aguirre, A. Dujella, J. C. Peral, "On the rank of elliptic curves coming from rational
\\   Diophantine triples", Rocky Mountain J. Math. 42 (2012) 1759-1776
\\   (preprint https://web.math.pmf.unizg.hr/~duje/pdf/ADP.pdf), Sec. 3.3, n = -7/3: A, B and
\\   sections X1..X5 eq (7), Theorem 2.  Rank-3 two-parameter family over Q(a,n): Sec. 3.1 eq (3).
\\ ADP found rank 10 and 11 fibres in this family (search r/s, |r|,s <= 1000 and 9|s, 10|r restricted search).
\\ Rank-6 subfamilies: Sec. 3.4 (parametrised by points of rank-3 / rank-2 elliptic curves).
\\ Model: y^2 = x^3 + A x^2 + B x.
famA(m) = -2*(167772160000000+1323093196800000*m-32195543040000*m^2-14929920000000*m^3-1863701913600*m^4 \
   +285400350720*m^5+5952139200*m^6-2908045152*m^7+43046721*m^8);
famB(m) = 81*m*(-640+9*m)*(-160+9*m)*(-80+9*m)^2*(32+9*m)*(80+9*m)^2*(-2240+96*m+27*m^2) \
   *(3200+240*m+27*m^2)*(-1600-4320*m+81*m^2)*(4480-720*m+81*m^2);
fam(m) = [0, famA(m), 0, famB(m), 0];
secx(m) = { [27*m*(-640+9*m)*(-80+9*m)^2*(80+9*m)^2*(-2240+96*m+27*m^2),
  3/49*m*(-640+9*m)*(-2240+96*m+27*m^2)*(-108800-11520*m+1539*m^2)^2,
  27*m*(-640+9*m)*(-80+9*m)*(80+9*m)*(-2240+96*m+27*m^2)*(-1600-4320*m+81*m^2),
  27*m*(-640+9*m)*(3200+240*m+27*m^2)*(-1600-4320*m+81*m^2)*(4480-720*m+81*m^2),
  9*(-2240+96*m+27*m^2)*(3200+240*m+27*m^2)*(4480-720*m+81*m^2)*(-1600-4320*m+81*m^2)]; }
\\ two-parameter rank-3 family (Sec 3.1), y^2 = x^3 + A2 x^2 + B2 x
famAn(a,n) = 81+108*a+108*a^2-96*a^3-32*a^4-180*n^2-84*a*n^2-120*a^2*n^2-32*a^3*n^2+118*n^4-28*a*n^4+12*a^2*n^4-20*n^6+4*a*n^6+n^8;
famBn(a,n) = 4*a^2*(9+2*a-n^2)*(3+2*a-4*n+n^2)*(3+2*a+4*n+n^2)*(-3+2*a+3*n^2)*(-9+4*a^2+10*n^2-n^4);
fam2(a,n) = [0, famAn(a,n), 0, famBn(a,n), 0];
secx2(a,n) = { my(p=(3+2*a-4*n+n^2)*(3+2*a+4*n+n^2)); [4*a^2*p, p*(9-6*a-10*n^2-2*a*n^2+n^4)^2/(16*n^2), 2*a*p*(-3+2*a+3*n^2)]; }
