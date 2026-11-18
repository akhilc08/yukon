\\ Torsion Z/2, generic rank 9 (8 over Q(u)[t] for every u; 9 when 5-u^2 is a square).
\\ Source: N. D. Elkies, Z. Klagsbrun, "New rank records for elliptic curves having rational torsion",
\\   ANTS XIV (Open Book Series 4, 2020) 233-250, arXiv:2003.00077, Section 9 eqs (3),(4);
\\   generators Appendix A eqs (7),(8), Gram matrix (9) det 115/16.
\\ K3 surface family from the genus-0 Shimura curve X/w230 (quaternion alg. ramified at {2,5}, Gamma0(23)).
\\ This is the family that produced the Z/2 rank 19 and 20 records (they used u = 2/5 and u = 11/5).
\\ Model: E_u : y^2 = x^3 + 2A x^2 + B x,  B = prod_{i=1..8} B_i.  E_u ~ E_{-u}.
\\ Use: fam(t,u) on Q(u)[t]; famm(t,m) uses u = 2(m^2-m-1)/(m^2+1) so the 9th section is rational.

ekA(t,u) = (u^8-18*u^6+163*u^4-1152*u^2+4096)*t^4 + (3*u^7-35*u^5-120*u^3+1536*u)*t^3 \
         + (u^8-13*u^6+32*u^4-152*u^2+1536)*t^2 + (u^7+3*u^5-156*u^3+672*u)*t + (3*u^6-33*u^4+112*u^2-80);
ekBodd(t,u) = [(u^2+u-8)*t+(-u+2), (u^2-u-8)*t+(u^2+u-10), (u^2-7*u+8)*t+(-u^2+u+2), (u^2+5*u+8)*t+(u^2+3*u+2)];
ekBi(t,u) = { my(o=ekBodd(t,u), e=ekBodd(-t,-u)); [o[1],-e[1],o[2],-e[2],o[3],-e[3],o[4],-e[4]]; }
fam(t,u) = { my(B=ekBi(t,u)); [0, 2*ekA(t,u), 0, prod(i=1,8,B[i]), 0]; }
\\ eight sections over Q(u)[t] (eq. 7)
secx(t,u) = { my(B=ekBi(t,u));
  [-B[1]*B[2]*B[3]*B[6], -B[1]*B[2]*B[4]*B[5], 4*B[1]*B[2]*B[5]*B[6], B[1]*B[3]*B[4]*B[6],
   -B[1]*B[3]*B[4]*B[7], B[1]*B[3]*B[4]*B[8], B[1]*B[3]*B[5]*B[6], -B[1]*B[5]*B[6]*B[7]]; }
um(m) = 2*(m^2-m-1)/(m^2+1);
famm(t,m) = fam(t, um(m));
\\ nine sections (eq. 8) on the m-model
secxm(t,m) = { my(u=um(m), B=ekBi(t,u)); concat(secx(t,u), [-(m-1)^2*B[1]*B[2]*B[3]*B[8]]); }
