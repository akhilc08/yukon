\\ Torsion Z/4, generic rank 5 (two one-parameter base changes) and rank 6 (one, over Q(r)),
\\ obtained from the Elkies disc -163 K3 surface E1 (see z4_elkies163_E1.gp) by imposing quadratic sections.
\\ Source: A. Dujella, J. C. Peral, "An elliptic curve over Q(u) with torsion Z/4Z and rank 6",
\\   Rad HAZU Mat. 28 (2024), arXiv:2207.08206v2, Sections 3-4 (all formulas below are from there).
\\ Model (Sec. 3): Y^2 = X^3 + (a^2-8b) X^2 + 16 b^2 X, a=(8t-1)(32t+7), b=8(t+1)(15t-8)(31t-7).
\\ Generic sections on E1 in this model: X1..X4 (Sec. 3).
dpA(t) = 65536*t^4-17472*t^3-10176*t^2+18672*t-3535;
dpfam(t) = [0, dpA(t), 0, 1024*(t+1)^2*(15*t-8)^2*(31*t-7)^2, 0];
dpsec4(t) = [-361*(t+1)*(31*t-7), -4*(t+1)*(15*t-8)*(16*t-7)^2, -16*(t+1)*(8*t+7)^2*(15*t-8), 4*(15*t-8)*(16*t+1)^2*(31*t-7)];
\\ quadratic sections (Sec. 4)
qx1(t) = -64*(1+t)^2*(-4+7*t)*(4+17*t)/(1+4*t)^2;                  \\ rational iff -(7t-4)(17t+4) = square
qx2(t) = 576*(-4+7*t)*(-8+15*t)^2*(-1324+5551*t)/(49*(-39+28*t)^2);  \\ rational iff (7t-4)(5551t-1324) = square
\\ rank-5 family A: t = 4(u^2-1)/(17+7u^2)
tA(u) = 4*(u^2-1)/(17+7*u^2);
famA(u) = dpfam(tA(u));
secA(u) = concat(dpsec4(tA(u)), [qx1(tA(u))]);
\\ rank-5 family B: t = 4(u^2-331)/(7(u^2-793))
tB(u) = 4*(u^2-331)/(7*(u^2-793));
famB(u) = dpfam(tB(u));
secB(u) = concat(dpsec4(tB(u)), [qx2(tB(u))]);
\\ rank-6 family over Q(r): both conditions
tR(r) = 4*(3*r^2-14*r-5390)*(10*r^2-14*r-1617)/(7*(72*r^4-182*r^3-13279*r^2+98098*r+20917512));
famR(r) = dpfam(tR(r));
secR(r) = concat(dpsec4(tR(r)), [qx1(tR(r)), qx2(tR(r))]);
\\ generic families entry points
fam(u) = famA(u);
secx(u) = secA(u);
