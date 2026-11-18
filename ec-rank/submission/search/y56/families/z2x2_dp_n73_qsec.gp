\\ Torsion Z/2 x Z/2: rank-5 and rank-6 subfamilies of the Dujella-Peral rank-4 family (z2x2_dp_dtriple_n73.gp)
\\ obtained HERE (y56) by automated quadratic-section search (qsec.gp): x = c*prod(bad-fibre factors),
\\ x(x^2+Ax+B) = square * s(t) with deg s = 2 and z^2=s(t) solvable -> base change t=t_i(w) adds 1 to the rank.
\\ Q6 (10(9a-20)(9a-2)=square) is DP's own rank-5 step (arXiv:1712.02082 Sec 2.4); the other five were not listed there.
\\ Requires z2x2_dp_dtriple_n73.gp to be read first (fam(a), secx(a)).
\\ Usage: famq(i,w) / secq(i,w) for i=1..6 (rank 5);  fam6(r) / sec6(r) (rank 6 over Q(r), = pair Q5+Q6).
QS = [81*t^2 + 45*t + 400, 81*t^2 + 144*t + 1280, 29241*t^2 + 13680*t - 12800, 81*t^2 - 45*t + 20, -81*t^2 - 684*t + 156, 810*t^2 - 1980*t + 400];
QX = [6561*t^4 + 102060*t^3 + 414720*t^2 + 403200*t + 1536000, 6561*t^4 + 69984*t^3 + 129600*t^2 + 276480*t + 409600, 6561*t^4 + 25920*t^2, 6561*t^4 + 102060*t^3 + 388800*t^2, -729*t^4 - 4536*t^3 - 720*t^2 - 17920*t + 25600/3, 6561*t^4 + 55404*t^3 - 51840*t^2 - 230400*t];
QT = [(-w^2 + 400)/(18*w - 45), (-w^2 + 1280)/(18*w - 144), (-160*w^2 - 80)/(171*w^2 - 171), (-w^2 + 20)/(18*w + 45), (2*w^2 - 78)/(9*w^2 + 9), (20*w^2 - 20)/(9*w^2 - 90)];
famq(i,w0) = fam(subst(QT[i],'w,w0));
secq(i,w0) = my(t0=subst(QT[i],'w,w0)); concat(secx(t0), [subst(QX[i],'t,t0)]);
T6 = (2*r^4 - 3004*r^2 + 4802)/(9*r^4 - 558*r^2 + 21609);
fam6(r0) = fam(subst(T6,'r,r0));
sec6(r0) = my(t0=subst(T6,'r,r0)); concat(secx(t0), [subst(QX[5],'t,t0), subst(QX[6],'t,t0)]);
\\ Genus-1 pairs (Qi,Qj both square) have many small points -> sporadic rank>=6 fibres: use qpair(Q[i],Q[j],bnd) from qsec.gp.
