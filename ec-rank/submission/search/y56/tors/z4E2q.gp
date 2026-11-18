\\ z4E2q: own rank-5 base changes of the EK disc -163 Z/4 family E2 (fam/z4_ek163_E2.gp) by quadratic sections,
\\ found by y56 qsec.gp on E2 shifted to y^2 = X(X^2+A X+B) (2-torsion root x0 below).
\\ Q1: section X = Q1X(t) is rational iff 20800 t^2 - 28160 t + 1681 = z^2  ->  t = (82 m + 28160)/(20800 - m^2)
\\ Q2: section X = Q2X(t) is rational iff 102724 t^2 + 9900 t + 2025 = z^2 ->  t = (9900 - 90 m)/(m^2 - 102724)
read("fam/z4_ek163_E2.gp");
E2x0(t) = -(144*t^4-6200*t^3+15068*t^2-6650*t+324);
E2Q1X(t)=144*t^4 - 440*t^3 + 348*t^2 - 90*t + 4;
E2Q2X(t)=144*t^4 - 5984*t^3 + 6200*t^2 - 1784*t + 81;
tQ1(m) = (82*m+28160)/(20800-m^2);
tQ2(m) = (9900-90*m)/(m^2-102724);
famQ1(m) = fam(tQ1(m));
secQ1(m) = {my(t0=tQ1(m)); concat(secx(t0), [E2Q1X(t0)+E2x0(t0)]);}
famQ2(m) = fam(tQ2(m));
secQ2(m) = {my(t0=tQ2(m)); concat(secx(t0), [E2Q2X(t0)+E2x0(t0)]);}
