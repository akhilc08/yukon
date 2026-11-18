\\ Torsion Z/4, generic rank 4 over Q(t) (Mordell-Weil Z/4 x Z^4), K3 fibration of the singular
\\ K3 surface with NS rank 20, disc -163.
\\ Source: N. D. Elkies, "Three lectures on elliptic surfaces and curves of high rank",
\\   arXiv:0709.2908, Lecture II p.9 (equation, 4 generator x-coords; a2 there has a typo);
\\   corrected equation: Elkies-Klagsbrun, arXiv:2003.00077, Section 11 eq (5).
\\ Produced the Z/4 rank 12 (Elkies 2006, Dujella-Peral 2014) and rank 13 (Elkies-Klagsbrun 2020) records.
\\ E1 : y^2 + a xy + a b y = x^3 + b x^2,  a=(8t-1)(32t+7), b=8(t+1)(15t-8)(31t-7); (0,0) has order 4.
fam(t) = { my(a=(8*t-1)*(32*t+7), b=8*(t+1)*(15*t-8)*(31*t-7)); [a, b, a*b, 0, 0]; }
secx(t) = {[-15*(t+1)*(31*t-7)*(32*t+7)/4, (8*t-1)*(15*t-8)*(31*t-7)*(32*t+7),
           -(t+1)*(8*t-1)*(15*t-8)*(32*t+7), -4*(t+1)*(2*t+5)*(15*t-8)*(32*t+7)];}
