\\ Torsion Z/4, generic rank 4 claimed (second Z^4 x Z/4 fibration of the disc -163 singular K3).
\\ Source: Elkies-Klagsbrun, arXiv:2003.00077, Section 11 eq (6) (text garbled in the PDF; the
\\   a1 x y term is restored from the Tate shape a3 = a1*a2, which the formula satisfies);
\\   also quoted in Dujella-Peral arXiv:2207.08206 Sec. 3 as a=4(9+80t), b=(t-2)(2t-81)(2t-1)(18t-1)/2.
\\ EK searched it to height 2^22: 32 rank-12 fibres, no rank 13.
\\ Sections: NOT published. The x-coords below were found here by brute-force search for polynomial
\\ sections (c * product of bad-fibre factors * small quadratic, deg <= 4); their rank was then
\\ checked by ellheightmatrix at specialisations: only 3 are independent (tested rank 3; 4th generator not found).
fam(t) = { my(d=(t-2)*(2*t-1)*(18*t-1)*(2*t-81), a1=-8*(80*t+9)); [a1, 2*d, 2*a1*d, 0, 0]; }
secx(t) = { [5120*t^3-12224*t^2+3680*t+576, -144*t^4+6200*t^3-15068*t^2+6650*t-324,
  432*t^4-19464*t^3+81972*t^2-91974*t+4860, -108*t^4+4704*t^3-13518*t^2+6204*t-1215/4,
  138240*t^4-291648*t^3+28160*t^2+4496*t-288, 17280*t^4+29784*t^3-17188*t^2-1246*t+117]; }
