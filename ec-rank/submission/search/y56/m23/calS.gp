\\ calibration: S-values of board curves (only for setting sieve thresholds)
S(E,X)=my(s=0.);forprime(p=5,X,if(E.disc%p,s-=ellap(E,p)*log(p)/p));s;
H(E)=log(max(abs(E.c4)^3,E.c6^2))/12;
read("/private/tmp/claude-501/-Users-sickle-Coding-yukon/77853d4d-5d16-4cd0-9e71-3622b7a4cd80/scratchpad/curves.gp");
{foreach(CV,c,my(E=ellinit(c[2])); printf("%s H=%.3f S199=%.2f S3000=%.2f S10000=%.2f\n",c[1],H(E),S(E,199),S(E,3000),S(E,10000)));}
