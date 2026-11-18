read("twsec.gp");
a=famA(552); b=famB(552); D=54*T^2+18*T+1;
r=twfind(a,b,D); print("twfind: ",#r); if(#r, print(r[1]));
\\ also E' twist
r2=twfind(-2*a,a^2-4*b,D); print("twfind E': ",#r2); if(#r2, print(r2[1]));
