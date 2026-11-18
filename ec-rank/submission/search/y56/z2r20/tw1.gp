read("twsec.gp"); a=famA(552); b=famB(552); D=(T-22)*(7*T+17);
t0=getabstime(); r=twfind(a,b,D); print(#r," ",getabstime()-t0,"ms"); print(r);
