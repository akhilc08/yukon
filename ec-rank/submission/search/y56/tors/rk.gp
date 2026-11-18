fA(F,p,q)=if(F==1, q*(121*q-99*p)*(26*p-23*q)*(75*p-44*q), -q*(29*p+2*q)*(400*p+49*q)*(221*p+86*q));
fB(F,p,q)=if(F==1, p*(21*q-74*p)*(25*p-36*q)*(p+77*q), -p*(8*p+53*q)*(25*p+88*q)*(392*p+221*q));
E2(F,p,q)={my(A=fA(F,p,q),B=fB(F,p,q),g=core(gcd(A,B),1)[2]);ellminimalmodel(ellinit([0,(A+B)/g^2,0,A*B/g^4,0]))}
