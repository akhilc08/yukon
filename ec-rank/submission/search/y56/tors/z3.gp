Ez3(p,q)={my(a1=7*p^2+309*q^2,a3=-91*(p^2-q^2)*(p^2-25*q^2)*(p^2-121*q^2));ellminimalmodel(ellinit([a1,0,a3,0,0]))}
H(E)=log(max(abs(E.c4)^3,E.c6^2))/12;
