\\ m23 helpers (PARI/GP 2.17.x). Ten-root Mestre: p=prod(x-a_i)=q^2-r; E = Jac(y^2=r(x)).
mesqr(a)={my(p=prod(i=1,#a,'x-a[i]),n=#a/2); my(s=Ser(subst(p,'x,1/'x)*'x^(2*n),'x)); my(sq=sqrt(s)); my(q=sum(i=0,n,polcoef(sq,i,'x)*'x^(n-i))); [q, q^2-p]};
quartIJ(r)={my(a=polcoef(r,4),b=polcoef(r,3),c=polcoef(r,2),d=polcoef(r,1),e=polcoef(r,0)); [12*a*e-3*b*d+c^2, 72*a*c*e+9*b*c*d-27*a*d^2-27*e*b^2-2*c^3]};
mesE(a)={my(r=mesqr(a)[2],IJ=quartIJ(r)); ellinit(ellminimalmodel(ellinit([0,0,0,-27*IJ[1],-27*IJ[2]])))};
Hn(E)=log(max(abs(E.c4)^3,E.c6^2))/12;
Sx(E,X)=my(s=0.);forprime(p=5,X,if(E.disc%p,s-=ellap(E,p)*log(p)/p));s;
