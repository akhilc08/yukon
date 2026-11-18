\\ Twelve-root Mestre variant: roots b_i and b_i - t, with sum b = 0 and e2 e3 = 2 e5 (see fsol.c).
\\ p = f(x) f(x+t) = q^2 - r with deg r <= 4; E = Jac(y^2 = r).
read(Str(ROOT,"/submission/search/y56/m23/lib.gp"));
r12E(b,t)=mesE(concat(b,apply(z->z-t,b)));
