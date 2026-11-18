\\ Quadratic-section finder for y^2 = x^3 + A(t) x^2 + B(t) x  (A, B in Q[t]).
\\ Enumerates x = c * prod f_i^{e_i} over irreducible factors f_i of B and of A^2-4B,
\\ and reports those for which x(x^2+Ax+B) = k * s(t) * square with deg s <= 2:
\\   deg s = 0 and k square  -> rational section over Q(t)
\\   deg s in {1,2} and k*s(t) = z^2 has a rational point -> base change t = t(w) gives rank +1.
\\ This is the Elkies / Dujella-Peral "quadratic section" step, automated (y56 own tool).
\\ square-root-free helper: P/s is a square up to constant; return a primitive poly whose square times const = P/s
sqfsq(Q) = { my(F=factor(Q), r=1); for(i=1,#F~, r*=F[i,1]^F[i,2]); r; }
sqfp(P) = { my(F=factor(P), r=1); for(i=1,#F~, if(F[i,2]%2, r*=F[i,1])); r; }
\\ does k*s(t) = z^2 have a rational point? (deg s <= 2)
conicpt(s) = {
  my(d=poldegree(s,'t));
  if(d<=0, return(issquare(s)));
  if(d==1, return(1));
  my(a=polcoef(s,2,'t), b=polcoef(s,1,'t), c=polcoef(s,0,'t));
  \\ z^2 = a t^2 + b t + c, homogenise: z^2 = a X^2 + b X Y + c Y^2 -> quadratic form in (X,Y,Z)
  my(M=[a, b/2, 0; b/2, c, 0; 0, 0, -1]); M = M*denominator(M);
  my(sol = qfsolve(M)); type(sol) == "t_COL";
}
qsec(A, B, maxdeg, C) = {
  my(F=concat(factor(B)[,1]~, factor(A^2-4*B)[,1]~), out=List(), E);
  F = Set(apply(f->my(g=f/content(f)); if(pollead(g)<0,-g,g), F)); F = select(f->poldegree(f,'t)>0, F);
  my(n=#F, mult=vector(n,i, max(valuation(B,F[i]), valuation(A^2-4*B,F[i]))));
  forvec(e = vector(n,i,[0,mult[i]]),
    my(base=prod(i=1,n,F[i]^e[i])); if(poldegree(base,'t) > maxdeg, next);
    foreach(C, c, my(x=c*base, f=x*(x^2+A*x+B)); if(f==0, next);
      my(s=sqfp(f), k=simplify(f/sqfsq(f)));
      \\ k = constant cofactor; keep only its squarefree part
      my(kk=core(numerator(k)*denominator(k)));  s = s*kk;
      my(d=poldegree(s,'t));
      if(d<=2 && conicpt(s), listput(out, [d, x, s]))));
  Vec(out);
}

\\ parametrise z^2 = s(t) (deg s = 1 or 2, solvable) : returns t(w) in Q(w)
qparam(s) = {
  my(d=poldegree(s,'t));
  if(d==1, my(b=polcoef(s,1,'t), c=polcoef(s,0,'t)); return(('w^2-c)/b));
  my(a=polcoef(s,2,'t), b=polcoef(s,1,'t), c=polcoef(s,0,'t));
  my(M=[a, b/2, 0; b/2, c, 0; 0, 0, -1]); M = M*denominator(M);
  my(sol=qfsolve(M), P=qfparam(M, sol)); \\ P: 3x3, rows give X,Y,Z as forms in (w^2, w, 1)
  my(v=['w^2, 'w, 1]~, X=(P[1,]*v), Y=(P[2,]*v)); X/Y;
}
\\ all distinct quadratic-section conditions -> list of [s, x(t), t(w)]
qfams(A, B, maxdeg, C) = {
  my(r=qsec(A,B,maxdeg,C), seen=List(), out=List());
  foreach(r, v, if(v[1]==0, next); my(key=v[3]/content(v[3])); if(setsearch(Set(Vec(seen)), key), next);
    listput(seen, key); listput(out, [v[3], v[2], qparam(v[3])]));
  Vec(out);
}

\\ pair two quadratic-section conditions: on the rank-5 base change t=t_i(w), condition j becomes
\\ z^2 = quartic(w). Returns [quartic(w) squarefree part, small points]. If the quartic part has
\\ degree <= 2 -> genus 0 -> rank-6 family over Q(r). Else genus 1: points give rank-6 fibres.
qpair(qi, qj, bnd) = {
  my(tw=qi[3], num=numerator(subst(qj[1],'t,tw)), den=denominator(subst(qj[1],'t,tw)));
  my(P=num*den);  \\ same square class as num/den
  my(s=sqfp(P), k=simplify(P/sqfsq(P))); s = s*core(numerator(k)*denominator(k));
  my(pts = if(poldegree(s,'w)>=3 && poldegree(s,'w)<=4, hyperellratpoints(s, bnd), []));
  [poldegree(s,'w), s, pts];
}

\\ variant: x = c * base * L^2, L = p t + q with |p|<=P, |q|<=Qm (catches sections like -4(t+1)(15t-8)(16t-7)^2)
qsecL(A, B, maxdeg, C, P, Qm) = {
  my(F=concat(factor(B)[,1]~, factor(A^2-4*B)[,1]~), out=List());
  F = Set(apply(f->my(g=f/content(f)); if(pollead(g)<0,-g,g), F)); F = select(f->poldegree(f,'t)>0, F);
  my(n=#F, mult=vector(n,i, max(valuation(B,F[i]), valuation(A^2-4*B,F[i]))));
  forvec(e = vector(n,i,[0,mult[i]]),
    my(base=prod(i=1,n,F[i]^e[i])); if(poldegree(base,'t) > maxdeg-2, next);
    for(p=1,P, for(q=-Qm,Qm, if(gcd(p,q)!=1, next); my(L=p*'t+q);
      foreach(C, c, my(x=c*base*L^2, f=x*(x^2+A*x+B)); if(f==0, next);
        my(s=sqfp(f), k=simplify(f/sqfsq(f))); s = s*core(numerator(k)*denominator(k));
        my(d=poldegree(s,'t)); if(d<=2 && conicpt(s), listput(out, [d, x, s]))))));
  Vec(out);
}
