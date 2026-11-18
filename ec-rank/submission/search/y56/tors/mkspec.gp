\\ mkspec.gp: turn a library family fam(t) (a-invariants, rational functions in one parameter) into an integral
\\ weighted-homogeneous spec for gsieve.c:  a_i(p,q) = q^(w*i) * lam(p/q)^i * a_i(p/q) with lam clearing denominators.
\\ writespec(F, file): F = closure u -> [a1,a2,a3,a4,a6]. Spec file lines: "w", then for i in 1,2,3,4,6: "deg c0 c1 ... cdeg".
\\ Also prints the real roots of the discriminant (strip centres) on the last line.
mkpolys(F)={my(a=F('u),lam=1,d,b,w=0,ws=[1,2,3,4,6],den);
  for(k=1,5, if(a[k]!=0, lam=lcm(lam, denominator(a[k]))));
  b=vector(5,k, a[k]*lam^ws[k]);
  den=content(b[1..5]); \\ integer content; leave it (minimal model fixes)
  b=vector(5,k, my(v=b[k]); v=v*denominator(content(v)); v); \\ ensure integral (scales hurt nothing? no: must keep model) 
  [b,lam]}
mkpolys2(F)={my(a=F('u),lam=1,b,ws=[1,2,3,4,6],D=1,fs=List());
  for(k=1,5, if(a[k]!=0, my(fa=factor(denominator(a[k]))); for(j=1,#fa~, if(poldegree(fa[j,1])>0, listput(fs,fa[j,1])))));
  fs=Set(Vec(fs));
  foreach(fs,f, my(e=0); for(k=1,5, if(a[k]!=0, my(m=valuation(denominator(a[k]),f)); e=max(e,ceil(m/ws[k])))); lam*=f^e);
  b=vector(5,k, a[k]*lam^ws[k]);
  for(k=1,5, if(b[k]!=0, my(c=content(b[k]),fa=factor(denominator(c))); for(j=1,#fa~, my(p=fa[j,1],m=fa[j,2]); D=lcm(D,p^ceil(m/ws[k]))) ));
  b=vector(5,k, b[k]*D^ws[k]);
  b}
writespec(F,file)={my(b=mkpolys2(F),ws=[1,2,3,4,6],w=0,dsc,rr,s);
  for(k=1,5, if(b[k]!=0, w=max(w, ceil(poldegree(b[k])/ws[k]))));
  write1(file,""); system(Str("rm -f ",file));
  write(file,w);
  for(k=1,5, my(v=Pol(b[k],'u),dg=w*ws[k],s=Str(dg)); for(j=0,dg, s=Str(s," ",polcoef(v,j,'u))); write(file,s));
  dsc=discw(b); rr=List(); foreach(factor(dsc)[,1],f, if(poldegree(f)>0, foreach(polroots(f),z, if(abs(imag(z))<1e-20 && abs(z)<1e6, listput(rr,real(z))))));
  rr=Vec(rr);
  write(file, #rr, " ", strjoin(apply(x->Str(round(x*1e9)/1e9*1.),rr)," "));
  b}
discw(a)={my(b2=a[1]^2+4*a[2],b4=a[1]*a[3]+2*a[4],b6=a[3]^2+4*a[5],b8=a[1]^2*a[5]+4*a[2]*a[5]-a[1]*a[3]*a[4]+a[2]*a[3]^2-a[4]^2);-b2^2*b8-8*b4^3-27*b6^2+9*b2*b4*b6}
