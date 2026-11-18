\\ certify.gp: from a .pts record [p,q,H,ainvs,points] write an .ec file: independent saturated points on the minimal model.
\\ cert(rec, file, hdr): hdr = vector of comment strings. Uses ellrank again (seeded) to fill missing generators if needed.
H(E)=log(max(abs(E.c4)^3,E.c6^2))/12;
indep(E,P)={my(B=List(),M);foreach(P,pt, if(ellorder(E,pt)!=0,next); my(T=concat(Vec(B),[pt])); M=ellheightmatrix(E,T); if(matrank(M)==#T, listput(B,pt))); Vec(B)}
cert(rec,file,hdr)={my(E=ellinit(rec[4]),P=rec[5],r,B);
  r=ellrank(E,1,P); B=indep(E,concat(P,r[4]));
  B=ellsaturation(E,B,100); B=indep(E,B);
  system(Str("rm -f ",file));
  foreach(hdr,h, write(file,"# ",h));
  write(file,"# H = ",H(E)," ; certified points ",#B," ; ellrank bounds ",r[1],"..",r[2]," ; ellsaturation(E,L,100)");
  write(file,"E ",strjoin(apply(x->Str(x),E[1..5])," "));
  foreach(B,pt, write(file,"P ",pt[1]," ",pt[2]));
  write(file,"#CLAIM rank=",#B);
  #B}
