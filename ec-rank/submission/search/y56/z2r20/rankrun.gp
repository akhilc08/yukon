\\ env: FAMK (808|552), IN (survivor file: p q ...), OUT (results), HB (height bound), SB (min S(1e5) for ellrank), HITS (dir)
K=eval(getenv("FAMK")); read("rank.gp");
IN=getenv("IN"); OUT=getenv("OUT"); HB=eval(getenv("HB")); SB=eval(getenv("SB")); HITS=getenv("HITS");
\\ cells we care about: bound per rank (records.json baselines and current board heights)
bnd(r)=if(r>=20,27.547796, r==19,26.528271, r==18,21.507228, r==17,18.911354, r==16,17.507728, 0);
{
my(L=readstr(IN));
for(i=1,#L, my(w=strsplit(L[i]," "),t,a,b);
  if(#w<2,next); t=eval(w[1])/eval(w[2]);
  iferr(a=stageA(t), e, write(OUT,w[1]," ",w[2]," ERR"); next);
  if(a[1]>=HB, write(OUT,w[1]," ",w[2]," H ",a[1]); next);
  if(a[3]<SB, write(OUT,w[1]," ",w[2]," H ",a[1]," S5 ",a[3]); next);
  iferr(b=stageB(t,a[2],1), e, write(OUT,w[1]," ",w[2]," RKERR"); next);
  write(OUT,w[1]," ",w[2]," H ",a[1]," S5 ",a[3]," R ",b[1]," ",b[2]," j ",a[2]);
  if(b[1]>=16 && a[1] < bnd(b[1])+0.0, 
     write(HITS"/K"K"_"w[1]"_"w[2]".txt", [K,w[1],w[2],a[1],b[1],b[2],a[2]], "\n", b[4][1..5], "\n", b[3])));
}
