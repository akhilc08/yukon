K=552; read("rank.gp");
t=18371/22532; a=stageA(t);
for(e=1,4, my_t0=getabstime(); b=stageB(t,a[2],e); print(e," ",b[1..2]," #L=",#b[3]," ",getabstime()-my_t0,"ms"));
K=808; read("rank.gp"); t=30/769; a=stageA(t); print(a);
for(e=0,3, my_t0=getabstime(); b=stageB(t,a[2],e); print(e," ",b[1..2]," #L=",#b[3]," ",getabstime()-my_t0,"ms"));
