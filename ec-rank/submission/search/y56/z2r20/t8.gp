K=552; read("rank.gp");
foreach([-14/15,22/13,136/97],t, a=stageA(t); b=stageB(t,a[2],1); print(K," ",t," H=",a[1]," S5=",a[3]," rank ",b[1..2]));
K=808; read("rank.gp");
foreach([-71/1042,-82/965],t, a=stageA(t); b=stageB(t,a[2],1); print(K," ",t," H=",a[1]," S5=",a[3]," rank ",b[1..2]));
