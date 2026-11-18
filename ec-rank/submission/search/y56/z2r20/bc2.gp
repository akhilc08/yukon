K=552; read("rank.gp");
Ts(s)=(18-2*s)/(s^2-54);
{
my(L1=List(),L2=List());
foreach([13/5,17/6,-19/7,23/8,-29/9,31/10,37/11,-41/12,43/13,47/14,-9/17,53/19],s, my(t=Ts(s),a=stageA(t),b=stageB(t,a[2],1)); listput(L1,b[1]));
foreach([13/57,-17/61,19/67,-23/71,29/73,-31/79,37/83,41/89,-43/97,47/101,-53/103,59/107],t, my(a=stageA(t),b=stageB(t,a[2],1)); listput(L2,b[1]));
print("base change ranks ",Vec(L1)); print("random T ranks    ",Vec(L2));
}
