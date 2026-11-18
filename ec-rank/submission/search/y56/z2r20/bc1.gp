K=552; read("rank.gp");
Ts(s)=(18-2*s)/(s^2-54);
{
foreach([1,2,3,1/2,-1,-3,5/2,7/3,-2/5,11/4],s, my(t=Ts(s),a=stageA(t),b=stageB(t,a[2],1)); print("s=",s," T=",t," H=",round(a[1]*100)/100.," rank ",b[1..2]));
print("--- random T for comparison");
foreach([3/7,-5/11,2/13,17/19,-7/23],t, my(a=stageA(t),b=stageB(t,a[2],1)); print("T=",t," H=",round(a[1]*100)/100.," rank ",b[1..2]));
}
