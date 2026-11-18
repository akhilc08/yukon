// Ten-root sub-family p = f(x) f(x+t), f quintic with roots b1..b5 (sum 0): r4 = c t^2 (5T^2+40 e2 T+512 e2^2-1920 e4), T=t^2.
// r4 = 0 makes y^2 = r(x) a cubic (extra rational point at infinity). Find integer f with D^2 = 38400 e4 - 8640 e2^2,
// T = (-40 e2 +- D)/10 > 0 a rational square (10*(-40 e2 +- D) / 10 ... i.e. 10*T*10 = perfect square).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef long long ll;
static ll isq(ll n){ if(n<0) return -1; ll r=(ll)sqrt((double)n); while(r*r>n) r--; while((r+1)*(r+1)<=n) r++; return r*r==n?r:-1; }
int main(int c,char**v){ int R=atoi(v[1]);
  for(ll b1=-R;b1<=0;b1++) for(ll b2=b1+1;b2<=R;b2++) for(ll b3=b2+1;b3<=R;b3++) for(ll b4=b3+1;b4<=R;b4++){
    ll b5=-(b1+b2+b3+b4); if(b5<=b4||b5>R) continue; ll b[5]={b1,b2,b3,b4,b5}; ll e[6]={1,0,0,0,0,0};
    for(int i=0;i<5;i++) for(int k=i+1;k>=1;k--) e[k]+=e[k-1]*b[i];
    ll D2=38400*e[4]-8640*e[2]*e[2]; ll D=isq(D2); if(D<0) continue;
    for(int s=-1;s<=1;s+=2){ ll n=-40*e[2]+s*D; if(n<=0) continue; // T = n/10 ; t = sqrt(10 n)/10
      ll w=isq(10*n); if(w<0) continue; printf("%lld %lld %lld %lld %lld  w=%lld (t=w/10)\n",b1,b2,b3,b4,b5,w); if(D==0) break; } }
  return 0; }
