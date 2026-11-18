// Faster fsol: b1<b2<b3<b4 chosen, (b5,b6) from sigma=-(b1+..+b4) and pi solving the quadratic
// E1 pi^2 + (A E1 + Bc - 2 E3) pi + (A Bc - 2 sigma E4) = 0, A=E2+sigma E1, Bc=E3+sigma E2 (E_k = e_k(b1..b4)),
// which is e2 e3 = 2 e5 for the six roots with e1 = 0. Prints primitive, non-symmetric sets with b4<b5<b6<=R, b1+b6>=0.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef long long ll; typedef __int128 i128;
static ll isq(i128 n){ if(n<0) return -1; ll r=(ll)sqrtl((long double)n); while((i128)r*r>n) r--; while((i128)(r+1)*(r+1)<=n) r++; return (i128)r*r==n?r:-1; }
static ll g2(ll a,ll b){ a=llabs(a); b=llabs(b); while(b){ ll t=a%b; a=b; b=t; } return a; }
static void emit(ll b1,ll b2,ll b3,ll b4,ll p,ll s,int R){
  i128 D=(i128)s*s-4*(i128)p; ll d=isq(D); if(d<0) return; if(((s+d)&1)) return; ll b5=(s-d)/2, b6=(s+d)/2;
  if(!(b4<b5 && b5<b6 && b6<=R)) return; if(b1+b6<0) return;
  ll g=g2(g2(g2(b1,b2),g2(b3,b4)),g2(b5,b6)); if(g!=1) return;
  if(b1==-b6 && b2==-b5 && b3==-b4) return;
  printf("%lld %lld %lld %lld %lld %lld\n",b1,b2,b3,b4,b5,b6); }
int main(int c,char**v){ int R=atoi(v[1]);
  for(ll b1=-R;b1<=0;b1++) for(ll b2=b1+1;b2<=R;b2++) for(ll b3=b2+1;b3<=R;b3++) for(ll b4=b3+1;b4<=R;b4++){
    ll E1=b1+b2+b3+b4, E2=b1*b2+b1*b3+b1*b4+b2*b3+b2*b4+b3*b4, E3=b1*b2*b3+b1*b2*b4+b1*b3*b4+b2*b3*b4, E4=b1*b2*b3*b4;
    ll s=-E1; if(s<=2*b4) continue; // need b5,b6 > b4
    i128 A=E2+(i128)s*E1, Bc=E3+(i128)s*E2; i128 qa=E1, qb=A*E1+Bc-2*(i128)E3, qc=A*Bc-2*(i128)s*E4;
    if(qa==0){ if(qb==0) continue; if(qc%qb) continue; emit(b1,b2,b3,b4,(ll)(-qc/qb),s,R); continue; }
    i128 D=qb*qb-4*qa*qc; if(D<0) continue; long double sd=sqrtl((long double)D); ll r=(ll)sd; i128 rr=r; while(rr*rr>D) rr--; while((rr+1)*(rr+1)<=D) rr++; if(rr*rr!=D) continue;
    for(int sg=-1;sg<=1;sg+=2){ i128 num=-qb+sg*rr, den=2*qa; if(num%den) continue; emit(b1,b2,b3,b4,(ll)(num/den),s,R); if(rr==0) break; } }
  return 0; }
