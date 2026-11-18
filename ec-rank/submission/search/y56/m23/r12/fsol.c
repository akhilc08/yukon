// Twelve-root Mestre variant: p(x)=f(x)f(x+t), f=prod(x-b_i), 6 roots with e1=0.
// Then r5 = t^2 (s2 s3 - 2 s5) (s_k = coeff of x^(6-k) of f), so deg r<=4 for all t iff e2*e3 == 2*e5.
// This enumerates integer root sets b1<..<b6, sum 0, |b|<=R, with e2 e3 = 2 e5.
#include <stdio.h>
#include <stdlib.h>
int main(int c,char**v){ int R=atoi(v[1]); long n=0;
  for(int b1=-R;b1<=R;b1++) for(int b2=b1+1;b2<=R;b2++) for(int b3=b2+1;b3<=R;b3++) for(int b4=b3+1;b4<=R;b4++) for(int b5=b4+1;b5<=R;b5++){
    int b6=-(b1+b2+b3+b4+b5); if(b6<=b5||b6>R) continue;
    long b[6]={b1,b2,b3,b4,b5,b6}; long e[7]={1,0,0,0,0,0,0};
    for(int i=0;i<6;i++) for(int k=i+1;k>=1;k--) e[k]+=e[k-1]*b[i];
    if(e[2]*e[3]==2*e[5]){ long g=0; for(int i=0;i<6;i++){ long x=labs(b[i]); while(x){ long t=g%x; g=x; x=t; } } if(g!=1) continue; // primitive
      int sym=1; for(int i=0;i<6;i++) if(b[i]!=-b[5-i]) sym=0; if(sym) continue;
      printf("%ld %ld %ld %ld %ld %ld\n",b[0],b[1],b[2],b[3],b[4],b[5]); n++; } }
  fprintf(stderr,"found %ld\n",n); }
