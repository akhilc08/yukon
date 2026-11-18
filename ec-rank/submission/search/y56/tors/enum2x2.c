// enum2x2: enumerate T=p/q in rank-5 Z/2xZ/2 families y^2=x(x+A)(x+B) (X400a, X312; equations from board headers),
// reduce common square factor, approximate naive height of the Legendre model, print survivors with H<HMAX+SLACK.
// Also flags whether one of -A,-B,A,B,A-B,B-A is a square (then a 2-isogenous curve has a Z/4 point).
// Usage: enum2x2 FAM QMAX QFULL HMAX [Q0 W FLAGONLY]  (FAM 1=X400a 2=X312); q<=QFULL: |p|<=80q+500; QFULL<q<=QMAX: strips |p-rq|<=W around roots r
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef __int128 i128;
static int fam;
static void AB(long long p,long long q,i128*A,i128*B){
  i128 P=p,Q=q;
  if(fam==1){*A=Q*(121*Q-99*P)*(26*P-23*Q)*(75*P-44*Q); *B=P*(21*Q-74*P)*(25*P-36*Q)*(P+77*Q);}
  else {*A=-Q*(29*P+2*Q)*(400*P+49*Q)*(221*P+86*Q); *B=-P*(8*P+53*Q)*(25*P+88*Q)*(392*P+221*Q);}
}
static i128 ab(i128 x){return x<0?-x:x;}
static i128 g128(i128 a,i128 b){a=ab(a);b=ab(b);while(b){i128 t=a%b;a=b;b=t;}return a;}
static long long gcdl(long long a,long long b){if(a<0)a=-a;if(b<0)b=-b;while(b){long long t=a%b;a=b;b=t;}return a;}
static int issq(i128 x){ if(x<0) return 0; if(x==0) return 1; long double r=sqrtl((long double)x); i128 s=(i128)r; for(i128 t=s-2;t<=s+2;t++) if(t>=0&&t*t==x) return 1; return 0;}
static void pr(i128 x){ char b[64]; int n=0; int neg=x<0; if(neg)x=-x; if(x==0)b[n++]='0'; while(x){b[n++]='0'+(int)(x%10);x/=10;} if(neg)putchar('-'); while(n)putchar(b[--n]); }
int main(int argc,char**argv){
  fam=atoi(argv[1]); long long QM=atoll(argv[2]), PM=atoll(argv[3]) /* QFULL */; double HM=atof(argv[4]); long long Q0=argc>5?atoll(argv[5]):1;
  double lb = 2*(HM+0.75) - log(16.0)/2.0; // log max(|A|,|B|) bound (generous)
  long long cnt=0;
  long long W=argc>6?atoll(argv[6]):60; int FO=argc>7?atoi(argv[7]):0;
  double R1[]={121./99,23./26,44./75,36./25,21./74,-77,0}, R2[]={-2./29,-49./400,-86./221,-53./8,-88./25,-221./392,0};
  double *R=fam==1?R1:R2;
  for(long long q=Q0;q<=QM;q++){
   int nseg = q<=PM ? 1 : 7; 
   for(int sg=0;sg<nseg;sg++){
    long long lo,hi; if(q<=PM){lo=-80*q-500;hi=80*q+500;} else {long long c=llround(R[sg]*q); lo=c-W; hi=c+W;}
    for(long long p=lo;p<=hi;p++){
    if(p==0) continue;
    // fast double prefilter
    double P=p,Q=q,Ad,Bd;
    if(fam==1){Ad=Q*(121*Q-99*P)*(26*P-23*Q)*(75*P-44*Q); Bd=P*(21*Q-74*P)*(25*P-36*Q)*(P+77*Q);}
    else {Ad=-Q*(29*P+2*Q)*(400*P+49*Q)*(221*P+86*Q); Bd=-P*(8*P+53*Q)*(25*P+88*Q)*(392*P+221*Q);}
    double m=fmax(fabs(Ad),fabs(Bd)); if(m==0) continue;
    if(log(m) > lb + 2*log(1e6)) continue; // allow common square up to 1e6
    if(gcdl(p,q)!=1) continue;
    i128 A,B; AB(p,q,&A,&B); if(A==0||B==0||A==B) continue;
    i128 g=g128(A,B), s=1;
    for(long long d=2; d*d<= (long long)( g<((i128)1<<62)? (long long)g : (1LL<<62)) && d<200000; d++){ i128 dd=(i128)d*d; while(g%dd==0){g/=dd;s*=d;} if(g%d==0) g/=d; }
    i128 s2=s*s; A/=s2; B/=s2;
    long double a2=(long double)(A+B), a4=(long double)A*(long double)B;
    long double c4=16*(a2*a2-3*a4), c6=-32*a2*(2*a2*a2-9*a4);
    long double h=fmaxl(3*logl(fabsl(c4)+1),2*logl(fabsl(c6)+1))/12;
    if(h> HM+0.75) continue;
    int fl=0; i128 v[6]={-A,-B,A,B,A-B,B-A}; for(int k=0;k<6;k++) if(issq(v[k])) fl|=1<<k;
    if(FO && !fl) continue;
    printf("%lld %lld ",p,q); pr(A); putchar(' '); pr(B); printf(" %.4f %d\n",(double)h,fl); cnt++;
  }}}
  fprintf(stderr,"survivors %lld\n",cnt);
}
