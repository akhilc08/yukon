// psieve: height-first Mestre-Nagao sieve over two-integer-parameter (p,q) families whose a_ell depends only on
// t = p/q in P^1(F_ell) (weighted-homogeneous models: scaling (p,q) by lambda is an isomorphism up to an even twist).
// Per-ell tables a_ell[t] for t in F_ell u {inf}; inner loop is table lookups with a running residue per prime.
// Stages: approximate naive height (cheap, long double) -> S(P1) -> S(P2) -> S(P3), with early abort.
// S = sum_{3<=ell<=P} -a_ell log(ell)/ell over primes of good reduction (bad residues weight 0).
// Families: 1 X400a, 2 X312 (Z/2xZ/2, y^2=x(x+A)(x+B), A,B quartic forms from board headers; common square removed),
//           3 z3k3 (Z/3, y^2+a1xy+a3y=x^3, a1=7p^2+309q^2, a3=-91(p^2-q^2)(p^2-25q^2)(p^2-121q^2), from board 3-r12-z3k3.ec).
// Build: gcc -O3 -mcpu=native -fopenmp -o psieve psieve.c -lm
// Usage: psieve FAM QMAX QFULL HMAX T1 T2 T3 [Q0 W PFULLMUL P1 P2 P3]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
static double omp_get_wtime(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
#endif
typedef __int128 i128; typedef long double LD;
static int fam;
static int P1=200,P2=1000,P3=5000;
static int np, pr[1000]; static int n1,n2; static int8_t *TAB[1000]; static float wt[1000];
static double R[16]; static int nR;
static long long md(long long x,long long m){x%=m;return x<0?x+m:x;}
// integral a-invariants at (p,q) mod m (m>0) or as long double (m=0 path separate)
static void ainv_mod(long long p,long long q,long long m,long long a[5]){
  #define M(x) md(x,m)
  memset(a,0,5*sizeof(long long));
  if(fam==1||fam==2){ long long A,B;
    if(fam==1){A=M(M(q*M(121*q-99*p))*M(M(26*p-23*q)*M(75*p-44*q))); B=M(M(p*M(21*q-74*p))*M(M(25*p-36*q)*M(p+77*q)));}
    else {A=M(-M(q*M(29*p+2*q))*M(M(400*p+49*q)*M(221*p+86*q))); B=M(-M(p*M(8*p+53*q))*M(M(25*p+88*q)*M(392*p+221*q)));}
    a[1]=M(A+B); a[3]=M(A*B);
  } else if(fam==3){ a[0]=M(7*M(p*p)+309*M(q*q)); long long pp=M(p*p),qq=M(q*q);
    a[2]=M(-91*M(M(pp-qq)*M(pp-25*qq))); a[2]=M(a[2]*M(pp-121*qq)); }
}
static void ainv_ld(long long p,long long q,LD a[5]){
  LD P=p,Q=q; for(int i=0;i<5;i++)a[i]=0;
  if(fam==1||fam==2){ LD A,B;
    if(fam==1){A=Q*(121*Q-99*P)*(26*P-23*Q)*(75*P-44*Q); B=P*(21*Q-74*P)*(25*P-36*Q)*(P+77*Q);}
    else {A=-Q*(29*P+2*Q)*(400*P+49*Q)*(221*P+86*Q); B=-P*(8*P+53*Q)*(25*P+88*Q)*(392*P+221*Q);}
    a[1]=A+B; a[3]=A*B;
  } else if(fam==3){ a[0]=7*P*P+309*Q*Q; a[2]=-91*(P*P-Q*Q)*(P*P-25*Q*Q)*(P*P-121*Q*Q); }
}
static LD Hld(const LD a[5]){ LD b2=a[0]*a[0]+4*a[1], b4=a[0]*a[2]+2*a[3], b6=a[2]*a[2]+4*a[4];
  LD c4=b2*b2-24*b4, c6=-b2*b2*b2+36*b2*b4-216*b6; LD x=3*logl(fabsl(c4)+1), y=2*logl(fabsl(c6)+1); return (x>y?x:y)/12; }
static int apcount(const long long a[5], long long l, const int8_t*chi){
  long long b2=md(a[0]*a[0]+4*a[1],l), b4=md(a[0]*a[2]+2*a[3],l), b6=md(a[2]*a[2]+4*a[4],l);
  long long b8=md(a[0]*a[0]%l*a[4]+4*a[1]*a[4]%l-a[0]*a[2]%l*a[3]+a[1]*a[2]%l*a[2]-a[3]*a[3],l);
  long long d=md(-b2*b2%l*b8-8*b4%l*b4%l*b4-27*b6*b6%l+9*b2*b4%l*b6,l); if(!d) return 99;
  long long s=0; for(long long x=0;x<l;x++){ long long u=md(a[0]*x+a[2],l); long long f=md(((x+a[1])*x%l+a[3])%l*x+a[4],l); s+=chi[md(u*u+4*f,l)]; }
  return (int)(-s);
}
int main(int argc,char**argv){
  fam=atoi(argv[1]); long long QM=atoll(argv[2]),QF=atoll(argv[3]); double HM=atof(argv[4]),T1=atof(argv[5]),T2=atof(argv[6]),T3=atof(argv[7]);
  long long Q0=argc>8?atoll(argv[8]):1, W=argc>9?atoll(argv[9]):60; double PMUL=argc>10?atof(argv[10]):80;
  if(argc>11)P1=atoi(argv[11]); if(argc>12)P2=atoi(argv[12]); if(argc>13)P3=atoi(argv[13]);
  if(fam==1){double r[]={121./99,23./26,44./75,36./25,21./74,-77,0}; nR=7; memcpy(R,r,sizeof r);}
  if(fam==2){double r[]={-2./29,-49./400,-86./221,-53./8,-88./25,-221./392,0}; nR=7; memcpy(R,r,sizeof r);}
  if(fam==3){double r[]={1,5,11}; nR=3; memcpy(R,r,sizeof r);}
  char*sv=calloc(P3+1,1); np=0; for(int i=2;i<=P3;i++) if(!sv[i]){ if(i>3) pr[np++]=i; for(int j=2*i;j<=P3;j+=i) sv[j]=1; }
  // include 3 for non-Z/3 families
  if(fam!=3){ memmove(pr+1,pr,np*sizeof(int)); pr[0]=3; np++; }
  n1=n2=0; for(int k=0;k<np;k++){ if(pr[k]<=P1) n1=k+1; if(pr[k]<=P2) n2=k+1; }
  double t0=omp_get_wtime();
  #pragma omp parallel for schedule(dynamic,1)
  for(int k=0;k<np;k++){ long long l=pr[k]; int8_t*chi=malloc(l); for(long long x=0;x<l;x++) chi[x]=-1; chi[0]=0; for(long long x=1;x<l;x++) chi[x*x%l]=1;
    TAB[k]=malloc(l+1); long long a[5];
    for(long long t=0;t<=l;t++){ if(t<l) ainv_mod(t,1,l,a); else ainv_mod(1,0,l,a); int v=apcount(a,l,chi); TAB[k][t]=(v==99)?0:(int8_t)v; }
    wt[k]=(float)(log((double)l)/l); free(chi); }
  fprintf(stderr,"tables %d primes (n1=%d n2=%d) in %.1fs\n",np,n1,n2,omp_get_wtime()-t0);
  long long tot=0,c0=0,c1=0,c2=0,c3=0; t0=omp_get_wtime();
  #pragma omp parallel for schedule(dynamic,4) reduction(+:tot,c0,c1,c2,c3)
  for(long long q=Q0;q<=QM;q++){
    int *iq=malloc(np*sizeof(int)), *res=malloc(np*sizeof(int)); char line[256];
    for(int k=0;k<np;k++){ long long l=pr[k],x=q%l; if(!x){iq[k]=-1;continue;} long long r=1,b=x,e=l-2; while(e){if(e&1)r=r*b%l;b=b*b%l;e>>=1;} iq[k]=(int)r; }
    int nseg = q<=QF?1:nR;
    for(int sg=0;sg<nseg;sg++){
      long long lo,hi; if(q<=QF){lo=fam==3?1:-(long long)(PMUL*q)-500;hi=(long long)(PMUL*q)+500;} else {long long c=llroundl((LD)R[sg]*q); lo=c-W; hi=c+W;}
      for(int k=0;k<np;k++){ long long l=pr[k]; res[k]= iq[k]<0 ? (int)l : (int)(md(lo,l)*iq[k]%l); }
      for(long long p=lo;p<=hi;p++){
        if(p!=lo) for(int k=0;k<n1;k++){ if(iq[k]>=0){ int l=pr[k]; int r=res[k]+iq[k]; if(r>=l) r-=l; res[k]=r; } }
        tot++;
        long long g0=p<0?-p:p, h0=q; while(h0){long long tt=g0%h0;g0=h0;h0=tt;} if(g0!=1) continue;
        LD a[5]; ainv_ld(p,q,a);
        LD h;
        if(fam<=2){ // remove common square of A,B exactly
          i128 P=p,Q=q,A,B; if(fam==1){A=Q*(121*Q-99*P)*(26*P-23*Q)*(75*P-44*Q); B=P*(21*Q-74*P)*(25*P-36*Q)*(P+77*Q);}
          else {A=-Q*(29*P+2*Q)*(400*P+49*Q)*(221*P+86*Q); B=-P*(8*P+53*Q)*(25*P+88*Q)*(392*P+221*Q);}
          if(A==0||B==0||A==B) continue;
          LD m=fabsl((LD)A)>fabsl((LD)B)?fabsl((LD)A):fabsl((LD)B); if(logl(m) > 2*(HM+0.75)+2*logl(1e7)) continue;
          i128 x=A<0?-A:A,y=B<0?-B:B; while(y){i128 tt=x%y;x=y;y=tt;} long long s=1;
          for(long long d=2; d<100000 && (i128)d*d<=x; d++){ i128 dd=(i128)d*d; while(x%dd==0){x/=dd;s*=d;} if(x%d==0)x/=d; }
          LD s2=(LD)s*s; a[1]/=s2; a[3]/=s2*s2;
        }
        h=Hld(a); if(h>HM+0.75) continue; c0++;
        float S=0; int k=0;
        for(;k<n1;k++) S-=TAB[k][res[k]]*wt[k];
        if(S<T1) continue; c1++;
        for(int j=n1;j<np;j++){ long long l=pr[j]; res[j]= iq[j]<0 ? (int)l : (int)(md(p,l)*iq[j]%l); }
        for(;k<n2;k++) S-=TAB[k][res[k]]*wt[k];
        if(S<T2) continue; c2++;
        float S2=S;
        for(;k<np;k++) S-=TAB[k][res[k]]*wt[k];
        if(S<T3) continue; c3++;
        snprintf(line,sizeof line,"%lld %lld %.4f %.3f %.3f\n",p,q,(double)h,S2,S);
        #pragma omp critical
        fputs(line,stdout);
      }
    }
    free(iq); free(res);
  }
  double dt=omp_get_wtime()-t0;
  fprintf(stderr,"scanned %lld height-ok %lld S1 %lld S2 %lld S3 %lld in %.1fs (%.3g cand/s total)\n",tot,c0,c1,c2,c3,dt,tot/dt);
}
