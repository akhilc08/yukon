// gsieve: generic height-first Mestre-Nagao sieve over a one-parameter family given as an integral weighted-homogeneous
// spec (from mkspec.gp): a_i(p,q) = sum_j c_ij p^j q^(w i - j), so the curve at coprime (p,q) depends (up to an even-weight
// isomorphism) only on t=p/q in P^1(F_ell) -> per-ell tables a_ell[t]. Stages: approx naive height of the spec model
// (an upper bound for the minimal model; use slack) -> S(P1) -> S(P2) -> S(P3) with early abort.
// Build: gcc -O3 -mcpu=neoverse-v1 -fopenmp -o gsieve gsieve.c -lm
// Usage: gsieve SPEC QMAX QFULL HMAX T1 T2 T3 [PMUL P0 W SYM P1 P2 P3]
//   q<=QFULL: |p|<=PMUL*q+P0 ; QFULL<q<=QMAX: strips |p-r q|<=W around real discriminant roots r. SYM=1: p>0 only.
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
typedef long double LD;
static int W, DEG[5]; static char *CS[5][200]; static LD CL[5][200];
static int np, pr[2000], n1, n2; static int8_t *TAB[2000]; static float wt[2000];
static double R[64]; static int nR;
static long long md(long long x,long long m){x%=m;return x<0?x+m:x;}
static long long strmod(const char*s,long long m){int neg=0; long long r=0; if(*s=='-'){neg=1;s++;} for(;*s;s++) if(*s>='0'&&*s<='9') r=(r*10+(*s-'0'))%m; return neg?md(-r,m):r;}
static void ainv_mod(long long p,long long q,long long l,long long a[5],long long *C[5]){
  for(int i=0;i<5;i++){ int d=DEG[i]; long long s=0;
    // homogeneous Horner: sum c_j p^j q^(d-j)
    long long pp=1; long long *qp=malloc((d+1)*sizeof(long long)); qp[0]=1; for(int j=1;j<=d;j++) qp[j]=qp[j-1]*md(q,l)%l;
    for(int j=0;j<=d;j++){ s=(s+C[i][j]*pp%l*qp[d-j])%l; pp=pp*md(p,l)%l; }
    free(qp); a[i]=s; }
}
static void ainv_ld(long long p,long long q,LD a[5]){
  LD P=p,Q=q;
  for(int i=0;i<5;i++){ int d=DEG[i]; LD s=0, pp=1; LD qp[200]; qp[0]=1; for(int j=1;j<=d;j++) qp[j]=qp[j-1]*Q;
    for(int j=0;j<=d;j++){ s+=CL[i][j]*pp*qp[d-j]; pp*=P; } a[i]=s; }
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
  FILE*f=fopen(argv[1],"r"); if(!f){perror("spec");return 1;}
  static char buf[1<<20]; fgets(buf,sizeof buf,f); W=atoi(buf);
  for(int i=0;i<5;i++){ fgets(buf,sizeof buf,f); char*tok=strtok(buf," \n"); DEG[i]=atoi(tok);
    for(int j=0;j<=DEG[i];j++){ tok=strtok(NULL," \n"); CS[i][j]=strdup(tok); CL[i][j]=strtold(tok,NULL);} }
  fgets(buf,sizeof buf,f); { char*tok=strtok(buf," \n"); nR=atoi(tok); for(int k=0;k<nR;k++){tok=strtok(NULL," \n"); R[k]=atof(tok);} } fclose(f);
  long long QM=atoll(argv[2]),QF=atoll(argv[3]); double HM=atof(argv[4]),T1=atof(argv[5]),T2=atof(argv[6]),T3=atof(argv[7]);
  double PMUL=argc>8?atof(argv[8]):20; long long P0=argc>9?atoll(argv[9]):500, WS=argc>10?atoll(argv[10]):60; int SYM=argc>11?atoi(argv[11]):0;
  int P1=argc>12?atoi(argv[12]):200, P2=argc>13?atoi(argv[13]):1000, P3=argc>14?atoi(argv[14]):5000;
  char*sv=calloc(P3+1,1); np=0; for(int i=2;i<=P3;i++) if(!sv[i]){ if(i>2) pr[np++]=i; for(int j=2*i;j<=P3;j+=i) sv[j]=1; }
  n1=n2=0; for(int k=0;k<np;k++){ if(pr[k]<=P1) n1=k+1; if(pr[k]<=P2) n2=k+1; }
  double t0=omp_get_wtime();
  #pragma omp parallel for schedule(dynamic,1)
  for(int k=0;k<np;k++){ long long l=pr[k]; int8_t*chi=malloc(l); for(long long x=0;x<l;x++) chi[x]=-1; chi[0]=0; for(long long x=1;x<l;x++) chi[x*x%l]=1;
    long long *C[5]; for(int i=0;i<5;i++){ C[i]=malloc((DEG[i]+1)*sizeof(long long)); for(int j=0;j<=DEG[i];j++) C[i][j]=strmod(CS[i][j],l); }
    TAB[k]=malloc(l+1); long long a[5];
    for(long long t=0;t<=l;t++){ if(t<l) ainv_mod(t,1,l,a,C); else ainv_mod(1,0,l,a,C); int v=apcount(a,l,chi); TAB[k][t]=(v==99)?0:(int8_t)v; }
    wt[k]=(float)(log((double)l)/l); free(chi); for(int i=0;i<5;i++) free(C[i]); }
  fprintf(stderr,"tables %d primes (n1=%d n2=%d) in %.1fs; w=%d nroots=%d\n",np,n1,n2,omp_get_wtime()-t0,W,nR);
  long long tot=0,c0=0,c1=0,c2=0,c3=0; t0=omp_get_wtime();
  #pragma omp parallel for schedule(dynamic,4) reduction(+:tot,c0,c1,c2,c3)
  for(long long q=1;q<=QM;q++){
    int *iq=malloc(np*sizeof(int)), *res=malloc(np*sizeof(int)); char line[256];
    for(int k=0;k<np;k++){ long long l=pr[k],x=q%l; if(!x){iq[k]=-1;continue;} long long r=1,b=x,e=l-2; while(e){if(e&1)r=r*b%l;b=b*b%l;e>>=1;} iq[k]=(int)r; }
    int nseg = q<=QF?1:nR;
    for(int sg=0;sg<nseg;sg++){
      long long lo,hi; if(q<=QF){lo=SYM?1:-(long long)(PMUL*q)-P0; hi=(long long)(PMUL*q)+P0;} else {long long c=llroundl((LD)R[sg]*q); lo=c-WS; hi=c+WS; if(SYM&&lo<1) lo=1;}
      if(lo>hi) continue;
      for(int k=0;k<n1;k++){ long long l=pr[k]; res[k]= iq[k]<0 ? (int)l : (int)(md(lo,l)*iq[k]%l); }
      for(long long p=lo;p<=hi;p++){
        if(p!=lo) for(int k=0;k<n1;k++){ if(iq[k]>=0){ int l=pr[k]; int r=res[k]+iq[k]; if(r>=l) r-=l; res[k]=r; } }
        tot++;
        if(p==0 && q!=1) continue;
        long long g0=p<0?-p:p, h0=q; while(h0){long long tt=g0%h0;g0=h0;h0=tt;} if(g0!=1) continue;
        LD a[5]; ainv_ld(p,q,a); LD h=Hld(a); if(h>HM) continue; c0++;
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
