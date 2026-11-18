// m23: Mestre ten-root sieve, "small height first", with random and hill-climb modes.
// Derived from submission/search/mestre-own/msr.c (same quartic / lscale logic).
// p(x)=prod(x-a_i)=q^2-r, quartic y^2=r(x); its Jacobian is E (c4=1296 I, c6=23328 J).
// env: MODE (1 random, 3 hill-climb, 4 neighbourhood of stdin sets), M, HMAX, S1 (min S at p<=199),
//      T2 (min S at p<=3000 to print), NS (iterations), SEED, STEPS (hill-climb steps), XOBJ (objective prime bound)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_thread_num(void){return 0;}
#endif
typedef __int128 i128;
#define MAXP 1200
static int PR[MAXP], NPR=0, IDX199=0, IDX1000=0; static signed char *CHI[MAXP]; static double W[MAXP];
static int M=150; static double HMAX=17.88, S1=30.0, T2=85.0; static int XOBJ=1000;

static long long mm(i128 z,long long m){ long long r=(long long)(z%m); return r<0?r+m:r; }
static double lscale(i128 r4,i128 r3,i128 r2,i128 r1,i128 r0){
  static const int L[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
  double ls=0;
  for(int t=0;t<15;t++){ int l=L[t]; int e=0;
    while(1){ int ee=e+1; double lg=6.0*ee*log((double)l); if(lg>42.0) break;
      long long m4=1,m6=1; for(int k=0;k<4*ee;k++) m4*=l; for(int k=0;k<6*ee;k++) m6*=l;
      long long m=m6; i128 A=mm(r4,m),B=mm(r3,m),C=mm(r2,m),D=mm(r1,m),E=mm(r0,m);
      i128 I=( (12*A%m)*E%m - (3*B%m)*D%m + C*C%m )%m; I=(I%m+m)%m;
      i128 J=( (((72*A%m)*C%m)*E)%m + (((9*B%m)*C%m)*D)%m - (((27*A%m)*D%m)*D)%m - (((27*E%m)*B%m)*B)%m - (((2*C%m)*C%m)*C)%m )%m; J=(J%m+m)%m;
      i128 c4=(1296*I)%m, c6=(23328*J)%m;
      if(c4%m4!=0 || c6%m6!=0) break; e=ee; }
    ls+=e*log((double)l); }
  return ls; }

typedef struct { i128 r[5]; double H; int ok; } quart;
// build quartic from 10 roots; ok=0 if not integral / degenerate / H too large
static quart mkq(const int *a){
  quart Q; Q.ok=0;
  i128 p[11]; for(int i=0;i<11;i++) p[i]=0; p[0]=1; int deg=0;
  for(int i=0;i<10;i++){ for(int k=deg+1;k>=1;k--) p[k]=p[k-1]-(i128)a[i]*p[k]; p[0]=-(i128)a[i]*p[0]; deg++; }
  i128 t;
  if(p[9]&1) return Q; i128 q4=p[9]/2;
  t=p[8]-q4*q4; if(t&1) return Q; i128 q3=t/2;
  t=p[7]-2*q4*q3; if(t&1) return Q; i128 q2=t/2;
  t=p[6]-2*q4*q2-q3*q3; if(t&1) return Q; i128 q1=t/2;
  t=p[5]-2*q4*q1-2*q3*q2; if(t&1) return Q; i128 q0=t/2;
  i128 r4=2*q4*q0+2*q3*q1+q2*q2-p[4];
  i128 r3=2*q3*q0+2*q2*q1-p[3];
  i128 r2=2*q2*q0+q1*q1-p[2];
  i128 r1=2*q1*q0-p[1];
  i128 r0=q0*q0-p[0];
  if(r4==0) return Q;
  long double A=r4,B=r3,C=r2,D=r1,E=r0;
  long double I=12*A*E-3*B*D+C*C, J=72*A*C*E+9*B*C*D-27*A*D*D-27*E*B*B-2*C*C*C;
  long double disc=4*I*I*I-J*J; if(fabsl(disc)<1e-30L*(fabsl(4*I*I*I)+J*J+1)) return Q;
  long double c4=1296.0L*I, c6=23328.0L*J;
  long double mx=fmaxl(3*logl(fabsl(c4)+1),2*logl(fabsl(c6)+1)); double H=(double)(mx/12);
  if(H>HMAX+6.0) return Q;
  H-=lscale(r4,r3,r2,r1,r0);
  if(H>HMAX) return Q;
  Q.r[0]=r0;Q.r[1]=r1;Q.r[2]=r2;Q.r[3]=r3;Q.r[4]=r4; Q.H=H; Q.ok=1; return Q; }

// -a_p for prime index i (0 if bad)
static long long negap(const quart *Q,int i){
  int pp=PR[i]; long long R[5]; for(int k=0;k<5;k++){ R[k]=(long long)(Q->r[k]%pp); if(R[k]<0) R[k]+=pp; }
  long long R4=R[4],R3=R[3],R2=R[2],R1=R[1],R0=R[0];
  if(R4==0) return 0;
  { long long Ip=((12*R4%pp*R0 - 3*R3%pp*R1 + R2*R2)%pp+pp)%pp;
    long long Jp=((72*R4%pp*R2%pp*R0 + 9*R3%pp*R2%pp*R1 - 27*R4%pp*R1%pp*R1 - 27*R0%pp*R3%pp*R3 - 2*R2%pp*R2%pp*R2)%pp+pp)%pp;
    long long Dp=((4*Ip%pp*Ip%pp*Ip - Jp*Jp)%pp+pp)%pp; if(Dp==0) return 0; }
  long long s=0, F[5]; for(int k=0;k<5;k++){ long long x=k; F[k]=((((R4*x+R3)%pp*x+R2)%pp*x+R1)%pp*x+R0)%pp; }
  long long d0=F[0], d1=(F[1]-F[0]+pp)%pp, d2=((F[2]-2*F[1]+F[0])%pp+2*pp)%pp, d3=((F[3]-3*F[2]+3*F[1]-F[0])%pp+4*pp)%pp, d4=((F[4]-4*F[3]+6*F[2]-4*F[1]+F[0])%pp+8*pp)%pp;
  const signed char *ch=CHI[i];
  for(long long x=0;x<pp;x++){ s+=ch[d0]; d0+=d1; if(d0>=pp)d0-=pp; d1+=d2; if(d1>=pp)d1-=pp; d2+=d3; if(d2>=pp)d2-=pp; d3+=d4; if(d3>=pp)d3-=pp; }
  return s+CHI[i][R4]; }

// S up to prime index n, with early abort if below thresholds (return -1e9)
static double Sto(const quart *Q,int n,double s1){
  double S=0; for(int i=0;i<n;i++){ if(i==IDX199 && S<s1) return -1e9; S+=negap(Q,i)*W[i];
    if(i==12 && S<s1*0.296) return -1e9; if(i==24 && S<s1*0.61) return -1e9; }
  return S; }

static int cmpi(const void*x,const void*y){ return *(const int*)x-*(const int*)y; }
static int normset(int *a){ // sort, translate min to 0, reject repeats; also reflect so that a canonical orientation
  qsort(a,10,sizeof(int),cmpi); for(int i=0;i<9;i++) if(a[i]==a[i+1]) return 0;
  int m=a[0]; for(int i=0;i<10;i++) a[i]-=m; return 1; }

static void report(const int *a,const quart *Q){
  double S199=Sto(Q,IDX199,-1e9); if(S199<S1) return;
  double S=Sto(Q,NPR,-1e9); if(S<T2) return;
  #pragma omp critical
  { for(int i=0;i<10;i++) printf("%d ",a[i]); printf("%.4f %.3f %.3f\n",Q->H,S199,S); fflush(stdout); } }

static unsigned long long xr(unsigned long long *s){ *s^=*s<<13; *s^=*s>>7; *s^=*s<<17; return *s; }
static int NDM=0, DM[8], NCV=0, CV[8];
static int congok(const int *a){ // distinct residues mod DM[], all residues covered mod CV[]
  for(int t=0;t<NDM;t++){ int m=DM[t]; unsigned long long seen=0; for(int i=0;i<10;i++){ int r=((a[i]%m)+m)%m; if(seen>>r&1) return 0; seen|=1ULL<<r; } }
  for(int t=0;t<NCV;t++){ int m=CV[t]; unsigned long long seen=0; for(int i=0;i<10;i++){ int r=((a[i]%m)+m)%m; seen|=1ULL<<r; } if(seen!=(1ULL<<m)-1) return 0; }
  return 1; }
static void parselist(const char *s,int *v,int *n){ *n=0; if(!s) return; const char *p=s; while(*p){ v[(*n)++]=atoi(p); while(*p && *p!=',') p++; if(*p) p++; } }
static int PAIRS=0, G=10;
static void pairset(int *a,unsigned long long *st){ // 5 pairs (L_i, L_i+g_i), L_1=0, spread <= M
  for(;;){ int L[5],g[5],SP=M/4; for(int i=0;i<5;i++) g[i]=1+xr(st)%G;
    L[0]=0; for(int i=1;i<5;i++){ int sp=1+xr(st)%SP; L[i]=L[i-1]+g[i-1]+sp; }
    if(L[4]+g[4]>M) continue; for(int i=0;i<5;i++){ a[2*i]=L[i]; a[2*i+1]=L[i]+g[i]; }
    if(congok(a) && normset(a)) return; } }
static void randset(int *a,unsigned long long *st){
  if(PAIRS){ pairset(a,st); return; }
  for(;;){ a[0]=0; int k=1, guard=0;
    while(k<10 && guard<5000){ guard++; int c=1+xr(st)%M; int bad=0; for(int j=0;j<k && !bad;j++){ if(a[j]==c) bad=1; for(int t=0;t<NDM && !bad;t++) if((c-a[j])%DM[t]==0) bad=1; } if(!bad) a[k++]=c; }
    if(k<10) continue; if(!congok(a)) continue; if(normset(a)) return; } }

int main(int argc,char**argv){
  (void)argc;(void)argv;
  #define GE(n,d) (getenv(n)?getenv(n):d)
  int MODE=atoi(GE("MODE","1")); M=atoi(GE("M","150")); HMAX=atof(GE("HMAX","17.88")); S1=atof(GE("S1","30"));
  T2=atof(GE("T2","85")); long long NS=atoll(GE("NS","100000")); unsigned long long seed=atoll(GE("SEED","1"));
  int STEPS=atoi(GE("STEPS","200")); parselist(getenv("DM"),DM,&NDM); PAIRS=atoi(GE("PAIRS","0")); G=atoi(GE("G","10")); parselist(getenv("CV"),CV,&NCV); XOBJ=atoi(GE("XOBJ","1000")); int PMAX=atoi(GE("PMAX","3000"));
  for(int pp=5;pp<=PMAX && NPR<MAXP;pp++){int ok=1;for(int d=2;d*d<=pp;d++) if(pp%d==0){ok=0;break;} if(!ok) continue; PR[NPR]=pp; CHI[NPR]=malloc(pp); for(int x=0;x<pp;x++) CHI[NPR][x]=-1; CHI[NPR][0]=0; for(long long x=1;x<pp;x++) CHI[NPR][x*x%pp]=1; W[NPR]=log(pp)/pp; NPR++;}
  for(IDX199=0;IDX199<NPR && PR[IDX199]<=199;IDX199++);
  for(IDX1000=0;IDX1000<NPR && PR[IDX1000]<=XOBJ;IDX1000++);
  if(MODE==4){ // neighbourhood of stdin sets: all single-root moves within +-R and pair moves
    int R=atoi(GE("R","6")); int a0[10]; static int sets[100000][10]; int ns=0;
    while(ns<100000){ int k=0; for(k=0;k<10;k++) if(scanf("%d",&a0[k])!=1) break; if(k<10) break; char buf[512]; if(!fgets(buf,512,stdin)) {} memcpy(sets[ns++],a0,sizeof a0); }
    #pragma omp parallel for schedule(dynamic,1)
    for(int s=0;s<ns;s++){ for(int i=0;i<10;i++) for(int d=-R;d<=R;d++){ if(!d) continue; int a[10]; memcpy(a,sets[s],sizeof a); a[i]+=d; if(!normset(a)) continue;
        quart Q=mkq(a); if(!Q.ok) continue; report(a,&Q); } }
    return 0; }
  #pragma omp parallel
  { unsigned long long st=seed*0x9E3779B97F4A7C15ULL+(unsigned long long)omp_get_thread_num()*7919+1; for(int w=0;w<10;w++) xr(&st);
    #pragma omp for schedule(dynamic,100)
    for(long long it=0;it<NS;it++){
      int a[10];
      if(MODE==1){ randset(a,&st); quart Q=mkq(a); if(!Q.ok) continue; report(a,&Q); continue; }
      if(MODE==3){ // hill-climb on S(XOBJ) subject to H<=HMAX
        quart Q; int tries=0; do{ randset(a,&st); Q=mkq(a); tries++; }while(!Q.ok && tries<1000); if(!Q.ok) continue;
        double best=Sto(&Q,IDX1000,-1e9); int stall=0;
        for(int s=0;s<50*STEPS && stall<STEPS;s++){ int b[10]; memcpy(b,a,sizeof b); int i=xr(&st)%10;
          int r=xr(&st)%4; if(r==0) b[i]=xr(&st)%(M+1); else { int d=1+xr(&st)%5; b[i]+= (xr(&st)&1)?d:-d; }
          if(!normset(b) || b[9]>M || !congok(b)) { continue; } quart Q2=mkq(b); if(!Q2.ok){ continue;}
          double v=Sto(&Q2,IDX1000,-1e9); if(v>=best){ if(v>best) stall=0; else stall++; best=v; memcpy(a,b,sizeof b); Q=Q2; } else stall++; }
        report(a,&Q); continue; }
    }
  }
  return 0; }
