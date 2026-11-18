// m23f: fast Mestre ten-root sieve (table lookups, no i128 in the hot loop).
// p(x)=prod(x-a_i)=q^2-r, quartic y^2=r(x); Jacobian E: y^2=x^3-27I x-27J (I,J invariants of r).
// For each prime l: table AP[l][I*l+J] = -a_l(E) (0 if singular). For a random base of 9 roots,
// per tabled prime l<=PB we build w_l[t] = -a_l * log l / l for the 10th root t mod l, then the
// S(PB) of every 10th root a in [0,M] is a sum of periodic lookups. Survivors of S(PB) continue
// prime by prime (base poly cached mod l) with checkpoints, then exact i128 height of the minimal
// model estimate (lscale), then S up to PMAX.
// All root sets are used (no parity filter): q is computed mod l with 1/2; for the exact height the
// roots are centred and scaled by 1,2,4 until q is integral (same curve up to isomorphism).
// env: M HMAX NB(bases) SEED PB C97 C199 C1000 T2 PMAX NPAIR G
// Build: cc -O3 -mcpu=neoverse-v1 -fopenmp -o m23f m23f.c -lm
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_thread_num(void){return 0;}
#include <time.h>
static double omp_get_wtime(void){struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec;}
#endif
typedef __int128 i128;
#define MAXP 1200
static int NPR=0, PR[MAXP]; static int8_t *AP[MAXP]; static double W[MAXP]; static uint64_t FM[MAXP];
static int IX97, IX199, IX1000, IXB;
static inline uint32_t fmod_(uint64_t x,int i){ uint64_t lowbits=FM[i]*x; return (uint32_t)(((unsigned __int128)lowbits*PR[i])>>64); } // x < 2^32

static int8_t chi(int64_t v,int l){ v%=l; if(v<0) v+=l; if(v==0) return 0; int64_t r=1,b=v,e=(l-1)/2; while(e){ if(e&1) r=r*b%l; b=b*b%l; e>>=1; } return r==1?1:-1; }
static void buildtab(int i){
  int l=PR[i]; int8_t *T=malloc((size_t)l*l); memset(T,127,(size_t)l*l); int8_t *ch=malloc(l); for(int x=0;x<l;x++) ch[x]=chi(x,l);
  int64_t u4[4096],u6[4096]; for(int u=1;u<l;u++){ int64_t a=(int64_t)u*u%l; a=a*a%l; u4[u]=a; u6[u]=a*u%l*u%l; }
  for(int I=0;I<l;I++) for(int J=0;J<l;J++){ if(T[I*l+J]!=127) continue;
    int64_t d=((4LL*I%l*I%l*I-(int64_t)J*J)%l+l)%l; int8_t v=0;
    if(d!=0){ int64_t A=((-27LL*I)%l+l)%l, B=((-27LL*J)%l+l)%l; int s=0; for(int64_t x=0;x<l;x++){ s+=ch[((x*x%l+A)*x+B)%l]; } v=(int8_t)s; } // -a_l = sum chi
    for(int u=1;u<l;u++){ int I2=(int)(u4[u]*I%l), J2=(int)(u6[u]*J%l); T[I2*l+J2]=v; } }
  AP[i]=T; free(ch); }

// poly mod l: c[0..deg]
static inline void mulroot(uint32_t *c,int deg,uint32_t na,int i){ // c*(x - a), na = l - a mod l
  int l=PR[i]; c[deg+1]=c[deg];
  for(int k=deg;k>=1;k--) c[k]=fmod_((uint64_t)c[k]*na + c[k-1],i);
  c[0]=fmod_((uint64_t)c[0]*na,i); (void)l; }
static inline int8_t val10(const uint32_t *p,int i){ // p monic deg 10 mod l -> -a_l
  int l=PR[i]; uint64_t h=(l+1)/2, L=l, L2=(uint64_t)l*l;
  #define R(x) fmod_((x),i)
  uint32_t q4=R(p[9]*h);
  uint32_t q3=R((p[8]+L2-(uint64_t)q4*q4)%L*h);
  uint32_t q2=R(R(p[7]+L2-2*(uint64_t)R((uint64_t)q4*q3))*h);
  uint32_t q1=R(R(p[6]+3*L2-2*(uint64_t)R((uint64_t)q4*q2)-(uint64_t)q3*q3)*h);
  uint32_t q0=R(R(p[5]+4*L2-2*(uint64_t)R((uint64_t)q4*q1)-2*(uint64_t)R((uint64_t)q3*q2))*h);
  uint32_t r4=R(2*(uint64_t)R((uint64_t)q4*q0)+2*(uint64_t)R((uint64_t)q3*q1)+(uint64_t)R((uint64_t)q2*q2)+L-p[4]);
  uint32_t r3=R(2*(uint64_t)R((uint64_t)q3*q0)+2*(uint64_t)R((uint64_t)q2*q1)+L-p[3]);
  uint32_t r2=R(2*(uint64_t)R((uint64_t)q2*q0)+(uint64_t)R((uint64_t)q1*q1)+L-p[2]);
  uint32_t r1=R(2*(uint64_t)R((uint64_t)q1*q0)+L-p[1]);
  uint32_t r0=R((uint64_t)R((uint64_t)q0*q0)+L-p[0]);
  uint32_t I=R(12*(uint64_t)R((uint64_t)r4*r0)+3*(L-(uint64_t)R((uint64_t)r3*r1))+R((uint64_t)r2*r2));
  uint32_t r22=R((uint64_t)r2*r2);
  uint64_t J= 72*(uint64_t)R((uint64_t)R((uint64_t)r4*r2)*r0) + 9*(uint64_t)R((uint64_t)R((uint64_t)r3*r2)*r1)
    + 27*(L-R((uint64_t)R((uint64_t)r4*r1)*r1)) + 27*(L-R((uint64_t)R((uint64_t)r0*r3)*r3)) + 2*(L-R((uint64_t)r22*r2));
  uint32_t Jm=(uint32_t)(J%L);
  #undef R
  (void)h; return AP[i][(size_t)I*l+Jm]; }

static long long mmod(i128 z,long long m){ long long r=(long long)(z%m); return r<0?r+m:r; }
static double lscale(const i128 *r){
  static const int L[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47}; double ls=0;
  for(int t=0;t<15;t++){ int l=L[t]; int e=0;
    while(1){ int ee=e+1; if(6.0*ee*log((double)l)>42.0) break;
      long long m4=1,m6=1; for(int k=0;k<4*ee;k++) m4*=l; for(int k=0;k<6*ee;k++) m6*=l; long long m=m6;
      i128 A=mmod(r[4],m),B=mmod(r[3],m),C=mmod(r[2],m),D=mmod(r[1],m),E=mmod(r[0],m);
      i128 I=((12*A%m)*E%m-(3*B%m)*D%m+C*C%m)%m; I=(I%m+m)%m;
      i128 J=((((72*A%m)*C%m)*E)%m+(((9*B%m)*C%m)*D)%m-(((27*A%m)*D%m)*D)%m-(((27*E%m)*B%m)*B)%m-(((2*C%m)*C%m)*C)%m)%m; J=(J%m+m)%m;
      i128 c4=(1296*I)%m, c6=(23328*J)%m; if(c4%m4!=0||c6%m6!=0) break; e=ee; }
    ls+=e*log((double)l); }
  return ls; }
// exact height estimate; returns 1e9 on failure
static double exactH(const int *a0){
  int a[10]; long s=0; for(int i=0;i<10;i++) s+=a0[i]; long c=lround((double)s/10);
  for(int sc=1;sc<=4;sc*=2){
    for(int i=0;i<10;i++) a[i]=(int)((a0[i]-c)*sc);
    i128 p[11]; for(int i=0;i<11;i++) p[i]=0; p[0]=1; int deg=0;
    for(int i=0;i<10;i++){ for(int k=deg+1;k>=1;k--) p[k]=p[k-1]-(i128)a[i]*p[k]; p[0]=-(i128)a[i]*p[0]; deg++; }
    i128 t; if(p[9]&1) continue; i128 q4=p[9]/2;
    t=p[8]-q4*q4; if(t&1) continue; i128 q3=t/2;
    t=p[7]-2*q4*q3; if(t&1) continue; i128 q2=t/2;
    t=p[6]-2*q4*q2-q3*q3; if(t&1) continue; i128 q1=t/2;
    t=p[5]-2*q4*q1-2*q3*q2; if(t&1) continue; i128 q0=t/2;
    i128 r[5]; r[4]=2*q4*q0+2*q3*q1+q2*q2-p[4]; r[3]=2*q3*q0+2*q2*q1-p[3]; r[2]=2*q2*q0+q1*q1-p[2]; r[1]=2*q1*q0-p[1]; r[0]=q0*q0-p[0];
    long double A=r[4],B=r[3],C=r[2],D=r[1],E=r[0];
    long double I=12*A*E-3*B*D+C*C, J=72*A*C*E+9*B*C*D-27*A*D*D-27*E*B*B-2*C*C*C;
    long double mx=fmaxl(3*logl(fabsl(1296.0L*I)+1),2*logl(fabsl(23328.0L*J)+1));
    return (double)(mx/12)-lscale(r); }
  return 1e9; }

static inline uint64_t xr(uint64_t *s){ *s^=*s<<13; *s^=*s>>7; *s^=*s<<17; return *s; }
static int cmpi(const void*x,const void*y){ return *(const int*)x-*(const int*)y; }


// ---- coordinate-descent climb (MODE=2) ----
static int CM; static double CHMAX;
static int nidx(int X){ int i=0; while(i<NPR && PR[i]<=X) i++; return i; }
// Sv[k] = S_{<n}(b U {lo+k}) for k in [0,len)
static void sweep(const int *b,int n,int lo,int len,float *Sv,float *wl){
  for(int k=0;k<len;k++) Sv[k]=0;
  for(int i=0;i<n;i++){ int l=PR[i]; float wt=(float)W[i]; uint32_t base[11]; base[0]=1; for(int j=0;j<9;j++) mulroot(base,j,(uint32_t)(((l-b[j]%l)%l+l)%l),i);
    for(int t=0;t<l;t++){ uint32_t c[11]; memcpy(c,base,sizeof(uint32_t)*10); mulroot(c,9,(uint32_t)((l-t)%l),i); wl[t]=wt*val10(c,i); }
    int t=((lo%l)+l)%l; for(int k=0;k<len;k++){ Sv[k]+=wl[t]; if(++t==l) t=0; } } }
static double Sfull(const int *a,int n){ double S=0; for(int i=0;i<n;i++){ int l=PR[i]; uint32_t c[11]; c[0]=1; for(int j=0;j<10;j++) mulroot(c,j,(uint32_t)(((l-a[j]%l)%l+l)%l),i); S+=W[i]*val10(c,i); } return S; }
static int normH(int *a,double *H){ qsort(a,10,sizeof(int),cmpi); for(int i=0;i<9;i++) if(a[i]==a[i+1]) return 0; int m=a[0]; for(int i=0;i<10;i++) a[i]-=m; if(a[9]>CM) return 0; *H=exactH(a); return *H<=CHMAX; }
static double climb1(int *cur,int n,uint64_t *st,float *Sv,float *wl,int *nsw){
  double curS=Sfull(cur,n); int stale=0, j=xr(st)%10;
  while(stale<10){ j=(j+1)%10; int b[9],k=0; for(int i=0;i<10;i++) if(i!=j) b[k++]=cur[i];
    int lo=b[8]-CM, hi=b[0]+CM, len=hi-lo+1; sweep(b,n,lo,len,Sv,wl); (*nsw)++;
    int found=0; for(int tries=0;tries<12 && !found;tries++){ int bk=-1; float bv=(float)curS+1e-4f; for(int q=0;q<len;q++) if(Sv[q]>bv){ bv=Sv[q]; bk=q; }
      if(bk<0) break; Sv[bk]=-1e30f; int a[10]; for(int i=0;i<9;i++) a[i]=b[i]; a[9]=lo+bk; double H; if(!normH(a,&H)) continue;
      memcpy(cur,a,sizeof a); curS=bv; found=1; }
    if(found) stale=0; else stale++; }
  return curS; }
static void climbmain(long long NB,uint64_t seed,int PX,int KICK,double T2,int NPAIR,int G){
  int n=nidx(PX); long long tot=0, nsw_all=0; double t0=omp_get_wtime();
  #pragma omp parallel reduction(+:tot,nsw_all)
  { uint64_t st=seed*0x9E3779B97F4A7C15ULL+(uint64_t)omp_get_thread_num()*0x632BE59BD9B4E019ULL+999; for(int w=0;w<20;w++) xr(&st);
    float *Sv=malloc(sizeof(float)*(2*CM+2)); float *wl=malloc(sizeof(float)*4096);
    #pragma omp for schedule(dynamic,1)
    for(long long it=0;it<NB;it++){
      int cur[10]; double H;
      for(;;){ int k=0; for(int pr=0;pr<NPAIR && k+1<10;pr++){ int L0=xr(&st)%(CM+1); int g=1+xr(&st)%G; if(L0+g>CM) L0=CM-g; cur[k++]=L0; cur[k++]=L0+g; } while(k<10) cur[k++]=xr(&st)%(CM+1); if(normH(cur,&H)) break; }
      int nsw=0; double best=climb1(cur,n,&st,Sv,wl,&nsw); int bestset[10]; memcpy(bestset,cur,sizeof cur);
      for(int kk=0;kk<KICK;kk++){ int c2[10]; memcpy(c2,bestset,sizeof c2); int ok=0;
        for(int tr=0;tr<50 && !ok;tr++){ memcpy(c2,bestset,sizeof c2); int nk=2+xr(&st)%2; for(int z=0;z<nk;z++) c2[xr(&st)%10]=xr(&st)%(CM+1); ok=normH(c2,&H); }
        if(!ok) continue; double v=climb1(c2,n,&st,Sv,wl,&nsw); if(v>best){ best=v; memcpy(bestset,c2,sizeof c2); } }
      tot++; nsw_all+=nsw;
      double Sf=Sfull(bestset,NPR); if(Sf<T2) continue; double Hh=exactH(bestset);
      double S199=Sfull(bestset,IX199), S1000=Sfull(bestset,IX1000);
      #pragma omp critical
      { for(int j=0;j<10;j++) printf("%d ",bestset[j]); printf("%.4f %.3f %.3f %.3f\n",Hh,S199,S1000,Sf); fflush(stdout); }
    }
    free(Sv); free(wl); }
  double dt=omp_get_wtime()-t0; fprintf(stderr,"climbs %lld sweeps %lld time %.2fs  %.3g climbs/s\n",tot,nsw_all,dt,tot/dt); }

// ---- 2D sweep (MODE=3): random base of 8 roots, all pairs (a9<a10) in [0,M] tabled per prime ----
static void sieve2d(int M,long long NB,uint64_t seed,int NPAIR,int G,double GATE,double C97,double C199,double C1000,double T2,double HMAX){
  long long cnt[4]={0}; double t0=omp_get_wtime(); int L1=M+1;
  #pragma omp parallel reduction(+:cnt[:4])
  { uint64_t st=seed*0x9E3779B97F4A7C15ULL+(uint64_t)omp_get_thread_num()*0x632BE59BD9B4E019ULL+777; for(int w=0;w<20;w++) xr(&st);
    uint32_t (*base)[9]=malloc(sizeof(uint32_t)*9*NPR); float *Sv=malloc(sizeof(float)*L1*L1); float *T=malloc(sizeof(float)*128*128); uint32_t (*c9)[10]=malloc(sizeof(uint32_t)*10*128);
    #pragma omp for schedule(dynamic,4)
    for(long long it=0;it<NB;it++){
      int b[8];
      for(;;){ int k=0;
        for(int pr=0;pr<NPAIR && k+1<8;pr++){ int L0=xr(&st)%(M+1); int g=1+xr(&st)%G; if(L0+g>M) L0=M-g; b[k++]=L0; b[k++]=L0+g; }
        while(k<8) b[k++]=xr(&st)%(M+1);
        qsort(b,8,sizeof(int),cmpi); int ok=1; for(int i=0;i<7;i++) if(b[i]==b[i+1]) ok=0; if(ok) break; }
      int have=0;
      #define NEED8(n) while(have<(n)){ int i_=have, l=PR[i_]; uint32_t c[11]; c[0]=1; for(int j=0;j<8;j++) mulroot(c,j,(uint32_t)((l-b[j]%l)%l),i_); memcpy(base[i_],c,sizeof(uint32_t)*9); have++; }
      NEED8(IXB);
      memset(Sv,0,sizeof(float)*L1*L1);
      for(int i=0;i<IXB;i++){ int l=PR[i]; float wt=(float)W[i];
        for(int t1=0;t1<l;t1++){ uint32_t *c=c9[t1]; memcpy(c,base[i],sizeof(uint32_t)*9); mulroot(c,8,(uint32_t)((l-t1)%l),i); }
        for(int t1=0;t1<l;t1++) for(int t2=t1;t2<l;t2++){ uint32_t c[11]; memcpy(c,c9[t1],sizeof(uint32_t)*10); mulroot(c,9,(uint32_t)((l-t2)%l),i); float v=wt*val10(c,i); T[t1*l+t2]=v; T[t2*l+t1]=v; }
        int r1=0; for(int x=0;x<=M;x++){ float *row=Sv+(size_t)x*L1; const float *tr=T+r1*l; int r2=(x+1)%l; for(int y=x+1;y<=M;y++){ row[y]+=tr[r2]; if(++r2==l) r2=0; } if(++r1==l) r1=0; } }
      for(int x=0;x<=M;x++){ int bad=0; for(int j=0;j<8;j++) if(b[j]==x) bad=1; if(bad) continue; float *row=Sv+(size_t)x*L1;
        for(int y=x+1;y<=M;y++){ cnt[0]++; if(row[y]<GATE) continue; int bad2=0; for(int j=0;j<8;j++) if(b[j]==y) bad2=1; if(bad2) continue;
          cnt[1]++; double S=row[y]; int okk=1;
          for(int i=IXB;i<NPR;i++){ if(i==IX97&&S<C97){okk=0;break;} if(i==IX199&&S<C199){okk=0;break;} if(i==IX1000&&S<C1000){okk=0;break;}
            NEED8(i+1); int l=PR[i]; uint32_t c[11]; memcpy(c,base[i],sizeof(uint32_t)*9); mulroot(c,8,(uint32_t)((l-x%l)%l),i); mulroot(c,9,(uint32_t)((l-y%l)%l),i); S+=W[i]*val10(c,i); }
          if(!okk||S<T2) continue; cnt[2]++;
          int s10[10]; for(int j=0;j<8;j++) s10[j]=b[j]; s10[8]=x; s10[9]=y; qsort(s10,10,sizeof(int),cmpi); int m0=s10[0]; for(int j=0;j<10;j++) s10[j]-=m0;
          double H=exactH(s10); if(H>HMAX) continue; cnt[3]++;
          double S199=0,S1000=0,Sf=0; for(int ii=0;ii<NPR;ii++){ int l=PR[ii]; uint32_t c[11]; c[0]=1; for(int j=0;j<10;j++) mulroot(c,j,(uint32_t)((l-s10[j]%l)%l),ii); Sf+=W[ii]*val10(c,ii); if(ii==IX199-1) S199=Sf; if(ii==IX1000-1) S1000=Sf; }
          #pragma omp critical
          { for(int j=0;j<10;j++) printf("%d ",s10[j]); printf("%.4f %.3f %.3f %.3f\n",H,S199,S1000,Sf); fflush(stdout); } } } }
    free(base); free(Sv); free(T); free(c9); }
  double dt=omp_get_wtime()-t0;
  fprintf(stderr,"2D cand %lld  passPB %lld  passS %lld  passH %lld  time %.2fs  rate %.3g cand/s\n",cnt[0],cnt[1],cnt[2],cnt[3],dt,cnt[0]/(dt>0?dt:1)); }

int main(void){
  #define GE(n,d) (getenv(n)?getenv(n):d)
  int M=atoi(GE("M","600")); double HMAX=atof(GE("HMAX","17.88")); long long NB=atoll(GE("NB","1000")); uint64_t seed=atoll(GE("SEED","1"));
  int PB=atoi(GE("PB","97")); double GATE=atof(GE("GATE","10")), C97=atof(GE("C97","14")), C199=atof(GE("C199","27")), C1000=atof(GE("C1000","56")), T2=atof(GE("T2","82"));
  int PMAX=atoi(GE("PMAX","3000")); int NPAIR=atoi(GE("NPAIR","0")); int G=atoi(GE("G","10")); int VERB=atoi(GE("VERB","0"));
  for(int l=5;l<=PMAX && NPR<MAXP;l++){ int ok=1; for(int d=2;d*d<=l;d++) if(l%d==0){ok=0;break;} if(!ok) continue; PR[NPR]=l; W[NPR]=log(l)/l; FM[NPR]=UINT64_C(0xFFFFFFFFFFFFFFFF)/l+1; NPR++; }
  for(IX97=0;IX97<NPR&&PR[IX97]<=97;IX97++); for(IX199=0;IX199<NPR&&PR[IX199]<=199;IX199++); for(IX1000=0;IX1000<NPR&&PR[IX1000]<=1000;IX1000++); for(IXB=0;IXB<NPR&&PR[IXB]<=PB;IXB++);
  double t0=omp_get_wtime();
  #pragma omp parallel for schedule(dynamic,1)
  for(int i=NPR-1;i>=0;i--) buildtab(i);
  if(VERB) fprintf(stderr,"tables %d primes %.1fs\n",NPR,omp_get_wtime()-t0);
  if(atoi(GE("MODE","1"))==2){ CM=M; CHMAX=HMAX; climbmain(NB,seed,atoi(GE("PX","500")),atoi(GE("KICK","0")),T2,NPAIR,G); return 0; }
  if(atoi(GE("MODE","1"))==3){ sieve2d(M,NB,seed,NPAIR,G,GATE,C97,C199,C1000,T2,HMAX); return 0; }
  long long cnt[6]={0}; t0=omp_get_wtime();
  #pragma omp parallel reduction(+:cnt[:6])
  { uint64_t st=seed*0x9E3779B97F4A7C15ULL+(uint64_t)omp_get_thread_num()*0x632BE59BD9B4E019ULL+12345; for(int w=0;w<20;w++) xr(&st);
    uint32_t (*base)[10]=malloc(sizeof(uint32_t)*10*NPR); float *Sv=malloc(sizeof(float)*(M+1)); float *wl=malloc(sizeof(float)*4096);
    #pragma omp for schedule(dynamic,16)
    for(long long it=0;it<NB;it++){
      int b[9];
      for(;;){ int k=0;
        for(int pr=0;pr<NPAIR && k+1<9;pr++){ int L0=xr(&st)%(M+1); int g=1+xr(&st)%G; if(L0+g>M) L0=M-g; b[k++]=L0; b[k++]=L0+g; }
        while(k<9) b[k++]=xr(&st)%(M+1);
        qsort(b,9,sizeof(int),cmpi); int ok=1; for(int i=0;i<8;i++) if(b[i]==b[i+1]) ok=0; if(ok) break; }
      // base polys mod l for all primes
      int have=0;
      #define NEED(n) while(have<(n)){ int i_=have, l=PR[i_]; uint32_t c[11]; c[0]=1; for(int j=0;j<9;j++) mulroot(c,j,(uint32_t)((l-b[j]%l)%l),i_); memcpy(base[i_],c,sizeof(uint32_t)*10); have++; }
      NEED(IXB);
      for(int a=0;a<=M;a++) Sv[a]=0;
      for(int i=0;i<IXB;i++){ int l=PR[i]; float wt=(float)W[i];
        for(int t=0;t<l;t++){ uint32_t c[11]; memcpy(c,base[i],sizeof(uint32_t)*10); mulroot(c,9,(uint32_t)((l-t)%l),i); wl[t]=wt*val10(c,i); }
        int t=0; for(int a=0;a<=M;a++){ Sv[a]+=wl[t]; if(++t==l) t=0; } }
      cnt[0]+=M+1-9;
      float gate=(float)GATE;
      for(int a=0;a<=M;a++){ if(Sv[a]<gate) continue;
        int dup=0; for(int j=0;j<9;j++) if(b[j]==a) dup=1; if(dup) continue;
        cnt[1]++;
        double S=Sv[a]; int i=IXB; int okk=1;
        for(;i<NPR;i++){ if(i==IX97 && S<C97){okk=0;break;} if(i==IX199 && S<C199){okk=0;break;} if(i==IX1000 && S<C1000){okk=0;break;}
          NEED(i+1); int l=PR[i]; uint32_t c[11]; memcpy(c,base[i],sizeof(uint32_t)*10); mulroot(c,9,(uint32_t)((l-a%l)%l),i); S+=W[i]*val10(c,i);
          if(i==IX199-1 && S<C199){okk=0;break;} }
        if(!okk) continue; if(S<T2) continue; cnt[2]++;
        int s10[10]; for(int j=0;j<9;j++) s10[j]=b[j]; s10[9]=a; qsort(s10,10,sizeof(int),cmpi); int m0=s10[0]; for(int j=0;j<10;j++) s10[j]-=m0;
        double H=exactH(s10); if(H>HMAX) continue; cnt[3]++;
        // recompute S199,S1000 for the record
        double S199=0,S1000=0,Sf=0; for(int ii=0;ii<NPR;ii++){ int l=PR[ii]; uint32_t c[11]; c[0]=1; for(int j=0;j<10;j++) mulroot(c,j,(uint32_t)((l-s10[j]%l)%l),ii); Sf+=W[ii]*val10(c,ii); if(ii==IX199-1) S199=Sf; if(ii==IX1000-1) S1000=Sf; }
        #pragma omp critical
        { for(int j=0;j<10;j++) printf("%d ",s10[j]); printf("%.4f %.3f %.3f %.3f\n",H,S199,S1000,Sf); fflush(stdout); }
      }
    }
    free(base); free(Sv); free(wl); }
  double dt=omp_get_wtime()-t0;
  fprintf(stderr,"cand %lld  passPB %lld  passS %lld  passH %lld  time %.2fs  rate %.3g cand/s\n",cnt[0],cnt[1],cnt[2],cnt[3],dt,cnt[0]/(dt>0?dt:1));
  return 0; }
