// m12f: sieve for the twelve-root Mestre variant p(x) = F(x) F(x+t), F = prod(x - b_i), sum b = 0, e2 e3 = 2 e5.
// For such F, r = q^2 - p has degree <= 4 for every t (r5 = t^2 (s2 s3 - 2 s5) = 0), so y^2 = r(x) is a genus-one
// quartic with the 12 rational points (root, q(root)). Parameter t = u/v: roots v b_i and v b_i - u.
// Same table method as ../m23f.c: AP[l][I*l+J] = -a_l of the Jacobian y^2 = x^3 - 27 I x - 27 J.
// stdin: lines "b1 .. b6" (from fsol). env: V (max v), UMAX, PB, GATE, C97, C199, C1000, T2, PMAX
// out: "b1..b6 v u S199 S1000 S3000"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
static double omp_get_wtime(void){struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec;}
#endif
#define MAXP 1200
static int NPR=0, PR[MAXP]; static int8_t *AP[MAXP]; static double W[MAXP]; static uint64_t FM[MAXP];
static int IX97, IX199, IX1000, IXB;
static inline uint32_t fmod_(uint64_t x,int i){ uint64_t lowbits=FM[i]*x; return (uint32_t)(((unsigned __int128)lowbits*PR[i])>>64); }
static int8_t chi(int64_t v,int l){ v%=l; if(v<0) v+=l; if(v==0) return 0; int64_t r=1,b=v,e=(l-1)/2; while(e){ if(e&1) r=r*b%l; b=b*b%l; e>>=1; } return r==1?1:-1; }
static void buildtab(int i){
  int l=PR[i]; int8_t *T=malloc((size_t)l*l); memset(T,127,(size_t)l*l); int8_t *ch=malloc(l); for(int x=0;x<l;x++) ch[x]=chi(x,l);
  int64_t u4[4096],u6[4096]; for(int u=1;u<l;u++){ int64_t a=(int64_t)u*u%l; a=a*a%l; u4[u]=a; u6[u]=a*u%l*u%l; }
  for(int I=0;I<l;I++) for(int J=0;J<l;J++){ if(T[I*l+J]!=127) continue;
    int64_t d=((4LL*I%l*I%l*I-(int64_t)J*J)%l+l)%l; int8_t v=0;
    if(d!=0){ int64_t A=((-27LL*I)%l+l)%l, B=((-27LL*J)%l+l)%l; int s=0; for(int64_t x=0;x<l;x++){ s+=ch[((x*x%l+A)*x+B)%l]; } v=(int8_t)s; } // -a_l = sum chi
    for(int u=1;u<l;u++){ int I2=(int)(u4[u]*I%l), J2=(int)(u6[u]*J%l); T[I2*l+J2]=v; } }
  AP[i]=T; free(ch); }
static inline void mulroot(uint32_t *c,int deg,uint32_t na,int i){ c[deg+1]=c[deg]; for(int k=deg;k>=1;k--) c[k]=fmod_((uint64_t)c[k]*na + c[k-1],i); c[0]=fmod_((uint64_t)c[0]*na,i); }
// p monic degree 12 mod l (p[0..11]) -> -a_l of Jac(y^2 = r), r = q^2 - p
static inline int8_t val12(const uint32_t *p,int i){
  int l=PR[i]; uint64_t L=l, h=(l+1)/2; uint32_t q[7]; q[6]=1;
  for(int m=5;m>=0;m--){ uint64_t s=0; for(int j=m+1;j<=5;j++){ int k=6+m-j; if(k<=m||k>5) continue; s+=(uint64_t)q[j]*q[k]; }
    s=fmod_(s,i); q[m]=fmod_((uint64_t)fmod_(p[6+m]+L-s,i)*h,i); }
  uint32_t r[5];
  for(int k=0;k<=4;k++){ uint64_t s=0; for(int a=0;a<=k;a++) s+=(uint64_t)q[a]*q[k-a]; s=fmod_(fmod_(s,i)+L-p[k],i); r[k]=(uint32_t)s; }
  #define R(x) fmod_((x),i)
  uint32_t r4=r[4],r3=r[3],r2=r[2],r1=r[1],r0=r[0];
  uint32_t I=R(12*(uint64_t)R((uint64_t)r4*r0)+3*(L-(uint64_t)R((uint64_t)r3*r1))+R((uint64_t)r2*r2));
  uint32_t r22=R((uint64_t)r2*r2);
  uint64_t J= 72*(uint64_t)R((uint64_t)R((uint64_t)r4*r2)*r0) + 9*(uint64_t)R((uint64_t)R((uint64_t)r3*r2)*r1)
    + 27*(L-R((uint64_t)R((uint64_t)r4*r1)*r1)) + 27*(L-R((uint64_t)R((uint64_t)r0*r3)*r3)) + 2*(L-R((uint64_t)r22*r2));
  #undef R
  return AP[i][(size_t)I*l+(uint32_t)(J%L)]; }
static inline uint32_t md(long x,int l){ long r=x%l; return (uint32_t)(r<0?r+l:r); }
static int gcdi(int a,int b){ while(b){int t=a%b;a=b;b=t;} return a<0?-a:a; }
int main(void){
  #define GE(n,d) (getenv(n)?getenv(n):d)
  int V=atoi(GE("V","1")), V0=atoi(GE("V0","1")); int UMAX=atoi(GE("UMAX","2000")); int PB=atoi(GE("PB","47"));
  double GATE=atof(GE("GATE","10.5")), C97=atof(GE("C97","16")), C199=atof(GE("C199","27")), C1000=atof(GE("C1000","56")), T2=atof(GE("T2","80"));
  int PMAX=atoi(GE("PMAX","3000")); int VERB=atoi(GE("VERB","0"));
  for(int l=5;l<=PMAX && NPR<MAXP;l++){ int ok=1; for(int d=2;d*d<=l;d++) if(l%d==0){ok=0;break;} if(!ok) continue; PR[NPR]=l; W[NPR]=log(l)/l; FM[NPR]=UINT64_C(0xFFFFFFFFFFFFFFFF)/l+1; NPR++; }
  for(IX97=0;IX97<NPR&&PR[IX97]<=97;IX97++); for(IX199=0;IX199<NPR&&PR[IX199]<=199;IX199++); for(IX1000=0;IX1000<NPR&&PR[IX1000]<=1000;IX1000++); for(IXB=0;IXB<NPR&&PR[IXB]<=PB;IXB++);
  double t0=omp_get_wtime();
  #pragma omp parallel for schedule(dynamic,1)
  for(int i=NPR-1;i>=0;i--) buildtab(i);
  if(VERB) fprintf(stderr,"tables %d primes %.1fs\n",NPR,omp_get_wtime()-t0);
  static int B[2000000][6]; int nf=0; while(nf<2000000){ int k; for(k=0;k<6;k++) if(scanf("%d",&B[nf][k])!=1) break; if(k<6) break; nf++; }
  long long cnt[4]={0}; t0=omp_get_wtime(); long long njobs=(long long)nf*(V-V0+1);
  #pragma omp parallel reduction(+:cnt[:4])
  { uint32_t (*base)[7]=malloc(sizeof(uint32_t)*7*NPR); float *Sv=malloc(sizeof(float)*(UMAX+1)); float *wl=malloc(sizeof(float)*4096);
    #pragma omp for schedule(dynamic,1)
    for(long long job=0;job<njobs;job++){
      int f=(int)(job%nf), v=V0+(int)(job/nf); const int *b=B[f]; int vb[6]; for(int k=0;k<6;k++) vb[k]=v*b[k];
      int have=0;
      #define NEED(n) while(have<(n)){ int i_=have, l=PR[i_]; uint32_t c[13]; c[0]=1; for(int j=0;j<6;j++) mulroot(c,j,md(-vb[j],l),i_); memcpy(base[i_],c,sizeof(uint32_t)*7); have++; }
      NEED(IXB);
      for(int u=0;u<=UMAX;u++) Sv[u]=0;
      for(int i=0;i<IXB;i++){ int l=PR[i]; float wt=(float)W[i];
        for(int t=0;t<l;t++){ uint32_t c[13]; memcpy(c,base[i],sizeof(uint32_t)*7); for(int j=0;j<6;j++) mulroot(c,6+j,md(t-vb[j],l),i); wl[t]=wt*val12(c,i); }
        int t=0; for(int u=0;u<=UMAX;u++){ Sv[u]+=wl[t]; if(++t==l) t=0; } }
      cnt[0]+=UMAX;
      for(int u=1;u<=UMAX;u++){ if(Sv[u]<GATE) continue; if(gcdi(u,v)!=1) continue;
        int dup=0; for(int j=0;j<6&&!dup;j++) for(int k=0;k<6;k++) if(vb[j]-u==vb[k]) {dup=1;break;} if(dup) continue;
        cnt[1]++; double S=Sv[u]; int okk=1;
        for(int i=IXB;i<NPR;i++){ if(i==IX97&&S<C97){okk=0;break;} if(i==IX199&&S<C199){okk=0;break;} if(i==IX1000&&S<C1000){okk=0;break;}
          NEED(i+1); int l=PR[i]; uint32_t c[13]; memcpy(c,base[i],sizeof(uint32_t)*7); for(int j=0;j<6;j++) mulroot(c,6+j,md((long)u-vb[j],l),i); S+=W[i]*val12(c,i); }
        if(!okk||S<T2) continue; cnt[2]++;
        double S199=0,S1000=0,Sf=0; for(int ii=0;ii<NPR;ii++){ int l=PR[ii]; uint32_t c[13]; c[0]=1; for(int j=0;j<6;j++) mulroot(c,j,md(-vb[j],l),ii); for(int j=0;j<6;j++) mulroot(c,6+j,md((long)u-vb[j],l),ii); Sf+=W[ii]*val12(c,ii); if(ii==IX199-1) S199=Sf; if(ii==IX1000-1) S1000=Sf; }
        #pragma omp critical
        { printf("%d %d %d %d %d %d %d %d %.3f %.3f %.3f\n",b[0],b[1],b[2],b[3],b[4],b[5],v,u,S199,S1000,Sf); fflush(stdout); }
      } }
    free(base); free(Sv); free(wl); }
  double dt=omp_get_wtime()-t0;
  fprintf(stderr,"cand %lld passPB %lld passS %lld time %.2fs rate %.3g cand/s\n",cnt[0],cnt[1],cnt[2],dt,cnt[0]/(dt>0?dt:1));
  return 0; }
