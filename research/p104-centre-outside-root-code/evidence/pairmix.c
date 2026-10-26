/* One mixed rotation replacing the pure rotations of directions j, k (R_8 = 0 profile).
   P_j subset h_j^perp, P_k subset h_k^perp meeting in a line l; Pq one of the 5 other planes of P_j + P_k.
   Mj = N_j \ W(P_j), Mk = N_k \ W(P_k).  Need Mj, Mk inside W(Pq); t in Mj sees P_k or is in Mk; t in Mk sees P_j or is in Mj;
   all double-missers share one ratio r_tj : r_tk.  Count solutions; report best partial (fewest failing ports). */
#include <stdio.h>
#include <stdint.h>
typedef unsigned __int128 u128;
int P, Q[128], H[8], R[128][8]; u128 pm[512], N[8];
static int par(int x){ return __builtin_parity(x); }
static int popc(u128 x){ return __builtin_popcountll((uint64_t)x) + __builtin_popcountll((uint64_t)(x>>64)); }
static u128 W2(int a, int b){ return pm[a] & pm[b]; }
int main(){
  FILE*f=fopen("mixdata.txt","r"); if(fscanf(f,"%d",&P)!=1) return 1;
  for(int i=0;i<P;i++) if(fscanf(f,"%d",&Q[i])!=1) return 1;
  for(int j=0;j<8;j++) if(fscanf(f,"%d",&H[j])!=1) return 1;
  for(int t=0;t<P;t++) for(int j=0;j<8;j++) if(fscanf(f,"%d",&R[t][j])!=1) return 1;
  for(int w=0;w<512;w++){ pm[w]=0; for(int p=0;p<P;p++) if(!par(Q[p]&w)) pm[w]|=((u128)1)<<p; }
  for(int j=0;j<8;j++){ N[j]=0; for(int t=0;t<P;t++) if(R[t][j]) N[j]|=((u128)1)<<t; }
  long long sols=0, configs=0; int bestfail=1000;
  for(int j=0;j<8;j++) for(int k=j+1;k<8;k++){
    long long ps=0; int pbest=1000;
    for(int l=1;l<512;l++){ if(par(l&H[j])||par(l&H[k])) continue;
      for(int a=1;a<512;a++){ if(a==l||par(a&H[j])) continue; if((a^l)<a) continue;   /* plane {l,a,a^l} once */
        u128 WPj=W2(l,a); u128 Mj=N[j]&~WPj;
        for(int b=1;b<512;b++){ if(b==l||par(b&H[k])) continue; if((b^l)<b) continue;
          if(b==a||b==(a^l)) continue;                 /* P_k != P_j */
          u128 WPk=W2(l,b); u128 Mk=N[k]&~WPk;
          /* F = span(l,a,b), 7 planes; P_j={l,a}, P_k={l,b}; others: {a,b},{a,b^l},{l^a... } enumerate via pairs */
          int Fv[7]={l,a,a^l,b,b^l,a^b,a^b^l};
          for(int x=0;x<7;x++) for(int y=x+1;y<7;y++){ int u=Fv[x], v=Fv[y]; int w=u^v;
            /* plane {u,v,w}; skip P_j and P_k and duplicates (require u<v<w ordering by index trick) */
            int z=-1; for(int q=0;q<7;q++) if(Fv[q]==w) z=q; if(z<y) continue;
            int isPj=((u==l||v==l||w==l)&&(u==a||v==a||w==a)); int isPk=((u==l||v==l||w==l)&&(u==b||v==b||w==b));
            if(isPj||isPk) continue;
            configs++;
            u128 Wq=W2(u,v);
            u128 fail=(Mj&~Wq)|(Mk&~Wq);
            fail |= Mj & ~WPk & ~Mk;     /* misses j, cannot separate (does not see P_k, not a double-misser) */
            fail |= Mk & ~WPj & ~Mj;
            u128 D=Mj&Mk; int nf=popc(fail);
            /* ratio consistency among double-missers */
            int ra=0, rb=0, bad=0;
            for(int t=0;t<P;t++) if((D>>t)&1){ if(!ra&&!rb){ ra=R[t][j]; rb=R[t][k]; } else if(ra*R[t][k]!=rb*R[t][j]) bad++; }
            nf+=bad;
            if(nf<pbest) pbest=nf;
            if(nf==0){ ps++; if(ps<=3) printf("SOLUTION j=%d k=%d l=%d a=%d b=%d Pq={%d,%d,%d} |Mj|=%d |Mk|=%d |D|=%d\n",j,k,l,a,b,u,v,w,popc(Mj),popc(Mk),popc(D)); }
          }
        }
      }
    }
    printf("pair %d %d: solutions %lld, fewest failing ports %d\n", j, k, ps, pbest); fflush(stdout);
    sols+=ps; if(pbest<bestfail) bestfail=pbest;
  }
  printf("configs %lld solutions %lld best %d\n", configs, sols, bestfail);
  return 0;
}
