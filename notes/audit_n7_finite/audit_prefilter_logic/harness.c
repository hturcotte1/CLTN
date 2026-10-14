#define main filter_main
#include "filter7p_copy.c"
#undef main
/* prints: d6 domfree | coefs per i (m! a_i(2+u), low->high) | a_eval(i,p,64) for listed p | filter verdict */
static void pr(i128 v){ char b[64]; int k=0; int neg=v<0; unsigned __int128 u = neg? -(unsigned __int128)v : (unsigned __int128)v; if(!u) b[k++]='0'; while(u){b[k++]='0'+(int)(u%10); u/=10;} if(neg) putchar('-'); while(k) putchar(b[--k]); }
int main(int argc, char **argv){
  char line[256]; int ps[64]; int np=0;
  for (int a=1;a<argc;a++) ps[np++]=atoi(argv[a]);
  while (fgets(line,sizeof line,stdin)){
    size_t L=strlen(line); while(L&&(line[L-1]=='\n'||line[L-1]=='\r')) line[--L]=0; if(!L) continue;
    char *sp=strchr(line,' '); if(sp)*sp=0;
    if (decode_d6(line)) { fprintf(stderr,"bad\n"); return 2; }
    int df = domination_free();
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) N[i][j]=(i!=j && A[j][i]==0)?1:0;
    printf("%s %d", line, df);
    int drop=0;
    for (int i=0;i<n;i++){ int r=descartes_nonpositive(i); if(r) drop=1; printf(" |"); for(int j=0;j<=n;j++){ putchar(' '); pr(g_coef[j]); } }
    printf(" | E");
    for (int s=0;s<np;s++) for (int i=0;i<n;i++){ putchar(' '); pr(a_eval(i, ps[s], 64)); }
    char flag='U';
    if (!df) flag='D'; else if (drop) flag='N'; else {
      for (int s=0;s<NSAMPLES && flag=='U';s++){ i128 p=2*SAMPLES[s][1]+SAMPLES[s][0], q=SAMPLES[s][1]; int allpos=1; for(int i=0;i<n;i++) if(a_eval(i,p,q)<=0){allpos=0;break;} if(allpos) flag='P'; }
    }
    printf(" | %c\n", flag);
  }
  return 0;
}
