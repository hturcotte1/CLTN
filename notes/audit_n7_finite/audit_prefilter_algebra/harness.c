#include "filter_instr.c"
/* input lines: <d6> nbeta p1 q1 p2 q2 ... ; output: per i: coefs of m!*a_i(2+u) low->high ; then a_eval(i,p,q) for each sample ; plus domfree flag, final flag of filter logic */
static void pr128(i128 v){ char buf[64]; int k=0; int neg=v<0; unsigned __int128 u = neg? -(unsigned __int128)v : (unsigned __int128)v; if(!u) buf[k++]='0'; while(u){buf[k++]='0'+(int)(u%10); u/=10;} if(neg) putchar('-'); while(k) putchar(buf[--k]); }
int main(void){
  char line[4096];
  while(fgets(line,sizeof line,stdin)){
    char d6[64]; int off=0; if(sscanf(line,"%63s%n",d6,&off)!=1) continue;
    decode_d6(d6);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) N[i][j] = (i != j && A[j][i] == 0) ? 1 : 0;
    printf("%s df=%d", d6, domination_free());
    int drop=0;
    for(int i=0;i<n;i++){ int r=descartes_nonpositive(i); if(r) drop=1; printf(" #c%d:",i); for(int j=0;j<=n;j++){putchar(' '); pr128(G_coef[j]);} printf(" np=%d",r);} 
    printf(" #drop=%d", drop);
    char *p=line+off; int cnt; long long pp,qq; int o2;
    if(sscanf(p,"%d%n",&cnt,&o2)==1){ p+=o2; printf(" #ev"); for(int s=0;s<cnt;s++){ sscanf(p,"%lld %lld%n",&pp,&qq,&o2); p+=o2; for(int i=0;i<n;i++){putchar(' '); pr128(a_eval(i,pp,qq));} printf(" ;"); } }
    putchar('\n');
  }
  return 0;
}
