/*
 * orbitenum.c -- enumerate, one representative per isomorphism class, the digraphs whose underlying simple
 * graph is a given graph G (every edge of G carries one arc or both) with at least min_arcs arcs.
 *
 *   orbitenum <graph6 G | @file> <min_arcs> [p/P] [-u]       (@file: G is the first line of file)
 *
 * Method (notes/k8_decomposition.md, relative form): a digraph on G is a pair (B, O): B = the set of edges
 * of G carrying both arcs, O = an orientation of G - B.  Two such digraphs are isomorphic iff some
 * automorphism of G maps one to the other, so the classes are in bijection with the pairs
 *   (B up to Aut(G),  orbit of Stab(B) = {g in Aut(G) : g(B) = B} on the orientations of G - B).
 * Aut(G) is found by brute force over the n! vertex permutations.  A subset B is a representative iff it is
 * the numerically least member of its Aut(G)-orbit (edge bitmask with the edges of G in lexicographic
 * order); an orientation O of G - B is a representative iff it is the least member of its Stab(B)-orbit
 * (bit t = 1 iff the t-th free edge (u,v), u < v, is oriented v -> u).  Arcs = |E(G)| + |B|, so only
 * |B| >= min_arcs - |E(G)| is enumerated.  Representatives B are numbered in increasing bitmask order and
 * p/P keeps those with index == p (mod P), so the parts partition the output.  -u counts without output.
 * Output: digraph6 (nauty: '&', N(n), n*n bits row-major, bit (i,j) = arc i -> j).  Stats on stderr:
 * "orbitenum <G> min_arcs <a> part <p/P> aut <|Aut(G)|> subsets <k> total <N>" (k = representatives B in THIS
 * part, N = digraphs in this part) plus one "B <mask> <count>" line per B of this part, in enumeration order.
 * Integers only; no isomorphism tests beyond the explicit group action.  Refuses |Aut(G)| > MAXAUT.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 10
#define MAXE 45
#define MAXAUT 8192
typedef uint64_t u64;

static int n, m;
static int eu[MAXE], ev[MAXE];          /* edge k = {eu[k], ev[k]}, eu < ev, lexicographic order */
static int pidx[MAXN][MAXN];           /* pair -> edge index or -1 */
static unsigned char adj[MAXN][MAXN];
static int naut;
static unsigned char perm[MAXAUT][MAXN];
static unsigned char emap[MAXAUT][MAXE];   /* image of edge k under automorphism g */
static unsigned char eflip[MAXAUT][MAXE];  /* 1 iff g reverses the (u<v) order of edge k */
static u64 *tbl;                           /* [g][byte][256] -> image mask of that byte */
static int nbytes;

static int decode_g6(const char *s) {
    int nn = s[0] - 63; if (nn < 1 || nn > MAXN) return -1;
    memset(adj, 0, sizeof adj); const char *p = s + 1; int k = 0;
    for (int j = 1; j < nn; j++) for (int i = 0; i < j; i++) {
        int c = (unsigned char)p[k / 6]; if (c < 63 || c > 126) return -1;
        adj[i][j] = adj[j][i] = (unsigned char)(((c - 63) >> (5 - k % 6)) & 1); k++; }
    return nn;
}
static int is_aut(const unsigned char *p) {
    for (int k = 0; k < m; k++) if (!adj[p[eu[k]]][p[ev[k]]]) return 0;
    return 1;
}
static void gen_perms(void) {           /* Heap's algorithm over all n! permutations */
    unsigned char a[MAXN], c[MAXN]; for (int i = 0; i < n; i++) { a[i] = (unsigned char)i; c[i] = 0; }
    naut = 0;
    if (is_aut(a)) memcpy(perm[naut++], a, n);
    int i = 0;
    while (i < n) {
        if (c[i] < i) {
            if (i % 2 == 0) { unsigned char t = a[0]; a[0] = a[i]; a[i] = t; } else { unsigned char t = a[c[i]]; a[c[i]] = a[i]; a[i] = t; }
            if (is_aut(a)) { if (naut >= MAXAUT) { fprintf(stderr, "orbitenum: |Aut(G)| > %d, refusing\n", MAXAUT); exit(3); } memcpy(perm[naut++], a, n); }
            c[i]++; i = 0;
        } else { c[i] = 0; i++; }
    }
}
static inline u64 img(int g, u64 mask) {
    u64 r = 0; const u64 *t = tbl + (size_t)g * nbytes * 256;
    for (int b = 0; b < nbytes; b++) r |= t[b * 256 + ((mask >> (8 * b)) & 255)];
    return r;
}
static char outbuf[2 + (MAXN * MAXN + 5) / 6 + 1];
static void emit(u64 Bmask, u64 O, const int *F, int f) {
    unsigned char D[MAXN][MAXN]; memset(D, 0, sizeof D);
    for (int k = 0; k < m; k++) if (Bmask >> k & 1) D[eu[k]][ev[k]] = D[ev[k]][eu[k]] = 1;
    for (int t = 0; t < f; t++) { int k = F[t]; if (O >> t & 1) D[ev[k]][eu[k]] = 1; else D[eu[k]][ev[k]] = 1; }
    int pos = 0, cur = 0, nb = 0; outbuf[pos++] = '&'; outbuf[pos++] = (char)(63 + n);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { cur = (cur << 1) | D[i][j]; if (++nb == 6) { outbuf[pos++] = (char)(63 + cur); cur = 0; nb = 0; } }
    if (nb) { cur <<= (6 - nb); outbuf[pos++] = (char)(63 + cur); }
    outbuf[pos++] = '\n'; fwrite(outbuf, 1, pos, stdout);
}
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: orbitenum <graph6 G> <min_arcs> [p/P] [-u]\n"); return 1; }
    int min_arcs = atoi(argv[2]); int p = 0, P = 1, count_only = 0;
    for (int a = 3; a < argc; a++) { if (!strcmp(argv[a], "-u")) count_only = 1; else if (sscanf(argv[a], "%d/%d", &p, &P) != 2 || P < 1 || p < 0 || p >= P) { fprintf(stderr, "bad argument %s\n", argv[a]); return 1; } }
    char gbuf[256]; const char *gstr = argv[1];
    if (gstr[0] == '@') {                   /* @file: read the graph6 string from the first line of the file */
        FILE *fh = fopen(gstr + 1, "r"); if (!fh || !fgets(gbuf, sizeof gbuf, fh)) { fprintf(stderr, "cannot read %s\n", gstr + 1); return 1; }
        fclose(fh); gbuf[strcspn(gbuf, "\r\n")] = 0; gstr = gbuf; }
    n = decode_g6(gstr); if (n < 0) { fprintf(stderr, "bad graph6\n"); return 1; }
    m = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) pidx[i][j] = -1;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (adj[i][j]) { eu[m] = i; ev[m] = j; pidx[i][j] = pidx[j][i] = m; m++; }
    if (m > MAXE || m > 62) { fprintf(stderr, "too many edges\n"); return 1; }
    gen_perms();
    for (int g = 0; g < naut; g++) for (int k = 0; k < m; k++) {
        int a = perm[g][eu[k]], b = perm[g][ev[k]]; emap[g][k] = (unsigned char)pidx[a][b]; eflip[g][k] = (unsigned char)(a > b); }
    nbytes = (m + 7) / 8; tbl = calloc((size_t)naut * nbytes * 256, sizeof(u64)); if (!tbl) { fprintf(stderr, "oom\n"); return 1; }
    for (int g = 0; g < naut; g++) for (int b = 0; b < nbytes; b++) for (int v = 0; v < 256; v++) {
        u64 r = 0; for (int k = 0; k < 8; k++) if ((v >> k & 1) && 8 * b + k < m) r |= (u64)1 << emap[g][8 * b + k];
        tbl[((size_t)g * nbytes + b) * 256 + v] = r; }
    int lo = min_arcs - m; if (lo < 0) lo = 0;
    u64 total = 0; long long nsub = 0, idx = 0;
    static int stab[MAXAUT]; static unsigned char posmap[MAXAUT][MAXE], pflip[MAXAUT][MAXE]; int F[MAXE], posOf[MAXE];
    for (u64 B = 0; B < ((u64)1 << m); B++) {
        if (__builtin_popcountll(B) < lo) continue;
        int canon = 1; for (int g = 1; g < naut; g++) if (img(g, B) < B) { canon = 0; break; }
        if (!canon) continue;
        long long my = idx++; if (my % P != p) continue;
        nsub++;
        int ns = 0; for (int g = 1; g < naut; g++) if (img(g, B) == B) stab[ns++] = g;
        int f = 0; for (int k = 0; k < m; k++) if (!(B >> k & 1)) { posOf[k] = f; F[f++] = k; }
        for (int s = 0; s < ns; s++) for (int t = 0; t < f; t++) { int g = stab[s]; posmap[s][t] = (unsigned char)posOf[emap[g][F[t]]]; pflip[s][t] = eflip[g][F[t]]; }
        u64 cnt = 0;
        for (u64 O = 0; O < ((u64)1 << f); O++) {
            int keep = 1;
            for (int s = 0; s < ns && keep; s++) {
                u64 io = 0; for (int t = 0; t < f; t++) io |= (u64)(((O >> t) & 1) ^ pflip[s][t]) << posmap[s][t];
                if (io < O) keep = 0;
            }
            if (!keep) continue;
            cnt++; if (!count_only) emit(B, O, F, f);
        }
        total += cnt; fprintf(stderr, "B %llu %llu\n", (unsigned long long)B, (unsigned long long)cnt);
    }
    fprintf(stderr, "orbitenum %s min_arcs %d part %d/%d aut %d subsets %lld total %llu\n", gstr, min_arcs, p, P, naut, nsub, (unsigned long long)total);
    return 0;
}
