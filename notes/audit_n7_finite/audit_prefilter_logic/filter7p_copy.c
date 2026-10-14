/*
 * filter7p.c -- exact C prefilter for the FINITE-eps full-support map (n <= 8).
 *
 * Keeps every digraph that can possibly carry a stable positive full-support
 * fixed point at SOME legal (eps, delta), by two proven necessary conditions
 * that hold at every parameter value (not only in the eps -> 0 limit):
 *
 *   (F1) domination-free -- Curto 2510.05098 Def 1.1 / Lemma D.
 *   (F3) positivity somewhere -- Lemma A.2: a stable full-support fixed point
 *        has det K > 0 and x = adj(K) 1 / det K, and adj(K) 1 = a(beta) =
 *        adj(I + beta N) 1 is t-free (Lemma A), so stability at legal (beta, t)
 *        forces a_i(beta) > 0 for all i at some beta > 2 (legal <=> beta > 2 + 1/t).
 *        Here N[i][j] = 1 iff j -/-> i (i != j), the receiver-first nonedge
 *        matrix, exactly as fastregion.beta_data.
 *
 * F3 is decided in C only when the decision is exact and one-sided:
 *   drop 'N'  if for some i every monomial coefficient of a_i(2+u) is <= 0
 *             (then a_i <= 0 for all u > 0, i.e. for all beta > 2);
 *   keep 'P'  if at some rational beta > 2 every a_i(beta) > 0 (exact integer
 *             determinant evaluation; this certifies positivity somewhere);
 *   keep 'U'  otherwise -- undecided here; Python (fastregion.region with
 *             positivity_only=True, exact root isolation) decides.
 * The coefficients of a_i(2+u) are obtained by exact evaluation at
 * u = 0..m (integer determinants, Bareiss in __int128), Newton forward
 * differences and a check at u = m+1 (degree bound verified).
 *
 * Column 4 is informational: 1 iff the graph has two non-adjacent in-twins
 * (Remark J.4': such a support is never stable); it is NOT used to drop.
 *
 * Output line per kept graph:  <digraph6> <edges> <P|U> <twin 0|1>
 * Stats on stderr. No floating point anywhere.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef __int128 i128;
#define MAXN 8

static int n;
static int A[MAXN][MAXN];
static int N[MAXN][MAXN];

static int decode_d6(const char *s) {
    if (s[0] != '&') return -1;
    int nn = s[1] - 63;
    if (nn < 1 || nn > MAXN) return -1;
    n = nn;
    const char *p = s + 2;
    int bit = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        int c = p[bit / 6] - 63;
        A[i][j] = (c >> (5 - bit % 6)) & 1;
        bit++;
    }
    return 0;
}

/* F1: k dominates j iff j->k, not k->j, and every i not in {j,k} with i->j has i->k */
static int domination_free(void) {
    for (int j = 0; j < n; j++) for (int k = 0; k < n; k++) {
        if (j == k) continue;
        if (!(A[j][k] == 1 && A[k][j] == 0)) continue;
        int dom = 1;
        for (int i = 0; i < n; i++) {
            if (i == j || i == k) continue;
            if (A[i][j] == 1 && A[i][k] == 0) { dom = 0; break; }
        }
        if (dom) return 0;
    }
    return 1;
}

/* Remark J.4': u != v non-adjacent with identical in-neighbourhoods among the other nodes */
static int nonadjacent_intwins(void) {
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) {
        if (A[u][v] || A[v][u]) continue;
        int tw = 1;
        for (int k = 0; k < n; k++) {
            if (k == u || k == v) continue;
            if (A[k][u] != A[k][v]) { tw = 0; break; }
        }
        if (tw) return 1;
    }
    return 0;
}

static i128 det128(i128 M[MAXN][MAXN], int k) {
    i128 prev = 1; int sign = 1;
    for (int c = 0; c < k - 1; c++) {
        if (M[c][c] == 0) {
            int sw = -1;
            for (int r = c + 1; r < k; r++) if (M[r][c] != 0) { sw = r; break; }
            if (sw < 0) return 0;
            for (int j = 0; j < k; j++) { i128 t = M[c][j]; M[c][j] = M[sw][j]; M[sw][j] = t; }
            sign = -sign;
        }
        for (int r = c + 1; r < k; r++)
            for (int j = c + 1; j < k; j++)
                M[r][j] = (M[r][j] * M[c][c] - M[r][c] * M[c][j]) / prev;
        prev = M[c][c];
    }
    return sign * M[k - 1][k - 1];
}

/* q^m * a_i(p/q): determinant of (q I + p N) with column i replaced by q * ones (Cramer) */
static i128 a_eval(int i, i128 p, i128 q) {
    i128 M[MAXN][MAXN];
    for (int r = 0; r < n; r++) for (int c = 0; c < n; c++)
        M[r][c] = (r == c ? q : 0) + p * N[r][c];
    for (int r = 0; r < n; r++) M[r][i] = q;
    return det128(M, n);
}

/* monomial coefficients of m! * a_i(2+u), degree <= m, from values at u = 0..m; checked at u = m+1.
   returns 1 if all coefficients <= 0 (so a_i <= 0 on u > 0), else 0 */
static i128 g_coef[MAXN + 2];
static int descartes_nonpositive(int i) {
    int m = n;
    i128 f[MAXN + 2];
    for (int k = 0; k <= m + 1; k++) f[k] = a_eval(i, 2 + k, 1);
    /* forward differences D[k] = Delta^k f(0) */
    i128 D[MAXN + 2], w[MAXN + 2];
    for (int k = 0; k <= m; k++) w[k] = f[k];
    for (int k = 0; k <= m; k++) {
        D[k] = w[0];
        for (int j = 0; j < m - k; j++) w[j] = w[j + 1] - w[j];
    }
    /* m! * p(u) = sum_k D[k] * (m!/k!) * FF_k(u), FF_k = u(u-1)...(u-k+1) */
    i128 fact = 1;
    for (int k = 2; k <= m; k++) fact *= k;
    i128 coef[MAXN + 2]; for (int j = 0; j <= m; j++) coef[j] = 0;
    i128 FF[MAXN + 2]; for (int j = 0; j <= m; j++) FF[j] = 0; FF[0] = 1;   /* FF_0 = 1 */
    i128 kfact = 1;
    for (int k = 0; k <= m; k++) {
        if (k > 0) {
            kfact *= k;
            /* FF_k = FF_{k-1} * (u - (k-1)) */
            for (int j = k; j >= 1; j--) FF[j] = FF[j - 1] - (i128)(k - 1) * FF[j];
            FF[0] = -(i128)(k - 1) * FF[0];
        }
        i128 scale = D[k] * (fact / kfact);
        for (int j = 0; j <= k; j++) coef[j] += scale * FF[j];
    }
    for (int jj = 0; jj <= m; jj++) g_coef[jj] = coef[jj];
    /* degree-bound check at u = m+1 */
    i128 val = 0, x = m + 1;
    for (int j = m; j >= 0; j--) val = val * x + coef[j];
    if (val != fact * f[m + 1]) { fprintf(stderr, "interpolation check failed\n"); exit(5); }
    for (int j = 0; j <= m; j++) if (coef[j] > 0) return 0;
    return 1;
}

static const int SAMPLES[][2] = { {1,1},{1,2},{1,4},{1,8},{2,1},{1,16},{4,1},{1,32},{8,1},{16,1},{48,1},{1,64},{3,8},{3,4} };
#define NSAMPLES (int)(sizeof(SAMPLES)/sizeof(SAMPLES[0]))

int main(int argc, char **argv) {
    char line[256];
    unsigned long long total = 0, domfree = 0, dropN = 0, keptP = 0, keptU = 0, twins = 0;
    int statsonly = (argc > 1 && strcmp(argv[1], "--stats-only") == 0);
    while (fgets(line, sizeof line, stdin)) {
        size_t L = strlen(line);
        while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
        if (!L) continue;
        char *sp = strchr(line, ' '); if (sp) *sp = 0;   /* accept "<d6> ..." lines too */
        if (decode_d6(line) != 0) { fprintf(stderr, "bad line: %s\n", line); exit(2); }
        total++;
        if (!domination_free()) continue;
        domfree++;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            N[i][j] = (i != j && A[j][i] == 0) ? 1 : 0;      /* N[i][j] = 1 iff j -/-> i */
        int drop = 0;
        for (int i = 0; i < n && !drop; i++) if (descartes_nonpositive(i)) drop = 1;
        if (drop) { dropN++; continue; }
        char flag = 'U';
        for (int s = 0; s < NSAMPLES && flag == 'U'; s++) {
            i128 p = 2 * SAMPLES[s][1] + SAMPLES[s][0], q = SAMPLES[s][1];   /* beta = 2 + p'/q */
            int allpos = 1;
            for (int i = 0; i < n; i++) if (a_eval(i, p, q) <= 0) { allpos = 0; break; }
            if (allpos) flag = 'P';
        }
        if (flag == 'P') keptP++; else keptU++;
        int tw = nonadjacent_intwins(); twins += tw;
        int edges = 0;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) edges += A[i][j];
        if (!statsonly) printf("%s %d %c %d\n", line, edges, flag, tw);
    }
    fprintf(stderr, "total %llu domfree %llu dropN %llu keptP %llu keptU %llu twinflag %llu\n",
            total, domfree, dropN, keptP, keptU, twins);
    return 0;
}
