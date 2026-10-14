/*
 * filter7.c -- exact C filter for the n-node full-support sweep (n <= 8).
 *
 * Reads digraph6 lines (nauty's directg output: bit (i,j) = 1 iff arc i -> j,
 * row-major; the convention is verified against showg in the build script)
 * and prints the graphs that pass BOTH proven necessary conditions for a
 * LIMIT (eps -> 0) stable positive full-support fixed point:
 *
 *   (F1) domination-free  -- Curto 2510.05098 Def 1.1 / our Lemma D: a
 *        dominated node forbids any full-support fixed point at any parameter.
 *   (F2) spectral abscissa of A|perp (adjacency compressed to the hyperplane
 *        orthogonal to the all-ones vector) is < -1/2  -- Lemma I.3/I.5: the
 *        limit-stable beta-set is (2, 1/(1+alpha)) and is nonempty iff
 *        alpha < -1/2.  Computed EXACTLY: m*chi_perp(lam) = chi_A(lam) -
 *        chi_{A+J}(lam) (integer Faddeev-LeVerrier), then g(mu) = 2^(m-1) *
 *        m*chi_perp(mu - 1/2) is an integer polynomial and alpha < -1/2 iff g
 *        is Hurwitz (all leading Hurwitz minors > 0, __int128 Bareiss).  A
 *        vanishing minor (alpha == -1/2 possible) is a BOUNDARY case and the
 *        graph is KEPT (conservative), flagged with 'B'.
 *
 * Output line per survivor:  <digraph6> <edges> <flag S|B>
 * Stats on stderr: total, domination-free, survivors, boundary.
 * No floating point anywhere.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef __int128 i128;
#define MAXN 8

static int n;
static int A[MAXN][MAXN];

/* decode a digraph6 line into A (bit (i,j) row-major, 1 iff i->j) */
static int decode_d6(const char *s) {
    if (s[0] != '&') return -1;
    int nn = s[1] - 63;
    if (nn < 1 || nn > MAXN) return -1;
    n = nn;
    const char *p = s + 2;
    int bit = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        int c = p[bit / 6] - 63;
        int b = (c >> (5 - bit % 6)) & 1;
        A[i][j] = b;
        bit++;
    }
    return 0;
}

/* F1: domination-free.  k dominates j iff j->k, not k->j, and every i not in {j,k} with i->j has i->k */
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

/* integer Faddeev-LeVerrier: coefficients c[0..n] of det(lam I - M), c[0]=1 */
static void charpoly(long long M[MAXN][MAXN], long long c[MAXN + 1]) {
    long long Mk[MAXN][MAXN], T[MAXN][MAXN];
    c[0] = 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Mk[i][j] = M[i][j];
    for (int k = 1; k <= n; k++) {
        if (k > 1) {
            for (int i = 0; i < n; i++) Mk[i][i] += c[k - 1];
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
                long long s = 0;
                for (int l = 0; l < n; l++) s += M[i][l] * Mk[l][j];
                T[i][j] = s;
            }
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Mk[i][j] = T[i][j];
        }
        long long tr = 0;
        for (int i = 0; i < n; i++) tr += Mk[i][i];
        if (tr % k != 0) { fprintf(stderr, "FL division not exact\n"); exit(3); }
        c[k] = -tr / k;
    }
}

/* Bareiss determinant of a k x k __int128 matrix (copy) */
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

/* Hurwitz test on integer polynomial g[0..d] (g[0] leading, must be > 0):
   returns 1 stable (all minors > 0), 0 unstable (some minor < 0 or coefficient <= 0 with none zero), 2 boundary (a minor or coefficient == 0 before any negative) */
static int hurwitz(i128 g[], int d) {
    if (g[0] <= 0) return 0;
    int zero = 0;
    for (int k = 1; k <= d; k++) { if (g[k] < 0) return 0; if (g[k] == 0) zero = 1; }
    if (zero) return 2;
    for (int k = 1; k <= d; k++) {
        i128 H[MAXN][MAXN];
        for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) {
            int idx = 2 * (i + 1) - (j + 1);
            H[i][j] = (idx >= 0 && idx <= d) ? g[idx] : 0;
        }
        i128 m = det128(H, k);
        if (m < 0) return 0;
        if (m == 0) return 2;
    }
    return 1;
}

/* F2: alpha(A|perp) < -1/2 ?  returns 1 yes, 0 no, 2 boundary */
static int spectral_filter(void) {
    long long M[MAXN][MAXN], MJ[MAXN][MAXN], cA[MAXN + 1], cAJ[MAXN + 1];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { M[i][j] = A[i][j]; MJ[i][j] = A[i][j] + 1; }
    charpoly(M, cA); charpoly(MJ, cAJ);
    /* p(lam) = chi_A - chi_{A+J} = m * chi_perp(lam), degree m-1: coefficients cA[k]-cAJ[k], k=1..n (k=0 cancels) */
    int d = n - 1;
    long long p[MAXN + 1];
    for (int k = 0; k <= d; k++) p[k] = cA[k + 1] - cAJ[k + 1];
    if (p[0] <= 0) { fprintf(stderr, "leading coefficient of chi_perp not positive\n"); exit(4); }
    /* g(mu) = 2^d * p(mu - 1/2): expand p(x) with x = mu - 1/2.  p(x) = sum_k p[k] x^(d-k).
       (mu - 1/2)^e = sum_j C(e,j) mu^j (-1/2)^(e-j).  Multiply by 2^d: 2^d (-1/2)^(e-j) = (-1)^(e-j) 2^(d-e+j). */
    i128 g[MAXN + 1];               /* g[idx] coefficient of mu^(d-idx) (leading first) */
    for (int k = 0; k <= d; k++) g[k] = 0;
    for (int k = 0; k <= d; k++) {
        int e = d - k;               /* power of x */
        i128 binom = 1;
        for (int j = 0; j <= e; j++) {
            if (j > 0) binom = binom * (e - j + 1) / j;
            i128 term = (i128)p[k] * binom;
            int sh = d - e + j;
            term = term * ((i128)1 << sh);
            if ((e - j) % 2) term = -term;
            g[d - j] += term;        /* mu^j -> index d-j */
        }
    }
    return hurwitz(g, d);
}

int main(int argc, char **argv) {
    char line[256];
    unsigned long long total = 0, domfree = 0, surv = 0, bnd = 0;
    int print_all_domfree = (argc > 1 && strcmp(argv[1], "--domfree") == 0);
    while (fgets(line, sizeof line, stdin)) {
        size_t L = strlen(line);
        while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
        if (!L) continue;
        if (decode_d6(line) != 0) { fprintf(stderr, "bad line: %s\n", line); exit(2); }
        total++;
        if (!domination_free()) continue;
        domfree++;
        int edges = 0;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) edges += A[i][j];
        if (print_all_domfree) { printf("%s %d D\n", line, edges); continue; }
        int r = spectral_filter();
        if (r == 0) continue;
        surv++;
        if (r == 2) bnd++;
        printf("%s %d %c\n", line, edges, r == 1 ? 'S' : 'B');
    }
    fprintf(stderr, "total %llu domfree %llu survivors %llu boundary %llu\n", total, domfree, surv, bnd);
    return 0;
}
