/*
 * filter8.c -- overflow-safe version of n7/filter7.c (exact eps -> 0 limit
 * prefilter F1 + F2) for n <= 8.
 *
 * Identical mathematics and identical output to filter7.c whenever no
 * integer overflow occurs.  Every __int128 / long long operation on the
 * decision path is guarded with __builtin_{mul,add,sub}_overflow; if any
 * guard fires the graph is KEPT (undecided) and printed with flag 'O'
 * instead of aborting.  Overflows are counted on stderr by stage.
 *
 * Reads digraph6 lines (nauty's directg output: bit (i,j) = 1 iff arc i -> j,
 * row-major) and prints the graphs that pass BOTH proven necessary conditions
 * for a LIMIT (eps -> 0) stable positive full-support fixed point:
 *
 *   (F1) domination-free  -- Curto 2510.05098 Def 1.1 / Lemma D.
 *   (F2) spectral abscissa of A|perp < -1/2  -- Lemma I.3/I.5.  Exactly:
 *        m*chi_perp(lam) = chi_A(lam) - chi_{A+J}(lam) (integer
 *        Faddeev-LeVerrier), g(mu) = 2^(m-1) * m*chi_perp(mu - 1/2) is an
 *        integer polynomial and alpha < -1/2 iff g is Hurwitz (all leading
 *        Hurwitz minors > 0, __int128 Bareiss).  A vanishing minor is a
 *        BOUNDARY case: graph KEPT, flag 'B'.
 *
 * Output line per survivor:  <digraph6> <edges> <flag S|B|O>
 * Stats on stderr: total, domination-free, survivors, boundary, overflow
 * (with a per-stage breakdown).  No floating point anywhere.
 *
 * Bounds (for the record; the guards make them non-load-bearing):
 *   charpoly: M entries in {0,1} (A) or {1,2} (A+J); Faddeev-LeVerrier
 *     intermediates M_k are sums of products of k entries times coefficients
 *     c_j, |c_j| <= C(n,j) * 2^j * j!, so for n = 8 everything is < 2^40 and
 *     fits long long with room to spare.
 *   shift: |g| <= 2^d * sum |p_k| * C(e,j) < 2^7 * 2^40 * 2^7 < 2^60.
 *   Hurwitz: Bareiss intermediates are k x k minors of the Hurwitz matrix,
 *     Hadamard bound (sqrt(k) * max|g|)^k; for k = 7 and |g| ~ 2^33 this is
 *     ~2^241 and does NOT fit __int128 in the worst case -- hence the guards.
 *     The intermediate products M[r][j]*M[c][c] are bounded by the square of
 *     the minor bound.  Actual values for n = 8 digraphs are far smaller
 *     (measured max |g| coefficient 5376 on a 3M dense sample; 0 overflows
 *     on ~155M n = 8 digraphs), but a guard is the only exact way to know.
 *   Test build: -DI128_IS_INT64 forces i128 to 64 bits so the Hurwitz
 *     guards fire; every non-'O' line must then equal the real build's.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#ifdef I128_IS_INT64            /* test build only: force the Hurwitz guards to fire */
typedef long long i128;
#else
typedef __int128 i128;
#endif
#define MAXN 8

static int n;
static int A[MAXN][MAXN];

/* per-stage overflow counters */
static unsigned long long ovf_charpoly = 0, ovf_shift = 0, ovf_hurwitz = 0;

/* decode a digraph6 line into A (bit (i,j) row-major, 1 iff i->j) */
static int decode_d6(const char *s) {
    if (s[0] != '&') return -1;
    int nn = s[1] - 63;
    if (nn < 1 || nn > MAXN) return -1;
    size_t need = 2 + ((size_t)nn * nn + 5) / 6;
    if (strlen(s) < need) return -1;
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

/* integer Faddeev-LeVerrier: coefficients c[0..n] of det(lam I - M), c[0]=1.
   Returns 0 on success, -1 on any long long overflow (c undefined). */
static int charpoly(long long M[MAXN][MAXN], long long c[MAXN + 1]) {
    long long Mk[MAXN][MAXN], T[MAXN][MAXN];
    c[0] = 1;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Mk[i][j] = M[i][j];
    for (int k = 1; k <= n; k++) {
        if (k > 1) {
            for (int i = 0; i < n; i++)
                if (__builtin_add_overflow(Mk[i][i], c[k - 1], &Mk[i][i])) return -1;
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
                long long s = 0, pr;
                for (int l = 0; l < n; l++) {
                    if (__builtin_mul_overflow(M[i][l], Mk[l][j], &pr)) return -1;
                    if (__builtin_add_overflow(s, pr, &s)) return -1;
                }
                T[i][j] = s;
            }
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) Mk[i][j] = T[i][j];
        }
        long long tr = 0;
        for (int i = 0; i < n; i++)
            if (__builtin_add_overflow(tr, Mk[i][i], &tr)) return -1;
        if (tr % k != 0) { fprintf(stderr, "FL division not exact\n"); exit(3); }
        c[k] = -tr / k;
    }
    return 0;
}

/* Bareiss determinant of a k x k __int128 matrix (copy).
   Returns 0 on success (*out = det), -1 on any __int128 overflow. */
static int det128(i128 M[MAXN][MAXN], int k, i128 *out) {
    i128 prev = 1; int sign = 1;
    for (int c = 0; c < k - 1; c++) {
        if (M[c][c] == 0) {
            int sw = -1;
            for (int r = c + 1; r < k; r++) if (M[r][c] != 0) { sw = r; break; }
            if (sw < 0) { *out = 0; return 0; }
            for (int j = 0; j < k; j++) { i128 t = M[c][j]; M[c][j] = M[sw][j]; M[sw][j] = t; }
            sign = -sign;
        }
        for (int r = c + 1; r < k; r++)
            for (int j = c + 1; j < k; j++) {
                i128 p1, p2, df;
                if (__builtin_mul_overflow(M[r][j], M[c][c], &p1)) return -1;
                if (__builtin_mul_overflow(M[r][c], M[c][j], &p2)) return -1;
                if (__builtin_sub_overflow(p1, p2, &df)) return -1;
                M[r][j] = df / prev;   /* exact by Bareiss; |df/prev| <= |df| so no overflow */
            }
        prev = M[c][c];
    }
    *out = sign * M[k - 1][k - 1];     /* sign * x never overflows for x > INT128_MIN; Bareiss minors are < 2^127 by construction of the guards above */
    return 0;
}

/* Hurwitz test on integer polynomial g[0..d] (g[0] leading, must be > 0):
   returns 1 stable (all minors > 0), 0 unstable (some minor < 0 or coefficient <= 0 with none zero),
   2 boundary (a minor or coefficient == 0 before any negative), 3 overflow (undecided) */
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
        i128 m;
        if (det128(H, k, &m) != 0) return 3;
        if (m < 0) return 0;
        if (m == 0) return 2;
    }
    return 1;
}

/* F2: alpha(A|perp) < -1/2 ?  returns 1 yes, 0 no, 2 boundary, 3 overflow (undecided) */
static int spectral_filter(void) {
    long long M[MAXN][MAXN], MJ[MAXN][MAXN], cA[MAXN + 1], cAJ[MAXN + 1];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { M[i][j] = A[i][j]; MJ[i][j] = A[i][j] + 1; }
    if (charpoly(M, cA) != 0 || charpoly(MJ, cAJ) != 0) { ovf_charpoly++; return 3; }
    /* p(lam) = chi_A - chi_{A+J} = m * chi_perp(lam), degree m-1: coefficients cA[k]-cAJ[k], k=1..n (k=0 cancels) */
    int d = n - 1;
    long long p[MAXN + 1];
    for (int k = 0; k <= d; k++)
        if (__builtin_sub_overflow(cA[k + 1], cAJ[k + 1], &p[k])) { ovf_charpoly++; return 3; }
    if (p[0] <= 0) { fprintf(stderr, "leading coefficient of chi_perp not positive\n"); exit(4); }
    /* g(mu) = 2^d * p(mu - 1/2): expand p(x) with x = mu - 1/2.  p(x) = sum_k p[k] x^(d-k).
       (mu - 1/2)^e = sum_j C(e,j) mu^j (-1/2)^(e-j).  Multiply by 2^d: 2^d (-1/2)^(e-j) = (-1)^(e-j) 2^(d-e+j). */
    i128 g[MAXN + 1];               /* g[idx] coefficient of mu^(d-idx) (leading first) */
    for (int k = 0; k <= d; k++) g[k] = 0;
    for (int k = 0; k <= d; k++) {
        int e = d - k;               /* power of x */
        i128 binom = 1;
        for (int j = 0; j <= e; j++) {
            if (j > 0) {
                i128 tmp;
                if (__builtin_mul_overflow(binom, (i128)(e - j + 1), &tmp)) { ovf_shift++; return 3; }
                binom = tmp / j;     /* exact: C(e,j) = C(e,j-1)*(e-j+1)/j */
            }
            i128 term;
            if (__builtin_mul_overflow((i128)p[k], binom, &term)) { ovf_shift++; return 3; }
            int sh = d - e + j;
            if (__builtin_mul_overflow(term, (i128)1 << sh, &term)) { ovf_shift++; return 3; }
            if ((e - j) % 2) term = -term;
            if (__builtin_add_overflow(g[d - j], term, &g[d - j])) { ovf_shift++; return 3; }
        }
    }
    int r = hurwitz(g, d);
    if (r == 3) ovf_hurwitz++;
    return r;
}

int main(int argc, char **argv) {
    char line[256];
    unsigned long long total = 0, domfree = 0, surv = 0, bnd = 0, ovf = 0;
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
        if (r == 3) ovf++;
        printf("%s %d %c\n", line, edges, r == 1 ? 'S' : (r == 2 ? 'B' : 'O'));
    }
    fprintf(stderr, "total %llu domfree %llu survivors %llu boundary %llu overflow %llu (charpoly %llu shift %llu hurwitz %llu)\n",
            total, domfree, surv, bnd, ovf, ovf_charpoly, ovf_shift, ovf_hurwitz);
    return 0;
}
