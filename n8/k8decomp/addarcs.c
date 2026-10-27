/*
 * addarcs.c -- OR a fixed set of bidirected arcs into every digraph of a digraph6 stream.
 *
 *   directg -q -o H | addarcs B        (H = complement of B on the same labelled vertex set)
 *
 * Reads digraph6 lines from stdin (nauty: '&', N(n), then the n*n adjacency bits row-major, bit (i,j)
 * = 1 iff arc i -> j), sets both arcs i -> j and j -> i for every edge {i,j} of the graph6 graph B given
 * as argv[1], and writes the result as digraph6.  Lines starting with '>' (format headers) are skipped;
 * a vertex-count mismatch or a malformed line is fatal (exit 2).  Prints "addarcs <B> lines <count>"
 * on stderr at the end so that the per-B counts can be audited.  No other transformation is made.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXN 60

static int nB;
static unsigned char B[MAXN][MAXN];
static unsigned char D[MAXN][MAXN];

static int decode_g6(const char *s, unsigned char adj[MAXN][MAXN]) {
    int n = s[0] - 63;
    if (n < 1 || n > MAXN) return -1;
    memset(adj, 0, (size_t)MAXN * MAXN);
    const char *p = s + 1; int k = 0;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++) {
            int c = (unsigned char)p[k / 6];
            if (c < 63 || c > 126) return -1;
            int bit = ((c - 63) >> (5 - k % 6)) & 1;
            adj[i][j] = adj[j][i] = (unsigned char)bit; k++;
        }
    return n;
}

static int decode_d6(const char *s, unsigned char adj[MAXN][MAXN]) {
    if (s[0] != '&') return -1;
    int n = s[1] - 63;
    if (n < 1 || n > MAXN) return -1;
    memset(adj, 0, (size_t)MAXN * MAXN);
    const char *p = s + 2; int k = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int c = (unsigned char)p[k / 6];
            if (c < 63 || c > 126) return -1;
            adj[i][j] = (unsigned char)(((c - 63) >> (5 - k % 6)) & 1); k++;
        }
    return n;
}

static void encode_d6(int n, unsigned char adj[MAXN][MAXN], char *out) {
    int pos = 0, cur = 0, nb = 0;
    out[pos++] = '&'; out[pos++] = (char)(63 + n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            cur = (cur << 1) | (adj[i][j] ? 1 : 0); nb++;
            if (nb == 6) { out[pos++] = (char)(63 + cur); cur = 0; nb = 0; }
        }
    if (nb > 0) { cur <<= (6 - nb); out[pos++] = (char)(63 + cur); }
    out[pos] = 0;
}

int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: addarcs <graph6 B>  < digraph6 stream\n"); return 1; }
    nB = decode_g6(argv[1], B);
    if (nB < 0) { fprintf(stderr, "addarcs: bad graph6 B %s\n", argv[1]); return 2; }
    char *line = NULL; size_t cap = 0; ssize_t len;
    char out[2 + (MAXN * MAXN) / 6 + 2];
    unsigned long long count = 0;
    while ((len = getline(&line, &cap, stdin)) != -1) {
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) line[--len] = 0;
        if (len == 0 || line[0] == '>') continue;
        int n = decode_d6(line, D);
        if (n < 0) { fprintf(stderr, "addarcs: bad digraph6 line %s\n", line); return 2; }
        if (n != nB) { fprintf(stderr, "addarcs: vertex count %d != %d of B\n", n, nB); return 2; }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (B[i][j]) D[i][j] = 1;
        encode_d6(n, D, out);
        puts(out); count++;
    }
    free(line);
    fprintf(stderr, "addarcs %s lines %llu\n", argv[1], count);
    return 0;
}
