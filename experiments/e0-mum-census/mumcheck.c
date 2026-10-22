// mumcheck.c -- pairwise MUM anchoring experiment for VNTR alleles.
// stdin, one pair per line:  id \t Lmin \t ca0 \t ca1 \t cb0 \t cb1 \t doDP \t A \t B
//   [ca0,ca1) is the "core" (array) zone of A, [cb0,cb1) of B (used only for MUMs crossing zones).
// stdout, one line per pair (tab-separated, see HEADER below).
// MUM = maximal exact match, length >= Lmin, occurring exactly once in A and once in B (MUMmer's
// definition, via a suffix array of A#B$ and its LCP array), left-maximal.
// With doDP: opt = global unit edit distance; for every top-level MUM, excess = ED(A[:p],B[:q]) +
// ED(A[p+l:],B[q+l:]) - opt (0 = the MUM lies on some optimal alignment); chain = heaviest collinear
// non-overlapping chain of top-level MUMs, chain_cost = sum of ED of the gaps; recursive partition =
// chain, then recompute MUMs (unique within the gap substrings) inside every gap, recursively;
// rec_cost = sum of ED of the leaves; rec_maxwin = largest leaf side; rec_cells = sum of leaf n*m.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { int pa, pb, len; } Mum;

static int *g_rank, g_k, g_n;
static int cmp_sa(const void *x, const void *y) {
    int i = *(const int *)x, j = *(const int *)y;
    if (g_rank[i] != g_rank[j]) return g_rank[i] < g_rank[j] ? -1 : 1;
    int ri = i + g_k < g_n ? g_rank[i + g_k] : -1;
    int rj = j + g_k < g_n ? g_rank[j + g_k] : -1;
    return ri < rj ? -1 : ri > rj;
}

static int code(char c) {
    switch (c) { case 'A': case 'a': return 2; case 'C': case 'c': return 3;
                 case 'G': case 'g': return 4; case 'T': case 't': return 5; default: return 6; }
}

// MUMs between A[a0,a1) and B[b0,b1), returned with absolute coordinates
static Mum *find_mums(const char *A, int a0, int a1, const char *B, int b0, int b1, int Lmin,
                      int *nout, int *lcs_out) {
    int na = a1 - a0, nb = b1 - b0, n = na + nb + 2;
    *nout = 0; *lcs_out = 0;
    if (na < Lmin || nb < Lmin) return NULL;
    int *s = malloc(sizeof(int) * n), *sa = malloc(sizeof(int) * n), *rk = malloc(sizeof(int) * n),
        *tmp = malloc(sizeof(int) * n), *lcp = malloc(sizeof(int) * (n + 1));
    for (int i = 0; i < na; i++) s[i] = code(A[a0 + i]);
    s[na] = 0;
    for (int i = 0; i < nb; i++) s[na + 1 + i] = code(B[b0 + i]);
    s[n - 1] = 1;
    // N never matches: give every N a unique rank by treating it as distinct below (simple: code 6+i)
    for (int i = 0; i < n; i++) { sa[i] = i; rk[i] = s[i] == 6 ? 7 + i : s[i]; }
    g_n = n; g_rank = rk;
    for (int k = 1;; k <<= 1) {
        g_k = k;
        qsort(sa, n, sizeof(int), cmp_sa);
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; i++) tmp[sa[i]] = tmp[sa[i - 1]] + (cmp_sa(&sa[i - 1], &sa[i]) < 0);
        memcpy(rk, tmp, sizeof(int) * n);
        if (rk[sa[n - 1]] == n - 1) break;
        if (k > n) break;
    }
    // Kasai: lcp[i] = lcp(sa[i-1], sa[i]); lcp[0] = lcp[n] = 0
    int h = 0;
    for (int i = 0; i < n; i++) {
        if (rk[i] > 0) {
            int j = sa[rk[i] - 1];
            while (i + h < n && j + h < n && s[i + h] == s[j + h] && s[i + h] != 6 && s[i + h] > 1) h++;
            lcp[rk[i]] = h;
            if (h > 0) h--;
        } else { lcp[0] = 0; h = 0; }
    }
    lcp[n] = 0;
    int cap = 64; Mum *m = malloc(sizeof(Mum) * cap); int k = 0;
    for (int i = 1; i < n; i++) {
        int x = sa[i - 1], y = sa[i];
        int xa = x < na, ya = y < na;
        if (xa == ya || x == na || y == na || x == n - 1 || y == n - 1) continue;
        int l = lcp[i];
        if (l > *lcs_out) *lcs_out = l;
        if (l < Lmin) continue;
        if (l <= lcp[i - 1] || l <= lcp[i + 1]) continue;  // shared by a third suffix: not unique
        int pa = xa ? x : y, pb = (xa ? y : x) - na - 1;
        if (pa > 0 && pb > 0 && s[pa - 1] == s[na + 1 + pb - 1] && s[pa - 1] != 6) continue;  // not left-maximal
        if (k == cap) { cap *= 2; m = realloc(m, sizeof(Mum) * cap); }
        m[k].pa = a0 + pa; m[k].pb = b0 + pb; m[k].len = l; k++;
    }
    free(s); free(sa); free(rk); free(tmp); free(lcp);
    *nout = k;
    return m;
}

static int ed(const char *A, int n, const char *B, int m) {  // global unit edit distance
    if (n == 0) return m;
    if (m == 0) return n;
    if (m > n) { const char *t = A; A = B; B = t; int u = n; n = m; m = u; }
    int *row = malloc(sizeof(int) * (m + 1));
    for (int j = 0; j <= m; j++) row[j] = j;
    for (int i = 1; i <= n; i++) {
        int diag = row[0]; row[0] = i; char a = A[i - 1];
        for (int j = 1; j <= m; j++) {
            int up = row[j];
            int v = diag + (a != B[j - 1]);
            int w = up + 1; if (w < v) v = w;
            w = row[j - 1] + 1; if (w < v) v = w;
            row[j] = v; diag = up;
        }
    }
    int r = row[m]; free(row); return r;
}

typedef struct { int p, q, idx; } Q;
static int cmp_q(const void *x, const void *y) { return ((const Q *)x)->p - ((const Q *)y)->p; }

// D(p,q) = ED(A[:p],B[:q]) for each query; out[idx]
static void ed_queries(const char *A, int n, const char *B, int m, Q *qs, int nq, int *out) {
    qsort(qs, nq, sizeof(Q), cmp_q);
    int *row = malloc(sizeof(int) * (m + 1));
    for (int j = 0; j <= m; j++) row[j] = j;
    int t = 0;
    while (t < nq && qs[t].p == 0) { out[qs[t].idx] = row[qs[t].q]; t++; }
    for (int i = 1; i <= n && t < nq; i++) {
        int diag = row[0]; row[0] = i; char a = A[i - 1];
        for (int j = 1; j <= m; j++) {
            int up = row[j];
            int v = diag + (a != B[j - 1]);
            int w = up + 1; if (w < v) v = w;
            w = row[j - 1] + 1; if (w < v) v = w;
            row[j] = v; diag = up;
        }
        while (t < nq && qs[t].p == i) { out[qs[t].idx] = row[qs[t].q]; t++; }
    }
    free(row);
}

static int cmp_mum(const void *x, const void *y) {
    const Mum *a = x, *b = y;
    if (a->pa != b->pa) return a->pa - b->pa;
    return a->pb - b->pb;
}

// heaviest chain; returns chain length (#mums) into idx[] (in order)
// score = sum of lengths - g_gapw * (sum of |diagonal shifts|, including from the segment's start
// diagonal and to its end diagonal); g_gapw = 0 is MUMmer-style heaviest chain.
static double g_gapw = 0.0; static int g_d0 = 0, g_d1 = 0;
static long iabs(long x) { return x < 0 ? -x : x; }
static int chain(Mum *m, int k, int *idx) {
    if (k == 0) return 0;
    qsort(m, k, sizeof(Mum), cmp_mum);
    double *best = malloc(sizeof(double) * k); int *prev = malloc(sizeof(int) * k);
    int bi = -1; double bs = 0;
    for (int j = 0; j < k; j++) {
        long dj = (long)m[j].pb - m[j].pa;
        best[j] = m[j].len - g_gapw * iabs(dj - g_d0); prev[j] = -1;
        for (int i = 0; i < j; i++)
            if (m[i].pa + m[i].len <= m[j].pa && m[i].pb + m[i].len <= m[j].pb) {
                long di = (long)m[i].pb - m[i].pa;
                double v = best[i] + m[j].len - g_gapw * iabs(dj - di);
                if (v > best[j]) { best[j] = v; prev[j] = i; }
            }
        double fin = best[j] - g_gapw * iabs((long)g_d1 - dj);
        if (bi < 0 || fin > bs) { bi = j; bs = fin; }
    }
    if (g_gapw > 0 && bs <= -g_gapw * iabs((long)g_d1 - g_d0)) { free(best); free(prev); return 0; }
    int c = 0;
    for (int j = bi; j >= 0; j = prev[j]) idx[c++] = j;
    for (int i = 0; i < c / 2; i++) { int t = idx[i]; idx[i] = idx[c - 1 - i]; idx[c - 1 - i] = t; }
    free(best); free(prev);
    return c;
}

static long rec_cells; static int rec_maxwin, rec_maxdepth, rec_leaves; static long rec_anchor_bp;

static int recurse(const char *A, int a0, int a1, const char *B, int b0, int b1, int Lmin, int depth) {
    if (depth > rec_maxdepth) rec_maxdepth = depth;
    int k, lcs;
    Mum *m = find_mums(A, a0, a1, B, b0, b1, Lmin, &k, &lcs);
    if (k == 0) {
        free(m);
        int w = (a1 - a0) > (b1 - b0) ? (a1 - a0) : (b1 - b0);
        if (w > rec_maxwin) rec_maxwin = w;
        rec_cells += (long)(a1 - a0) * (b1 - b0); rec_leaves++;
        return ed(A + a0, a1 - a0, B + b0, b1 - b0);
    }
    int *idx = malloc(sizeof(int) * k);
    g_d0 = b0 - a0; g_d1 = b1 - a1;
    int c = chain(m, k, idx);
    if (c == 0) {
        free(idx); free(m);
        int w = (a1 - a0) > (b1 - b0) ? (a1 - a0) : (b1 - b0);
        if (w > rec_maxwin) rec_maxwin = w;
        rec_cells += (long)(a1 - a0) * (b1 - b0); rec_leaves++;
        return ed(A + a0, a1 - a0, B + b0, b1 - b0);
    }
    int cost = 0, pa = a0, pb = b0;
    for (int t = 0; t < c; t++) {
        Mum u = m[idx[t]];
        cost += recurse(A, pa, u.pa, B, pb, u.pb, Lmin, depth + 1);
        rec_anchor_bp += u.len;
        pa = u.pa + u.len; pb = u.pb + u.len;
    }
    cost += recurse(A, pa, a1, B, pb, b1, Lmin, depth + 1);
    free(idx); free(m);
    return cost;
}

static char *read_line(FILE *f) {
    size_t cap = 1 << 16, len = 0; char *buf = malloc(cap); int ch;
    while ((ch = fgetc(f)) != EOF && ch != '\n') {
        if (len + 1 >= cap) { cap *= 2; buf = realloc(buf, cap); }
        buf[len++] = (char)ch;
    }
    if (ch == EOF && len == 0) { free(buf); return NULL; }
    buf[len] = 0; return buf;
}

int main(void) {
    if (getenv("CHAIN_GAPW")) g_gapw = atof(getenv("CHAIN_GAPW"));
    printf("id\tlenA\tlenB\topt\tlcs\tn_mum\tmum_bp\tn_core\tn_flank\tn_cross\tcore_bp\tflank_bp"
           "\tzero_core\tzero_flank\tzero_cross\texcess_core\texcess_flank\tchain_n\tchain_bp\tchain_core_n"
           "\tchain_core_bp\tchain_cost\trec_cost\trec_anchor_bp\trec_maxwin\trec_cells\trec_depth\trec_leaves\n");
    char *line;
    while ((line = read_line(stdin))) {
        char *f[9]; int nf = 0; char *p = line;
        f[nf++] = p;
        while (*p && nf < 9) { if (*p == '\t') { *p = 0; f[nf++] = p + 1; } p++; }
        if (nf < 9) { free(line); continue; }
        int Lmin = atoi(f[1]), ca0 = atoi(f[2]), ca1 = atoi(f[3]), doDP = atoi(f[6]);
        const char *A = f[7], *B = f[8]; int n = strlen(A), m = strlen(B);
        (void)atoi(f[4]); (void)atoi(f[5]);
        int k, lcs;
        Mum *mm = find_mums(A, 0, n, B, 0, m, Lmin, &k, &lcs);
        int ncore = 0, nflank = 0, ncross = 0; long corebp = 0, flankbp = 0, mumbp = 0;
        int *zone = malloc(sizeof(int) * (k + 1));
        for (int i = 0; i < k; i++) {
            int s = mm[i].pa, e = mm[i].pa + mm[i].len;
            mumbp += mm[i].len;
            if (s >= ca0 && e <= ca1) { zone[i] = 0; ncore++; corebp += mm[i].len; }
            else if (e <= ca0 || s >= ca1) { zone[i] = 1; nflank++; flankbp += mm[i].len; }
            else { zone[i] = 2; ncross++; }
        }
        int opt = -1, zero[3] = {0, 0, 0}; long exc[3] = {0, 0, 0};
        int chain_n = 0, chain_core_n = 0; long chain_bp = 0, chain_core_bp = 0; int chain_cost = -1;
        int rcost = -1;
        rec_cells = 0; rec_maxwin = 0; rec_maxdepth = 0; rec_leaves = 0; rec_anchor_bp = 0;
        if (doDP) {
            opt = ed(A, n, B, m);
            if (k > 0) {
                Q *qs = malloc(sizeof(Q) * k); int *pre = malloc(sizeof(int) * k), *suf = malloc(sizeof(int) * k);
                for (int i = 0; i < k; i++) { qs[i].p = mm[i].pa; qs[i].q = mm[i].pb; qs[i].idx = i; }
                ed_queries(A, n, B, m, qs, k, pre);
                char *rA = malloc(n + 1), *rB = malloc(m + 1);
                for (int i = 0; i < n; i++) rA[i] = A[n - 1 - i];
                for (int i = 0; i < m; i++) rB[i] = B[m - 1 - i];
                for (int i = 0; i < k; i++) { qs[i].p = n - mm[i].pa - mm[i].len; qs[i].q = m - mm[i].pb - mm[i].len; qs[i].idx = i; }
                ed_queries(rA, n, rB, m, qs, k, suf);
                for (int i = 0; i < k; i++) {
                    int x = pre[i] + suf[i] - opt;
                    exc[zone[i]] += x;
                    if (x == 0) zero[zone[i]]++;
                }
                free(qs); free(pre); free(suf); free(rA); free(rB);
                // chain cost (top level only)
                Mum *cp = malloc(sizeof(Mum) * k); memcpy(cp, mm, sizeof(Mum) * k);
                int *idx = malloc(sizeof(int) * k);
                g_d0 = 0; g_d1 = m - n;
                chain_n = chain(cp, k, idx);
                int pa = 0, pb = 0; chain_cost = 0;
                for (int t = 0; t < chain_n; t++) {
                    Mum u = cp[idx[t]];
                    chain_cost += ed(A + pa, u.pa - pa, B + pb, u.pb - pb);
                    chain_bp += u.len;
                    if (u.pa >= ca0 && u.pa + u.len <= ca1) { chain_core_n++; chain_core_bp += u.len; }
                    pa = u.pa + u.len; pb = u.pb + u.len;
                }
                chain_cost += ed(A + pa, n - pa, B + pb, m - pb);
                free(cp); free(idx);
            } else chain_cost = opt;
            rcost = recurse(A, 0, n, B, 0, m, Lmin, 0);
        }
        printf("%s\t%d\t%d\t%d\t%d\t%d\t%ld\t%d\t%d\t%d\t%ld\t%ld\t%d\t%d\t%d\t%ld\t%ld\t%d\t%ld\t%d\t%ld\t%d\t%d\t%ld\t%d\t%ld\t%d\t%d\n",
               f[0], n, m, opt, lcs, k, mumbp, ncore, nflank, ncross, corebp, flankbp, zero[0], zero[1], zero[2],
               exc[0], exc[1], chain_n, chain_bp, chain_core_n, chain_core_bp, chain_cost, rcost, rec_anchor_bp,
               rec_maxwin, rec_cells, rec_maxdepth, rec_leaves);
        fflush(stdout);
        free(mm); free(zone); free(line);
    }
    return 0;
}
