/* Exact complete central-piece rank in the characteristic-two homogeneous family.
 * The source crossing and quotient are attributed in research/modular-central-rank.md.
 * Use tools/check_modular_central_rank.partsy to generate the finite input and audit it.
 * This is a finite research computation, not an asymptotic factoring bound.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
typedef struct {
    int length;
    int parts[16];
} Partition;
static int determinant_two(uint16_t *rows, int n) {
    for (int c = 0; c < n; c++) {
        int p = c;
        while (p < n && !((rows[p] >> c) & 1))
            p++;
        if (p == n)
            return 0;
        uint16_t t = rows[p];
        rows[p] = rows[c];
        rows[c] = t;
        for (int i = c + 1; i < n; i++)
            if ((rows[i] >> c) & 1)
                rows[i] ^= rows[c];
    }
    return 1;
}
static int coefficient_two(Partition a, Partition out, Partition in) {
    int twice = out.length + in.length - a.length;
    if (twice < 0 || twice % 2)
        return 0;
    int remain = twice / 2, common[16], nc = 0;
    for (int i = 0; i < out.length; i++)
        for (int j = 0; j < in.length; j++)
            if (out.parts[i] == in.parts[j])
                common[nc++] = out.parts[i];
    if (remain > nc)
        return 0;
    int result = 0;
    for (int mask = 0; mask < (1 << nc); mask++) {
        if (__builtin_popcount((unsigned)mask) != remain)
            continue;
        int signedmode[16], create[16], n = 0;
        for (int which = 0; which < 2; which++) {
            Partition q = which ? in : out;
            for (int i = 0; i < q.length; i++) {
                int skip = 0;
                for (int j = 0; j < nc; j++)
                    if (((mask >> j) & 1) && q.parts[i] == common[j])
                        skip = 1;
                if (!skip) {
                    signedmode[n] = q.parts[i];
                    create[n++] = !which;
                }
            }
        }
        if (n != a.length)
            abort();
        uint16_t rows[16];
        for (int i = 0; i < a.length; i++) {
            rows[i] = 0;
            for (int j = 0; j < n; j++) {
                int m = a.parts[i], k = signedmode[j];
                int b = create[j] ? (((m - 1) & ~(k - 1)) == 0) : (((m - 1) & k) == 0);
                rows[i] |= b << j;
            }
        }
        result ^= determinant_two(rows, n);
    }
    return result;
}
static Partition *read_parts(FILE *f, int *n) {
    fscanf(f, "%d", n);
    Partition *p = calloc(*n, sizeof(Partition));
    for (int i = 0; i < *n; i++) {
        fscanf(f, "%d", &p[i].length);
        for (int j = 0; j < p[i].length; j++)
            fscanf(f, "%d", &p[i].parts[j]);
    }
    return p;
}
static int insert_basis(uint64_t *row, uint64_t **basis, int words) {
    for (int w = words - 1; w >= 0; w--) {
        while (row[w]) {
            int b = 63 - __builtin_clzll(row[w]);
            int p = 64 * w + b;
            if (!basis[p]) {
                basis[p] = malloc((w + 1) * sizeof(uint64_t));
                memcpy(basis[p], row, (w + 1) * sizeof(uint64_t));
                return 1;
            }
            for (int i = 0; i <= w; i++)
                row[i] ^= basis[p][i];
        }
    }
    return 0;
}

static uint64_t **orthogonal_complement(uint64_t **basis, int columns, int words,
                                        int *dimension) {
    uint64_t **out = calloc(columns, sizeof(uint64_t *));
    *dimension = 0;
    for (int f = 0; f < columns; f++)
        if (!basis[f]) {
            uint64_t *v = calloc(words, sizeof(uint64_t));
            v[f / 64] |= 1ULL << (f % 64);
            for (int p = 0; p < columns; p++)
                if (basis[p]) {
                    int parity = 0;
                    for (int j = 0; j <= p / 64; j++)
                        parity ^= __builtin_parityll(basis[p][j] & v[j]);
                    if (parity)
                        v[p / 64] ^= 1ULL << (p % 64);
                }
            out[(*dimension)++] = v;
        }
    return out;
}
static int kernel_pair_rank(uint64_t **left, int nl, uint64_t **right, int nr, int words) {
    int rw = (nr + 63) / 64, rank = 0;
    if (!nr || !nl)
        return 0;
    uint64_t **p = calloc(nr, sizeof(uint64_t *));
    uint64_t *row = calloc(rw, sizeof(uint64_t));
    for (int i = 0; i < nl; i++) {
        memset(row, 0, rw * 8);
        for (int j = 0; j < nr; j++) {
            int parity = 0;
            for (int w = 0; w < words; w++)
                parity ^= __builtin_parityll(left[i][w] & right[j][w]);
            row[j / 64] |= (uint64_t)parity << (j % 64);
        }
        rank += insert_basis(row, p, rw);
    }
    return rank;
}
int main(int argc, char **argv) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "usage: %s INPUT [--coefficients]\n", argv[0]);
        return 2;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f)
        return 2;
    int d, q;
    Partition *lab = read_parts(f, &d);
    Partition *mid = read_parts(f, &q);
    fclose(f);
    if (argc == 3) {
        if (strcmp(argv[2], "--coefficients"))
            return 2;
        int a, y, j;
        while (scanf("%d %d %d", &a, &y, &j) == 3) {
            if (a < 0 || a >= d || y < 0 || y >= d || j < 0 || j >= q)
                return 2;
            Partition empty = {0};
            printf("%d %d %d %d\n", coefficient_two(lab[a], mid[j], lab[y]),
                   coefficient_two(lab[a], lab[y], mid[j]),
                   coefficient_two(lab[a], empty, lab[y]), a == y);
        }
        return 0;
    }
    int columns = q + 1;
    int words = (columns + 63) / 64;
    uint64_t **lb = calloc(columns, sizeof(uint64_t *)),
             **rb = calloc(columns, sizeof(uint64_t *));
    uint64_t *l = calloc(words, sizeof(uint64_t)), *r = calloc(words, sizeof(uint64_t));
    int lr = 0, rr = 0;
    int target = columns;
    clock_t begin = clock();
    for (int i = 0; i < d * d; i++) {
        memset(l, 0, words * 8);
        memset(r, 0, words * 8);
        int a = i / d, y = i % d;
        for (int j = 0; j < q; j++) {
            l[j / 64] |= (uint64_t)coefficient_two(lab[a], mid[j], lab[y]) << (j % 64);
            r[j / 64] |= (uint64_t)coefficient_two(lab[a], lab[y], mid[j]) << (j % 64);
        }
        Partition empty = {0};
        l[q / 64] |= (uint64_t)coefficient_two(lab[a], empty, lab[y]) << (q % 64);
        r[q / 64] |= (uint64_t)(a == y) << (q % 64);
        lr += insert_basis(l, lb, words);
        rr += insert_basis(r, rb, words);
        if ((i + 1) % 128 == 0 || (lr >= target && rr >= target) || i + 1 == d * d) {
            printf("{\"pairs\":%d,\"labels\":%d,\"channels\":%d,\"left_rank\":%d,\"right_"
                   "rank\":%d,\"cpu_seconds\":%.2f}\n",
                   i + 1, d, q, lr, rr, (double)(clock() - begin) / CLOCKS_PER_SEC);
            fflush(stdout);
        }
        if (lr >= target && rr >= target)
            break;
    }
    int kl, kr;
    uint64_t **lk = orthogonal_complement(lb, columns, words, &kl),
             **rk = orthogonal_complement(rb, columns, words, &kr);
    int small = kernel_pair_rank(lk, kl, rk, kr, words);
    printf("{\"kernel_left\":%d,\"kernel_right\":%d,\"kernel_pair_rank\":%d,\"central_rank\":%"
           "d,\"total_pairs\":%d}\n",
           kl, kr, small, lr + rr - columns + small, d * d);
    return 0;
}
