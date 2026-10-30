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
    if (fscanf(f, "%d", n) != 1 || *n < 1)
        exit(2);
    Partition *p = calloc(*n, sizeof(Partition));
    for (int i = 0; i < *n; i++) {
        if (fscanf(f, "%d", &p[i].length) != 1 || p[i].length < 0 || p[i].length > 16)
            exit(2);
        for (int j = 0; j < p[i].length; j++)
            if (fscanf(f, "%d", &p[i].parts[j]) != 1)
                exit(2);
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
static int binary_prefix_rank(FILE *f, Partition *labels, int d, Partition *states, int q) {
    int selected;
    if (fscanf(f, "%d", &selected) != 1 || selected < 1 || selected > d)
        return 2;
    int *counts = calloc(selected, sizeof(int));
    int **dual = calloc(selected, sizeof(int *));
    for (int i = 0; i < selected; i++) {
        if (fscanf(f, "%d", &counts[i]) != 1 || counts[i] < 1 || counts[i] > d)
            return 2;
        dual[i] = calloc(counts[i], sizeof(int));
        for (int j = 0; j < counts[i]; j++)
            if (fscanf(f, "%d", &dual[i][j]) != 1 || dual[i][j] < 0 || dual[i][j] >= d)
                return 2;
    }
    int channels = q + 1, words = (channels + 63) / 64;
    uint64_t *left = calloc(d * words, 8), *right = calloc(d * words, 8),
             *scratch = calloc(words, 8);
    uint64_t **basis[4];
    int ranks[4] = {0};
    for (int k = 0; k < 4; k++)
        basis[k] = calloc(channels, sizeof(uint64_t *));
    clock_t begun = clock();
    Partition empty = {0};
    for (int y = 0; y < selected; y++) {
        memset(left, 0, d * words * 8);
        memset(right, 0, d * words * 8);
        for (int a = 0; a < d; a++) {
            uint64_t *l = left + a * words, *r = right + a * words;
            for (int j = 0; j < q; j++) {
                l[j / 64] |= (uint64_t)coefficient_two(labels[a], states[j], labels[y])
                             << (j % 64);
                r[j / 64] |= (uint64_t)coefficient_two(labels[a], labels[y], states[j])
                             << (j % 64);
            }
            l[q / 64] |= (uint64_t)coefficient_two(labels[a], empty, labels[y]) << (q % 64);
            r[q / 64] |= (uint64_t)(a == y) << (q % 64);
        }
        for (int s = 0; s < selected; s++) {
            memcpy(scratch, right + s * words, words * 8);
            ranks[1] += insert_basis(scratch, basis[1], words);
            memcpy(scratch, left + s * words, words * 8);
            ranks[3] += insert_basis(scratch, basis[3], words);
            memset(scratch, 0, words * 8);
            for (int i = 0; i < counts[s]; i++) {
                uint64_t *l = left + dual[s][i] * words;
                for (int j = 0; j < words; j++)
                    scratch[j] ^= l[j];
            }
            ranks[0] += insert_basis(scratch, basis[0], words);
            memset(scratch, 0, words * 8);
            for (int i = 0; i < counts[s]; i++) {
                uint64_t *r = right + dual[s][i] * words;
                for (int j = 0; j < words; j++)
                    scratch[j] ^= r[j];
            }
            ranks[2] += insert_basis(scratch, basis[2], words);
        }
        if ((y + 1) % 16 == 0 || y + 1 == selected) {
            printf("{\"states\":%d,\"selected\":%d,\"ranks\":[%d,%d,%d,%d],\"cpu_seconds\":%."
                   "2f}\n",
                   y + 1, selected, ranks[0], ranks[1], ranks[2], ranks[3],
                   (double)(clock() - begun) / CLOCKS_PER_SEC);
            fflush(stdout);
        }
    }
    for (int orientation = 0; orientation < 2; orientation++) {
        int a = 2 * orientation, b = a + 1, ka, kb;
        uint64_t **k1 = orthogonal_complement(basis[a], channels, words, &ka),
                 **k2 = orthogonal_complement(basis[b], channels, words, &kb);
        int r = kernel_pair_rank(k1, ka, k2, kb, words);
        printf("{\"orientation\":%d,\"kernel_dimensions\":[%d,%d],\"kernel_pair_rank\":%d,"
               "\"complete_rank\":%d}\n",
               orientation, ka, kb, r, ranks[a] + ranks[b] - channels + r);
    }
    return 0;
}
int main(int argc, char **argv) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "usage: %s INPUT [--coefficients|--binary-prefix]\n", argv[0]);
        return 2;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f)
        return 2;
    int d, q;
    Partition *lab = read_parts(f, &d);
    Partition *mid = read_parts(f, &q);
    if (argc == 3 && !strcmp(argv[2], "--binary-prefix")) {
        int status = binary_prefix_rank(f, lab, d, mid, q);
        fclose(f);
        return status;
    }
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
