/* Finite full-span search for the Hamming Pfaffian family at at most eight labels.
 * Exact F_2 ranks; threshold clipping cannot remove a better decomposition.
 * The uniform full-span construction is described in research/hamming-pfaffian.md.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct {
    int cost;
    uint32_t index;
} Record;
int d, n, rows, threshold;
uint64_t gen[64][64];
int rank_cap(uint64_t *a, int cap) {
    uint64_t p[64] = {0};
    int r = 0;
    for (int i = 0; i < rows; i++) {
        uint64_t v = a[i];
        while (v) {
            int b = 63 - __builtin_clzll(v);
            if (p[b])
                v ^= p[b];
            else {
                p[b] = v;
                if (++r == cap)
                    return r;
                break;
            }
        }
    }
    return r;
}
int cost(uint64_t *a) {
    int r = rank_cap(a, rows);
    if (!r)
        return 0;
    uint64_t c[64] = {0}, mask = (1ULL << d) - 1;
    for (int x = 0; x < d; x++)
        for (int u = 0; u < d; u++)
            for (int y = 0; y < d; y++)
                c[x * d + u] |= ((a[x * d + y] >> (u * d)) & mask) << (y * d);
    int cap = threshold / r + 1;
    if (cap > rows)
        cap = rows;
    int k = rank_cap(c, cap);
    return r * k <= threshold ? r * k : 0;
}
int cmp(const void *a, const void *b) {
    const Record *x = a, *y = b;
    return x->cost - y->cost;
}
void value_at(uint64_t *a, uint32_t i) {
    memset(a, 0, rows * 8);
    uint64_t g = i ^ (i >> 1);
    for (int j = 0; j < n; j++)
        if ((g >> j) & 1)
            for (int k = 0; k < rows; k++)
                a[k] ^= gen[j][k];
}
int high(uint64_t *a) {
    for (int j = rows - 1; j >= 0; j--)
        if (a[j])
            return j * rows + 63 - __builtin_clzll(a[j]);
    return -1;
}
int main(int argc, char **argv) {
    if (argc < 3 || argc > 4)
        return 2;
    FILE *f = fopen(argv[1], "r");
    if (!f)
        return 2;
    if (fscanf(f, "%d %d", &d, &n) != 2)
        return 2;
    rows = d * d;
    if (d < 2 || d > 8 || n < 1 || n > 31)
        return 2;
    threshold = atoi(argv[2]);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < rows; j++) {
            unsigned long long v;
            if (fscanf(f, "%llx", &v) != 1)
                return 2;
            gen[i][j] = v;
        }
    fclose(f);
    if (argc == 4) {
        if (strcmp(argv[3], "--classify"))
            return 2;
        uint32_t index;
        uint64_t a[64];
        while (scanf("%u", &index) == 1) {
            if (index >= (1U << n))
                return 2;
            value_at(a, index);
            printf("%d\n", cost(a));
        }
        return 0;
    }
    uint64_t v[64] = {0};
    int count = 0, capacity = 1024;
    Record *records = malloc(capacity * sizeof(Record));
    uint32_t limit = 1U << n;
    for (uint32_t i = 1; i < limit; i++) {
        int bit = __builtin_ctz(i);
        for (int j = 0; j < rows; j++)
            v[j] ^= gen[bit][j];
        int c = cost(v);
        if (c) {
            if (count == capacity) {
                capacity *= 2;
                records = realloc(records, capacity * sizeof(Record));
            }
            records[count++] = (Record){c, i};
        }
        if ((i & ((1U << 18) - 1)) == 0) {
            printf("tested %u accepted %d\n", i, count);
            fflush(stdout);
        }
    }
    qsort(records, count, sizeof(Record), cmp);
    uint64_t *basis[4096] = {0}, masks[4096] = {0}, selected[64];
    int rank = 0;
    uint64_t target[64] = {0};
    for (int x = 0; x < d; x++)
        for (int y = 0; y < d; y++)
            target[x * d + y] = 1ULL << (y * d + x);
    for (int i = 0; i < count; i++) {
        value_at(v, records[i].index);
        uint64_t mask = 1ULL << rank;
        int p;
        while ((p = high(v)) >= 0) {
            if (basis[p]) {
                for (int j = 0; j < rows; j++)
                    v[j] ^= basis[p][j];
                mask ^= masks[p];
            } else {
                basis[p] = malloc(rows * 8);
                memcpy(basis[p], v, rows * 8);
                masks[p] = mask;
                selected[rank++] = records[i].index;
                break;
            }
        }
        uint64_t a[64];
        memcpy(a, target, rows * 8);
        mask = 0;
        while ((p = high(a)) >= 0) {
            if (!basis[p])
                break;
            for (int j = 0; j < rows; j++)
                a[j] ^= basis[p][j];
            mask ^= masks[p];
        }
        if (p < 0) {
            printf("FOUND cutoff %d span %d accepted %d\n", records[i].cost, rank, count);
            for (int j = 0; j < rank; j++)
                if ((mask >> j) & 1)
                    printf("piece %llu\n", (unsigned long long)selected[j]);
            return 0;
        }
    }
    printf("NONE threshold %d span %d accepted %d\n", threshold, rank, count);
    return 0;
}
