/* Exhaustive finite screen of a specified F_2 subspace of 9 by 9 matrices.
 * Each row is a (x,y), each column a (u,v), with all labels in {0,1,2}.
 * Enumerates a Gray-code interval and retains only rank products <= 9.
 * This is a finite research computation, not an asymptotic theorem.
 * Compile: cc -O3 -std=c11 -o /tmp/search-three-label tools/search_three_label_span.c
 */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef __uint128_t Vector;

static int high_bit(Vector value) {
    uint64_t upper = (uint64_t)(value >> 64);
    return upper ? 127 - __builtin_clzll(upper)
                 : 63 - __builtin_clzll((uint64_t)value);
}

/* Return min(actual rank, cap). No floating-point arithmetic is used. */
static int rank_cap(const uint16_t rows[9], int cap) {
    uint16_t pivots[9] = {0};
    int rank = 0;
    for (int i = 0; i < 9; ++i) {
        unsigned row = rows[i];
        while (row) {
            int pivot = 31 - __builtin_clz(row);
            if (pivots[pivot]) {
                row ^= pivots[pivot];
            } else {
                pivots[pivot] = (uint16_t)row;
                if (++rank == cap) return rank;
                break;
            }
        }
    }
    return rank;
}

static int low_product(Vector vector) {
    uint16_t rows[9], columns[9];
    for (int i = 0; i < 9; ++i) rows[i] = (vector >> (9*i)) & 511;
    int row_rank = rank_cap(rows, 5);
    if (!row_rank) return 0;
    for (int x = 0; x < 3; ++x) {
        for (int u = 0; u < 3; ++u) {
            columns[3*x+u] = ((rows[3*x] >> (3*u)) & 7)
                | (((rows[3*x+1] >> (3*u)) & 7) << 3)
                | (((rows[3*x+2] >> (3*u)) & 7) << 6);
        }
    }
    /* A row rank >=5 permits only column rank 1; its actual row rank
     * is at most 9. For row ranks 2,3,4 the first failing column ranks
     * are 5,4,3. These caps classify product <=9 exactly. */
    int cap = row_rank >= 5 ? 2 : row_rank == 4 ? 3
            : row_rank == 3 ? 4 : row_rank == 2 ? 5 : 9;
    return row_rank * rank_cap(columns, cap) <= 9;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s BASIS_HEX START END\n", argv[0]);
        return 2;
    }
    FILE *input = fopen(argv[1], "r");
    if (!input) return 2;
    Vector generators[63];
    uint64_t lower, upper;
    int dimension = 0;
    while (fscanf(input, "%" SCNx64 " %" SCNx64, &lower, &upper) == 2) {
        if (dimension == 63 || upper >= (UINT64_C(1) << 17)) return 2;
        generators[dimension++] = ((Vector)upper << 64) | lower;
    }
    fclose(input);
    uint64_t start = strtoull(argv[2], NULL, 0);
    uint64_t end = strtoull(argv[3], NULL, 0);
    if (start >= end || end > (UINT64_C(1) << dimension)) return 2;
    Vector pivots[81] = {0};
    int accepted_rank = 0;
    Vector vector = 0;
    uint64_t gray = start ^ (start >> 1);
    for (int i = 0; i < dimension; ++i) {
        if ((gray >> i) & 1) vector ^= generators[i];
    }
    uint64_t accepted = 0;
    clock_t begun = clock();
    printf("{\"dimension\":%d,\"start\":%" PRIu64 ",\"end\":%" PRIu64 "}\n",
           dimension, start, end);
    for (uint64_t index = start; index < end; ++index) {
        if (index != start) vector ^= generators[__builtin_ctzll(index)];
        if (low_product(vector)) {
            ++accepted;
            Vector reduced = vector;
            while (reduced) {
                int pivot = high_bit(reduced);
                if (pivots[pivot]) {
                    reduced ^= pivots[pivot];
                } else {
                    pivots[pivot] = reduced;
                    ++accepted_rank;
                    printf("{\"new_basis\":\"%" PRIx64 "%016" PRIx64
                           "\",\"index\":%" PRIu64 "}\n",
                           (uint64_t)(vector >> 64), (uint64_t)vector, index);
                    break;
                }
            }
        }
        if (((index+1) & ((UINT64_C(1) << 28)-1)) == 0) {
            printf("{\"checked_through\":%" PRIu64 ",\"low_product_count\":%" PRIu64
                   ",\"span_rank\":%d,\"cpu_seconds\":%.2f}\n",
                   index+1, accepted, accepted_rank,
                   (double)(clock()-begun)/CLOCKS_PER_SEC);
            fflush(stdout);
        }
    }
    printf("{\"complete\":true,\"checked\":%" PRIu64
           ",\"low_product_count\":%" PRIu64 ",\"span_rank\":%d}\n",
           end-start, accepted, accepted_rank);
    return 0;
}
