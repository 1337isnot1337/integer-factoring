// Independent scalar GF(256) determinant checker for 16-by-16 witnesses.
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <utility>

static uint8_t multiply[256][256];

static uint8_t bit_multiply(unsigned a, unsigned b) {
    unsigned result = 0;
    while (b) {
        if (b & 1) result ^= a;
        b >>= 1;
        a <<= 1;
        if (a & 256) a ^= 0x11d;
    }
    return static_cast<uint8_t>(result);
}

static uint8_t power(uint8_t a, unsigned exponent) {
    uint8_t result = 1;
    while (exponent) {
        if (exponent & 1) result = multiply[result][a];
        a = multiply[a][a];
        exponent >>= 1;
    }
    return result;
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const unsigned count = static_cast<unsigned>(std::strtoul(argv[1], nullptr, 10));
    for (unsigned a = 0; a < 256; ++a)
        for (unsigned b = 0; b < 256; ++b)
            multiply[a][b] = bit_multiply(a, b);
    for (unsigned a = 1; a < 256; ++a)
        if (multiply[a][power(static_cast<uint8_t>(a), 254)] != 1) return 3;
    for (unsigned example = 0; example < count; ++example) {
        uint8_t matrix[16][16];
        if (std::fread(matrix, 1, sizeof(matrix), stdin) != sizeof(matrix)) return 4;
        for (unsigned column = 0; column < 16; ++column) {
            unsigned pivot = column;
            while (pivot < 16 && !matrix[pivot][column]) ++pivot;
            if (pivot == 16) {
                std::fprintf(stderr, "Singular witness %u at column %u\n", example, column);
                return 5;
            }
            if (pivot != column)
                for (unsigned j = 0; j < 16; ++j)
                    std::swap(matrix[pivot][j], matrix[column][j]);
            const uint8_t inverse = power(matrix[column][column], 254);
            for (unsigned j = column; j < 16; ++j)
                matrix[column][j] = multiply[inverse][matrix[column][j]];
            for (unsigned i = column + 1; i < 16; ++i) {
                const uint8_t factor = matrix[i][column];
                if (!factor) continue;
                for (unsigned j = column; j < 16; ++j)
                    matrix[i][j] ^= multiply[factor][matrix[column][j]];
            }
        }
    }
    if (std::getc(stdin) != EOF) return 6;
    std::printf("{\"independent_nonsingular_minors\":%u,\"minor_order\":16,\"field_order\":256}\n", count);
    return 0;
}
