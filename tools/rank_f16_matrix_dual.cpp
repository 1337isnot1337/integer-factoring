#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

uint8_t multiply(uint8_t a, uint8_t b) {
    uint8_t result = 0;
    while (b != 0) {
        if (b & 1) result ^= a;
        a <<= 1;
        if (a & 0x10) a ^= 0x13;
        b >>= 1;
    }
    return result;
}

struct Field {
    uint8_t product[16][16] = {};
    uint8_t inverse[16] = {};
    uint8_t table[16][256] = {};

    Field() {
        for (int a = 0; a < 16; ++a)
            for (int b = 0; b < 16; ++b) product[a][b] = multiply(a, b);
        for (int a = 1; a < 16; ++a)
            for (int b = 1; b < 16; ++b)
                if (product[a][b] == 1) inverse[a] = static_cast<uint8_t>(b);
        for (int a = 0; a < 16; ++a)
            for (int v = 0; v < 256; ++v)
                table[a][v] = static_cast<uint8_t>(
                    product[a][v & 15] | (product[a][v >> 4] << 4));
    }
};

bool read_exact(FILE* file, void* data, size_t size) {
    return fread(data, 1, size, file) == size;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr,
                     "usage: %s MATRIX TARGET DUAL_OUTPUT\n"
                     "MATRIX is uint64 rows, uint64 columns, row-major bytes.\n",
                     argv[0]);
        return 2;
    }
    const Field field;
    FILE* matrix_file = std::fopen(argv[1], "rb");
    if (matrix_file == nullptr) return 2;
    uint64_t rows = 0, columns = 0;
    if (!read_exact(matrix_file, &rows, sizeof(rows)) ||
        !read_exact(matrix_file, &columns, sizeof(columns))) {
        return 2;
    }
    std::vector<uint8_t> matrix(static_cast<size_t>(rows * columns));
    if (!read_exact(matrix_file, matrix.data(), matrix.size())) return 2;
    std::fclose(matrix_file);

    FILE* target_file = std::fopen(argv[2], "rb");
    if (target_file == nullptr) return 2;
    std::vector<uint8_t> target(static_cast<size_t>(columns));
    if (!read_exact(target_file, target.data(), target.size())) return 2;
    std::fclose(target_file);

    std::vector<size_t> pivots;
    size_t rank = 0;
    for (size_t column = 0; column < columns && rank < rows; ++column) {
        size_t pivot = rank;
        while (pivot < rows && matrix[pivot * columns + column] == 0) ++pivot;
        if (pivot == rows) continue;
        if (pivot != rank) {
            for (size_t j = 0; j < columns; ++j)
                std::swap(matrix[pivot * columns + j],
                          matrix[rank * columns + j]);
        }
        const uint8_t scale = field.inverse[matrix[rank * columns + column]];
        if (scale != 1) {
            const uint8_t* table = field.table[scale];
            for (size_t j = 0; j < columns; ++j)
                matrix[rank * columns + j] = table[matrix[rank * columns + j]];
        }
        const uint8_t* pivot_row = &matrix[rank * columns];
        for (size_t row = rank + 1; row < rows; ++row) {
            const uint8_t factor = matrix[row * columns + column];
            if (factor == 0) continue;
            const uint8_t* table = field.table[factor];
            uint8_t* row_data = &matrix[row * columns];
            for (size_t j = 0; j < columns; ++j)
                row_data[j] ^= table[pivot_row[j]];
        }
        pivots.push_back(column);
        ++rank;
    }

    std::vector<uint8_t> is_pivot(static_cast<size_t>(columns), 0);
    for (size_t column : pivots) is_pivot[column] = 1;
    std::vector<uint8_t> dual(static_cast<size_t>(columns), 0);
    size_t free_column = columns;
    uint8_t target_dot = 0;
    for (size_t candidate = 0; candidate < columns; ++candidate) {
        if (is_pivot[candidate]) continue;
        std::fill(dual.begin(), dual.end(), 0);
        dual[candidate] = 1;
        for (size_t p = rank; p-- > 0;) {
            const size_t column = pivots[p];
            const uint8_t* row = &matrix[p * columns];
            uint8_t value = 0;
            for (size_t j = column + 1; j < columns; ++j)
                if (dual[j]) value ^= field.product[row[j]][dual[j]];
            dual[column] = value;
        }
        uint8_t dot = 0;
        for (size_t j = 0; j < columns; ++j)
            if (dual[j]) dot ^= field.product[target[j]][dual[j]];
        if (dot != 0) {
            free_column = candidate;
            target_dot = dot;
            break;
        }
    }
    if (free_column == columns) {
        std::printf("no separating functional\n");
        return 1;
    }
    for (size_t row = 0; row < rows; ++row) {
        const uint8_t* row_data = &matrix[row * columns];
        uint8_t dot = 0;
        for (size_t j = 0; j < columns; ++j)
            if (row_data[j] && dual[j])
                dot ^= field.product[row_data[j]][dual[j]];
        if (dot != 0) {
            std::printf("verification failed at row %zu\n", row);
            return 2;
        }
    }
    FILE* output_file = std::fopen(argv[3], "wb");
    if (output_file == nullptr) return 2;
    if (fwrite(dual.data(), 1, dual.size(), output_file) != dual.size()) return 2;
    std::fclose(output_file);
    std::printf("rank=%zu free_column=%zu target_dot=%u\n", rank, free_column,
                static_cast<unsigned>(target_dot));
    return 0;
}
