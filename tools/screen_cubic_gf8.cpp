#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

struct Row { std::array<std::uint64_t, 3> p{}; };

static unsigned multiply(unsigned a, unsigned b) {
    unsigned result = 0;
    while (b) {
        if (b & 1U) result ^= a;
        b >>= 1U;
        a <<= 1U;
        if (a & 8U) a ^= 11U; // x^3+x+1
    }
    return result;
}

static unsigned inverse(unsigned a) {
    for (unsigned b = 1; b < 8; ++b) if (multiply(a, b) == 1U) return b;
    throw std::runtime_error("zero has no inverse");
}

static Row scaled(const Row& row, unsigned scalar) {
    Row result;
    for (unsigned i = 0; i < 3; ++i) {
        const unsigned value = multiply(scalar, 1U << i);
        for (unsigned j = 0; j < 3; ++j) if ((value >> j) & 1U) result.p[j] ^= row.p[i];
    }
    return result;
}

static void add(Row& target, const Row& source) {
    for (unsigned i = 0; i < 3; ++i) target.p[i] ^= source.p[i];
}

static unsigned rank(std::vector<Row> rows, unsigned cap) {
    std::array<std::array<Row, 8>, 64> basis{};
    std::uint64_t occupied = 0;
    unsigned result = 0;
    for (auto row : rows) {
        while (const auto support = row.p[0] | row.p[1] | row.p[2]) {
            const unsigned pivot = 63U - std::countl_zero(support);
            const unsigned value = ((row.p[0] >> pivot) & 1U)
                | (((row.p[1] >> pivot) & 1U) << 1U)
                | (((row.p[2] >> pivot) & 1U) << 2U);
            if (!(occupied & (std::uint64_t{1} << pivot))) {
                ++result;
                if (result > cap) return result;
                occupied |= std::uint64_t{1} << pivot;
                const auto inv = inverse(value);
                for (unsigned coefficient = 1; coefficient < 8; ++coefficient)
                    basis[pivot][coefficient] = scaled(row, multiply(coefficient, inv));
                break;
            }
            add(row, basis[pivot][value]);
        }
    }
    return result;
}

int main() {
    try {
        unsigned labels, dimension, gate;
        if (!(std::cin >> labels >> dimension >> gate)) throw std::runtime_error("missing header");
        const unsigned width = labels * labels;
        if (!labels || width > 64 || !dimension || dimension > 8 || !gate)
            throw std::runtime_error("invalid dimensions");
        using Matrix = std::vector<Row>;
        using Multiples = std::array<Matrix, 8>;
        std::array<std::vector<Multiples>, 2> generators;
        for (auto& cut : generators) {
            cut.resize(dimension);
            for (auto& generator : cut) {
                Matrix matrix(width);
                for (auto& row : matrix) for (auto& plane : row.p)
                    if (!(std::cin >> std::hex >> plane)) throw std::runtime_error("truncated matrix");
                for (unsigned scalar = 0; scalar < 8; ++scalar) {
                    generator[scalar].reserve(width);
                    for (const auto& row : matrix) generator[scalar].push_back(scaled(row, scalar));
                }
            }
        }
        std::uint64_t tested = 0, eligible = 0;
        const unsigned cap = gate * gate - 1;
        std::cout << "{\"event\":\"start\",\"labels\":" << labels
                  << ",\"dimension\":" << dimension << ",\"field_size\":8,\"field_polynomial\":11"
                  << ",\"integer_cost_gate\":" << gate << ",\"product_cap\":" << cap << "}\n";
        for (unsigned pivot = 0; pivot < dimension; ++pivot) {
            const std::uint64_t leading = std::uint64_t{1} << (3U * pivot);
            for (std::uint64_t tail = 0; tail < leading; ++tail) {
                const std::uint64_t code = leading | tail;
                auto combine = [&](unsigned cut) {
                    Matrix matrix(width);
                    for (unsigned i = 0; i < dimension; ++i) {
                        const unsigned scalar = (code >> (3U * i)) & 7U;
                        if (!scalar) continue;
                        for (unsigned j = 0; j < width; ++j) add(matrix[j], generators[cut][i][scalar][j]);
                    }
                    return matrix;
                };
                ++tested;
                const unsigned row_rank = rank(combine(0), width);
                const unsigned column_rank = rank(combine(1), row_rank ? cap / row_rank : width);
                const unsigned value = row_rank * column_rank;
                if (value <= cap) {
                    ++eligible;
                    std::cout << "{\"event\":\"eligible\",\"code\":" << code
                              << ",\"row_rank\":" << row_rank << ",\"column_rank\":" << column_rank
                              << ",\"product\":" << value << "}\n";
                }
            }
        }
        std::cout << "{\"event\":\"complete\",\"projective_tested\":" << tested
                  << ",\"eligible\":" << eligible << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
