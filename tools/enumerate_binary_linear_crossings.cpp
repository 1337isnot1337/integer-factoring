// Exact Smith profiles of every scalar binary degree-one crossing.
// Polynomial bits represent F_2[t]; every minor is evaluated literally.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using Matrix = std::array<unsigned, 16>;
using Orders = std::array<int, 4>;

unsigned polynomial_product(unsigned a, unsigned b) {
    unsigned result = 0;
    for (; b; b >>= 1, a <<= 1)
        if (b & 1) result ^= a;
    return result;
}

unsigned minor_determinant(const Matrix& matrix, const std::vector<int>& rows,
                           const std::vector<int>& columns) {
    std::vector<int> permutation = columns;
    unsigned result = 0;
    do {
        unsigned term = 1;
        for (unsigned i = 0; i < rows.size(); ++i)
            term = polynomial_product(term, matrix[4 * rows[i] + permutation[i]]);
        result ^= term;  // Determinant signs are all one over F_2.
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return result;
}

std::vector<int> support(unsigned mask) {
    std::vector<int> result;
    for (int i = 0; i < 4; ++i)
        if (mask >> i & 1) result.push_back(i);
    return result;
}

Orders smith_orders(const Matrix& matrix) {
    Orders orders{};
    int previous = 0;
    for (int size = 1; size <= 4; ++size) {
        int best = 99;
        for (unsigned row_mask = 1; row_mask < 16; ++row_mask) {
            if (__builtin_popcount(row_mask) != size) continue;
            const auto rows = support(row_mask);
            for (unsigned column_mask = 1; column_mask < 16; ++column_mask) {
                if (__builtin_popcount(column_mask) != size) continue;
                unsigned determinant = minor_determinant(matrix, rows, support(column_mask));
                if (determinant) best = std::min(best, __builtin_ctz(determinant));
                if (best == 0) break;
            }
            if (best == 0) break;
        }
        if (best == 99) {
            std::cerr << "Unexpected generic rank deficiency\n";
            std::exit(2);
        }
        orders[size - 1] = best - previous;
        previous = best;
    }
    if (!std::is_sorted(orders.begin(), orders.end()) || orders[0] < 0 || previous > 4)
        std::exit(3);
    return orders;
}

std::string profile_key(const Orders& first, const Orders& second) {
    std::string result;
    for (int order : first) result += std::to_string(order);
    result += ':';
    for (int order : second) result += std::to_string(order);
    return result;
}

struct Record { unsigned count = 0, witness = 0; };

int main() {
    std::map<std::string, Record> catalogue;
    std::vector<unsigned> minimizers;
    for (unsigned mask = 0; mask < 65536; ++mask) {
        Matrix infinity{}, zero{};
        for (int x = 0; x < 2; ++x)
        for (int y = 0; y < 2; ++y)
        for (int u = 0; u < 2; ++u)
        for (int v = 0; v < 2; ++v) {
            const int index = ((2 * x + y) * 2 + u) * 2 + v;
            const unsigned a = mask >> index & 1;
            const unsigned e = (x == v && y == u), b = a ^ e;
            infinity[4 * (2 * x + y) + 2 * u + v] = b | (a << 1);
            zero[4 * (2 * x + u) + 2 * y + v] = a | (b << 1);
        }
        const auto first = smith_orders(infinity), second = smith_orders(zero);
        const auto key = profile_key(first, second);
        auto& record = catalogue[key];
        if (!record.count) record.witness = mask;
        ++record.count;
        if (key == "0112:0112") minimizers.push_back(mask);
    }
    std::cout << "{\"pencils\":65536,\"profiles\":[";
    bool first = true;
    for (const auto& [key, record] : catalogue) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << "{\"profile\":\"" << key << "\",\"count\":" << record.count
                  << ",\"witness\":" << record.witness << '}';
    }
    std::cout << "],\"minimizers\":[";
    for (unsigned i = 0; i < minimizers.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << minimizers[i];
    }
    std::cout << "]}\n";
}
