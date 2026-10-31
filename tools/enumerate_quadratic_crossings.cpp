// Exact coefficient-space enumeration; polynomial bits are coefficients over F2.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

using Matrix = std::array<std::array<unsigned, 4>, 4>;
using Orders = std::array<int, 4>;

unsigned multiply(unsigned a, unsigned b) {
    unsigned result = 0;
    while (b) {
        if (b & 1) result ^= a;
        a <<= 1;
        b >>= 1;
    }
    return result;
}

Orders orders(const Matrix &a) {
    Orders result{};
    int previous = 0;
    for (int k = 1; k <= 4; ++k) {
        int minimum = 100;
        for (unsigned rows = 1; rows < 16; ++rows) {
            if (__builtin_popcount(rows) != k) continue;
            std::vector<int> r;
            for (int i = 0; i < 4; ++i) if (rows >> i & 1) r.push_back(i);
            for (unsigned columns = 1; columns < 16; ++columns) {
                if (__builtin_popcount(columns) != k) continue;
                std::vector<int> c;
                for (int i = 0; i < 4; ++i) if (columns >> i & 1) c.push_back(i);
                unsigned determinant = 0;
                do {
                    unsigned term = 1;
                    for (int i = 0; i < k; ++i) term = multiply(term, a[r[i]][c[i]]);
                    determinant ^= term;
                } while (std::next_permutation(c.begin(), c.end()));
                if (determinant) minimum = std::min(minimum, __builtin_ctz(determinant));
            }
        }
        if (minimum == 100) throw std::runtime_error("generic rank lost");
        result[k - 1] = minimum - previous;
        previous = minimum;
    }
    return result;
}

std::pair<Orders, Orders> profile(unsigned mask) {
    Matrix infinity{}, zero{};
    for (int x = 0; x < 2; ++x) for (int y = 0; y < 2; ++y)
    for (int u = 0; u < 2; ++u) for (int v = 0; v < 2; ++v) {
        int index = 8*x + 4*y + 2*u + v;
        unsigned identity = (x == u && y == v);
        unsigned crossing = (x == v && y == u);
        unsigned middle = mask >> index & 1;
        unsigned leading = crossing ^ identity ^ middle;
        infinity[2*x+y][2*u+v] = leading | middle << 1 | identity << 2;
        zero[2*x+u][2*y+v] = identity | middle << 1 | leading << 2;
    }
    return {orders(infinity), orders(zero)};
}

int main() {
    std::map<std::pair<Orders, Orders>, std::pair<unsigned, unsigned>> catalogue;
    for (unsigned mask = 0; mask < 65536; ++mask) {
        auto key = profile(mask);
        auto it = catalogue.find(key);
        if (it == catalogue.end()) catalogue[key] = {1, mask};
        else ++it->second.first;
    }
    std::cout << "{\"field\":2,\"degree\":2,\"tensors\":65536,\"profiles\":[";
    bool first = true;
    for (const auto &entry : catalogue) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << "{\"infinity\":[";
        for (int i = 0; i < 4; ++i) std::cout << (i ? "," : "") << entry.first.first[i];
        std::cout << "],\"zero\":[";
        for (int i = 0; i < 4; ++i) std::cout << (i ? "," : "") << entry.first.second[i];
        std::cout << "],\"count\":" << entry.second.first
                  << ",\"first_mask\":" << entry.second.second << '}';
    }
    std::cout << "]}\n";
}
