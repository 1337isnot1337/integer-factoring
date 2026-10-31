// Exact coefficient-space enumeration; polynomial bits are coefficients over F2.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
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

std::pair<Orders, Orders> profile(unsigned mask, int degree = 2,
                                unsigned a = 9, unsigned b = 9) {
    Matrix infinity{}, zero{};
    for (int x = 0; x < 2; ++x) for (int y = 0; y < 2; ++y)
    for (int u = 0; u < 2; ++u) for (int v = 0; v < 2; ++v) {
        int index = 8*x + 4*y + 2*u + v;
        unsigned identity = ((a >> (2*x+u)) & 1) * ((b >> (2*y+v)) & 1);
        unsigned crossing = (x == v && y == u);
        unsigned middle = mask >> index & 1;
        if (degree == 2) {
            unsigned leading = crossing ^ identity ^ middle;
            infinity[2*x+y][2*u+v] = leading | middle << 1 | identity << 2;
            zero[2*x+u][2*y+v] = identity | middle << 1 | leading << 2;
        } else {
            unsigned second = mask >> (16 + index) & 1;
            unsigned leading = crossing ^ identity ^ middle ^ second;
            infinity[2*x+y][2*u+v] = leading | second << 1 | middle << 2 | identity << 3;
            zero[2*x+u][2*y+v] = identity | middle << 1 | second << 2 | leading << 3;
        }
    }
    return {orders(infinity), orders(zero)};
}

std::pair<Orders, Orders> parity_profile(uint64_t packed, int degree) {
    std::array<unsigned, 8> entries{};
    const unsigned identity = 0xa5, crossing = 0x99;
    for (unsigned a = 0; a < 8; ++a) {
        unsigned polynomial = (identity >> a) & 1;
        unsigned endpoint = polynomial;
        for (int j = 1; j < degree; ++j) {
            unsigned coefficient = (packed >> (8*(j-1)+a)) & 1;
            polynomial |= coefficient << j;
            endpoint ^= coefficient;
        }
        polynomial |= (endpoint ^ ((crossing >> a) & 1)) << degree;
        entries[a] = polynomial;
    }
    auto block_orders = [](const std::array<unsigned, 4> &a) {
        int minimum = 100;
        for (unsigned entry : a) if (entry) minimum = std::min(minimum, __builtin_ctz(entry));
        unsigned det = multiply(a[0], a[3]) ^ multiply(a[1], a[2]);
        if (!det) throw std::runtime_error("parity block lost generic rank");
        return std::array<int, 2>{minimum, __builtin_ctz(det)-minimum};
    };
    auto flatten = [&](bool infinity) {
        auto work = entries;
        if (infinity) for (unsigned &entry : work) {
            unsigned reversed = 0;
            for (int j = 0; j <= degree; ++j) reversed |= ((entry >> j) & 1) << (degree-j);
            entry = reversed;
        }
        auto first = infinity ? block_orders({work[0], work[1], work[6], work[7]})
                              : block_orders({work[0], work[2], work[5], work[7]});
        auto second = infinity ? block_orders({work[2], work[3], work[4], work[5]})
                               : block_orders({work[1], work[3], work[4], work[6]});
        Orders result{first[0], first[1], second[0], second[1]};
        std::sort(result.begin(), result.end());
        return result;
    };
    return {flatten(true), flatten(false)};
}

int main(int argc, char **argv) {
    int degree = 2;
    uint64_t count = 65536;
    unsigned constant_a = 9, constant_b = 9;
    bool parity = false;
    bool fast_parity = false;
    if (argc == 3 && std::string(argv[1]) == "--cubic-sample") {
        degree = 3;
        count = std::stoull(argv[2]);
        if (count == 0 || count > (uint64_t(1) << 32)) throw std::runtime_error("invalid sample count");
    } else if (argc == 2 && std::string(argv[1]) == "--parity-cubic") {
        degree = 3;
        parity = true;
    } else if (argc == 3 && std::string(argv[1]) == "--parity-degree") {
        degree = std::stoi(argv[2]);
        if (degree < 2 || degree > 4) throw std::runtime_error("parity degree must be two through four");
        count = uint64_t(1) << (8*(degree-1));
        parity = fast_parity = true;
    } else if (argc == 4 && std::string(argv[1]) == "--constant") {
        constant_a = std::stoul(argv[2]);
        constant_b = std::stoul(argv[3]);
        if (constant_a == 0 || constant_a > 15 || constant_b == 0 || constant_b > 15)
            throw std::runtime_error("constants must be nonzero binary 2x2 matrices");
    } else if (argc != 1) {
        throw std::runtime_error("usage: enumerator [--cubic-sample COUNT | --parity-cubic | --parity-degree DEGREE | --constant A B]");
    }
    std::map<std::pair<Orders, Orders>, std::pair<uint64_t, uint64_t>> catalogue;
    for (uint64_t i = 0; i < count; ++i) {
        uint64_t mask = degree == 2 ? i : uint32_t(i * 0x9e3779b1u + 0x7f4a7c15u);
        if (parity) {
            mask = 0;
            int bit = 0;
            for (unsigned entry = 0; entry < 16; ++entry) {
                if (__builtin_popcount(entry) % 2) continue;
                for (int j = 0; j < degree-1; ++j)
                    mask |= ((i >> (8*j + bit)) & 1) << (16*j + entry);
                ++bit;
            }
        }
        auto key = fast_parity ? parity_profile(i, degree) : profile(mask, degree, constant_a, constant_b);
        auto it = catalogue.find(key);
        if (it == catalogue.end()) catalogue[key] = {1, mask};
        else ++it->second.first;
    }
    std::cout << "{\"field\":2,\"degree\":" << degree << ",\"tensors\":" << count;
    if (parity) std::cout << ",\"mask_schedule\":\"all parity-preserving " << (degree == 3 ? "cubic" : "degree-" + std::to_string(degree)) << " coefficients\"";
    else if (degree == 3) std::cout << ",\"mask_schedule\":\"(i*0x9e3779b1+0x7f4a7c15) mod 2^32\"";
    if (constant_a != 9 || constant_b != 9)
        std::cout << ",\"constant_factors\":[" << constant_a << ',' << constant_b << ']';
    std::cout << ",\"profiles\":[";
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
