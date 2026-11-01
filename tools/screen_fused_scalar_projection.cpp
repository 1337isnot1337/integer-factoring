// Exhaustive binary coefficient screens; see fused-scalar-projection-envelope.md.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

constexpr int THRESHOLD = 501;

int row_rank(const uint16_t* rows, int n) {
    uint16_t pivots[9] = {};
    int rank = 0;
    for (int i = 0; i < n; ++i) {
        uint16_t value = rows[i];
        while (value) {
            int p = __builtin_ctz(value);
            if (pivots[p]) value ^= pivots[p];
            else { pivots[p] = value; ++rank; break; }
        }
    }
    return rank;
}

// A separate algorithm for the dual affine hyperplane replay.
int column_rank(const uint16_t* rows, int n) {
    uint16_t matrix[9];
    std::copy(rows, rows+n, matrix);
    int rank = 0;
    for (int c = 0; c < n; ++c) {
        int p = rank;
        while (p < n && !(matrix[p] & (1 << c))) ++p;
        if (p == n) continue;
        std::swap(matrix[p], matrix[rank]);
        for (int i = rank+1; i < n; ++i)
            if (matrix[i] & (1 << c)) matrix[i] ^= matrix[rank];
        ++rank;
    }
    return rank;
}

int main(int argc, char** argv) {
    if (argc != 4 && argc != 6) return 1;
    const std::string mode = argv[3];
    if (mode != "full" && mode != "dual") return 1;
    std::ifstream input(argv[1], std::ios::binary);
    uint32_t dimension, blocks;
    uint64_t target, dual;
    input.read(reinterpret_cast<char*>(&dimension), 4);
    input.read(reinterpret_cast<char*>(&blocks), 4);
    input.read(reinterpret_cast<char*>(&target), 8);
    input.read(reinterpret_cast<char*>(&dual), 8);
    if (!input || dimension > 40 || blocks > 81 || !dual) return 2;
    std::vector<uint32_t> sizes(blocks);
    input.read(reinterpret_cast<char*>(sizes.data()), 4*blocks);
    int length = 0;
    for (auto size : sizes) {
        if (!size || size > 9) return 3;
        length += size;
    }
    std::vector<uint16_t> basis(2*dimension*length);
    input.read(reinterpret_cast<char*>(basis.data()), 2*basis.size());
    if (!input) return 4;

    std::vector<uint64_t> generators;
    uint64_t offset = 0;
    if (mode == "dual") {
        int p = __builtin_ctzll(dual);
        offset = 1ULL << p;
        for (unsigned i = 0; i < dimension; ++i)
            if (int(i) != p)
                generators.push_back((1ULL << i) ^ (((dual >> i) & 1) << p));
    } else {
        for (unsigned i = 0; i < dimension; ++i) generators.push_back(1ULL << i);
    }
    if (argc == 6) {
        unsigned bits = std::stoul(argv[4]);
        uint64_t shard = std::stoull(argv[5]);
        if (bits > generators.size() || shard >= (1ULL << bits)) return 5;
        for (unsigned i = 0; i < bits; ++i)
            if ((shard >> i) & 1) offset ^= generators[generators.size()-bits+i];
        generators.resize(generators.size()-bits);
    }
    std::vector<uint16_t> row(length), column(length);
    auto toggle = [&](uint64_t coordinate) {
        while (coordinate) {
            int bit = __builtin_ctzll(coordinate);
            coordinate &= coordinate-1;
            for (int j = 0; j < length; ++j) {
                row[j] ^= basis[bit*length+j];
                column[j] ^= basis[(dimension+bit)*length+j];
            }
        }
    };
    toggle(offset);
    std::ofstream output(argv[2]);
    if (!output) return 6;
    output << "coordinates,row_rank,column_rank\n";
    uint64_t cheap = 0, previous = 0, coordinate = offset;
    const uint64_t states = 1ULL << generators.size();
    for (uint64_t k = 0; k < states; ++k) {
        uint64_t gray = k ^ (k >> 1);
        if (k) {
            int bit = __builtin_ctzll(gray ^ previous);
            coordinate ^= generators[bit];
            toggle(generators[bit]);
        }
        previous = gray;
        if (!coordinate) continue;
        int rr = 0, rc = 0, position = 0;
        bool rejected = false;
        for (auto size : sizes) {
            auto rank = mode == "dual" ? column_rank : row_rank;
            rr += rank(row.data()+position, size);
            rc += rank(column.data()+position, size);
            position += size;
            if (rr*rc > THRESHOLD) { rejected = true; break; }
        }
        if (!rejected) {
            ++cheap;
            output << coordinate << ',' << rr << ',' << rc << '\n';
        }
    }
    output.close();
    if (!output) return 7;
    std::cout << "{\"dimension\":" << dimension << ",\"states\":" << states
              << ",\"cheap\":" << cheap << ",\"target\":" << target << "}\n";
}
