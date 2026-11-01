// Exhaustive projective F4 rank distance of an explicitly supplied binary space.
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

using Rows = std::array<std::uint32_t, 25>;
std::vector<Rows> basis;
Rows low{}, high{};
std::vector<unsigned> coefficients, witness;
unsigned best = 26;
std::uint64_t checked = 0;

unsigned rank_below(Rows a, Rows b, unsigned cutoff) {
    unsigned rank = 0;
    for (unsigned column = 0; column < 25 && rank < 25; ++column) {
        unsigned row = rank;
        while (row < 25 && !(((a[row] | b[row]) >> column) & 1U)) ++row;
        if (row == 25) continue;
        std::swap(a[row], a[rank]);
        std::swap(b[row], b[rank]);
        unsigned value = ((a[rank] >> column) & 1U) | (((b[rank] >> column) & 1U) << 1);
        if (value == 2) {
            auto old = a[rank];
            a[rank] ^= b[rank]; b[rank] = old;
        } else if (value == 3) {
            auto old = a[rank];
            a[rank] = b[rank]; b[rank] ^= old;
        }
        for (unsigned i = rank + 1; i < 25; ++i) {
            unsigned factor = ((a[i] >> column) & 1U) | (((b[i] >> column) & 1U) << 1);
            if (factor == 1) { a[i] ^= a[rank]; b[i] ^= b[rank]; }
            if (factor == 2) { a[i] ^= b[rank]; b[i] ^= a[rank] ^ b[rank]; }
            if (factor == 3) { a[i] ^= a[rank] ^ b[rank]; b[i] ^= a[rank]; }
        }
        ++rank;
        if (rank >= cutoff) return rank;
    }
    return rank;
}

void toggle(unsigned index, unsigned value) {
    for (unsigned row = 0; row < 25; ++row) {
        if (value & 1U) low[row] ^= basis[index][row];
        if (value & 2U) high[row] ^= basis[index][row];
    }
}

void enumerate(unsigned index) {
    if (index == basis.size()) {
        ++checked;
        unsigned rank = rank_below(low, high, best);
        if (rank < best) { best = rank; witness = coefficients; }
        return;
    }
    for (unsigned value = 0; value < 4; ++value) {
        coefficients[index] = value;
        toggle(index, value);
        enumerate(index + 1);
        toggle(index, value);
    }
    coefficients[index] = 0;
}

int main() {
    unsigned dimension;
    if (!(std::cin >> dimension) || dimension == 0 || dimension > 14) return 1;
    basis.resize(dimension);
    coefficients.assign(dimension, 0);
    for (auto &matrix : basis) for (auto &row : matrix) {
        if (!(std::cin >> row) || row >= (1U << 25)) return 2;
    }
    for (unsigned first = 0; first < dimension; ++first) {
        coefficients[first] = 1;
        toggle(first, 1);
        enumerate(first + 1);
        toggle(first, 1);
        coefficients[first] = 0;
    }
    std::cout << "{\"dimension\":" << dimension << ",\"projective_vectors_checked\":" << checked
              << ",\"minimum_rank\":" << best << ",\"witness_coefficients\":[";
    for (unsigned i = 0; i < witness.size(); ++i) {
        if (i) std::cout << ',';
        std::cout << witness[i];
    }
    std::cout << "]}\n";
}
