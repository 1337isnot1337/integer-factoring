// Exact finite binary-tensor visibility screen at packet length four.
// The two projections only filter candidates. Every surviving tensor pair
// is checked against every supplied polynomial coefficient.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using U32 = std::uint32_t;
using U64 = std::uint64_t;

int main(int argc, char** argv) {
    try {
        if (argc != 4) throw std::runtime_error("Expected data path, precision, and seed.");
        constexpr int states = 16;
        constexpr int buckets = 1024;
        int precision = std::stoi(argv[2]);
        if (precision < 1 || precision > 31) throw std::runtime_error("Precision must be 1 through 31.");
        U32 mask = (U32(1) << precision) - 1;
        U64 generator = std::stoull(argv[3]);
        auto random_word = [&]() {
            generator ^= generator << 13;
            generator ^= generator >> 7;
            generator ^= generator << 17;
            return U32(generator);
        };
        std::ifstream input(argv[1]);
        int count;
        if (!(input >> count) || count != 116) throw std::runtime_error("Expected all 116 local coefficients.");
        using Matrix = std::array<std::array<U32, states>, states>;
        std::vector<Matrix> coefficients(count);
        for (auto& matrix : coefficients) for (auto& row : matrix) for (auto& value : row)
            if (!(input >> value)) throw std::runtime_error("Incomplete coefficient input.");
        Matrix first{}, second{};
        for (const auto& matrix : coefficients) {
            U32 a = random_word(), b = random_word();
            for (int row = 0; row < states; ++row) for (int column = 0; column < states; ++column) {
                first[row][column] += a * matrix[row][column];
                second[row][column] += b * matrix[row][column];
            }
        }

        std::array<U32, states> first_column{}, second_column{};
        U32 previous = 0;
        U64 candidates = 0, queries = 0;
        for (U32 index = 1; index < 65536; ++index) {
            U32 column_tensor = index ^ (index >> 1);
            U32 changed = column_tensor ^ previous;
            int coordinate = __builtin_ctz(changed);
            bool added = column_tensor & changed;
            for (int row = 0; row < states; ++row) {
                if (added) {
                    first_column[row] += first[row][coordinate];
                    second_column[row] += second[row][coordinate];
                } else {
                    first_column[row] -= first[row][coordinate];
                    second_column[row] -= second[row][coordinate];
                }
            }
            previous = column_tensor;

            std::array<U32, 256> left_first{}, left_second{}, right_first{}, right_second{};
            std::array<int, 256> left_parity{}, right_parity{}, next{};
            std::array<int, buckets> heads;
            std::array<U64, buckets> keys{};
            heads.fill(-1);
            for (int subset = 1; subset < 256; ++subset) {
                int bit = __builtin_ctz(unsigned(subset));
                int rest = subset & (subset - 1);
                left_first[subset] = (left_first[rest] + first_column[bit]) & mask;
                left_second[subset] = (left_second[rest] + second_column[bit]) & mask;
                left_parity[subset] = left_parity[rest] ^ ((column_tensor >> bit) & 1);
                right_first[subset] = (right_first[rest] + first_column[bit + 8]) & mask;
                right_second[subset] = (right_second[rest] + second_column[bit + 8]) & mask;
                right_parity[subset] = right_parity[rest] ^ ((column_tensor >> (bit + 8)) & 1);
            }
            auto key = [&](U32 a, U32 b, int parity) {
                return (U64(a) << (precision + 1)) | (U64(b) << 1) | U64(parity);
            };
            auto slot_for = [&](U64 value) {
                int slot = (value * 11400714819323198485ull) >> 54;
                while (heads[slot] >= 0 && keys[slot] != value) slot = (slot + 1) & (buckets - 1);
                return slot;
            };
            for (int subset = 0; subset < 256; ++subset) {
                U64 value = key(left_first[subset], left_second[subset], left_parity[subset]);
                int slot = slot_for(value);
                if (heads[slot] < 0) keys[slot] = value;
                next[subset] = heads[slot];
                heads[slot] = subset;
            }
            for (int right = 0; right < 256; ++right) {
                ++queries;
                U64 value = key((-right_first[right]) & mask, (-right_second[right]) & mask,
                                1 ^ right_parity[right]);
                int slot = slot_for(value);
                if (heads[slot] < 0) continue;
                for (int left = heads[slot]; left >= 0; left = next[left]) {
                    ++candidates;
                    U32 row_tensor = left | (right << 8);
                    U32 aggregate = 0;
                    for (const auto& matrix : coefficients) {
                        U32 coefficient = 0;
                        for (U32 rows = row_tensor; rows; rows &= rows - 1) {
                            int row = __builtin_ctz(rows);
                            for (U32 columns = column_tensor; columns; columns &= columns - 1) {
                                int column = __builtin_ctz(columns);
                                coefficient += matrix[row][column];
                            }
                        }
                        aggregate |= coefficient;
                        if (aggregate & mask) break;
                    }
                    if (!(aggregate & mask)) {
                        std::cout << "{\"found\":true,\"precision\":" << precision
                                  << ",\"row_tensor_mask\":" << row_tensor
                                  << ",\"column_tensor_mask\":" << column_tensor
                                  << ",\"column_tensors_tested\":" << index
                                  << ",\"right_subset_queries\":" << queries
                                  << ",\"projection_survivors\":" << candidates << "}\n";
                        return 0;
                    }
                }
            }
        }
        std::cout << "{\"found\":false,\"precision\":" << precision
                  << ",\"column_tensors_tested\":65535,\"odd_original_pairs_covered\":2147450880"
                  << ",\"right_subset_queries\":" << queries
                  << ",\"projection_survivors\":" << candidates << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
