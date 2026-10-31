// Exact differences of canonical binary columns at packet length four.
// Each row tensor is binary; a column difference has entries -1, 0, 1.
// Three modular projections are necessary filters, followed by every
// original polynomial coefficient. No random exclusion is used.
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
        if (argc != 5) throw std::runtime_error("Expected data, precision, odd/nonzero mode, and seed.");
        constexpr int states = 16, half_states = 6561, buckets = 32768;
        int precision = std::stoi(argv[2]);
        if (precision < 1 || precision > 21) throw std::runtime_error("Precision must be 1 through 21.");
        bool odd_mode = std::string(argv[3]) == "odd";
        if (!odd_mode && std::string(argv[3]) != "nonzero") throw std::runtime_error("Invalid mode.");
        U32 mask = (U32(1) << precision) - 1;
        U64 generator = std::stoull(argv[4]);
        auto random_word = [&]() {
            generator ^= generator << 13;
            generator ^= generator >> 7;
            generator ^= generator << 17;
            return U32(generator);
        };
        std::ifstream input(argv[1]);
        int count;
        if (!(input >> count) || count != 116) throw std::runtime_error("Expected 116 coefficient matrices.");
        using Matrix = std::array<std::array<U32, states>, states>;
        std::vector<Matrix> coefficients(count);
        for (auto& matrix : coefficients) for (auto& row : matrix) for (auto& value : row)
            if (!(input >> value)) throw std::runtime_error("Incomplete coefficient input.");
        std::array<Matrix, 3> projections{};
        for (int projection = 0; projection < 3; ++projection) for (const auto& matrix : coefficients) {
            U32 factor = random_word();
            for (int row = 0; row < states; ++row) for (int column = 0; column < states; ++column)
                projections[projection][row][column] += factor * matrix[row][column];
        }

        std::array<std::array<U32, states>, 3> values{};
        std::array<std::array<U32, half_states>, 3> left_values{}, right_values{};
        std::array<int, half_states> left_sum{}, right_sum{}, next{};
        std::array<int, buckets> heads;
        std::array<U64, buckets> keys{};
        U32 previous = 0;
        U64 candidates = 0, queries = 0;
        for (U32 index = 1; index < 65536; ++index) {
            U32 row_tensor = index ^ (index >> 1);
            U32 change = row_tensor ^ previous;
            int coordinate = __builtin_ctz(change);
            bool added = row_tensor & change;
            previous = row_tensor;
            for (int projection = 0; projection < 3; ++projection) for (int column = 0; column < states; ++column) {
                if (added) values[projection][column] += projections[projection][coordinate][column];
                else values[projection][column] -= projections[projection][coordinate][column];
            }
            for (int projection = 0; projection < 3; ++projection)
                left_values[projection][0] = right_values[projection][0] = 0;
            left_sum[0] = right_sum[0] = 0;
            int size = 1;
            for (int bit = 0; bit < 8; ++bit) {
                for (int old = 0; old < size; ++old) for (int digit = 1; digit <= 2; ++digit) {
                    int current = old + digit * size;
                    int sign = digit == 1 ? 1 : -1;
                    left_sum[current] = left_sum[old] + sign * int((row_tensor >> bit) & 1);
                    right_sum[current] = right_sum[old] + sign * int((row_tensor >> (bit + 8)) & 1);
                    for (int projection = 0; projection < 3; ++projection) {
                        left_values[projection][current] = (left_values[projection][old]
                            + U32(sign) * values[projection][bit]) & mask;
                        right_values[projection][current] = (right_values[projection][old]
                            + U32(sign) * values[projection][bit + 8]) & mask;
                    }
                }
                size *= 3;
            }
            heads.fill(-1);
            auto key = [&](U32 a, U32 b, U32 c) {
                return (U64(a) << (2 * precision)) | (U64(b) << precision) | U64(c);
            };
            auto slot_for = [&](U64 value) {
                int slot = (value * 11400714819323198485ull) >> 49;
                while (heads[slot] >= 0 && keys[slot] != value) slot = (slot + 1) & (buckets - 1);
                return slot;
            };
            for (int left = 0; left < half_states; ++left) {
                U64 value = key(left_values[0][left], left_values[1][left], left_values[2][left]);
                int slot = slot_for(value);
                if (heads[slot] < 0) keys[slot] = value;
                next[left] = heads[slot];
                heads[slot] = left;
            }
            for (int right = 0; right < half_states; ++right) {
                ++queries;
                U64 value = key((-right_values[0][right]) & mask,
                                (-right_values[1][right]) & mask,
                                (-right_values[2][right]) & mask);
                int slot = slot_for(value);
                if (heads[slot] < 0) continue;
                for (int left = heads[slot]; left >= 0; left = next[left]) {
                    int original = left_sum[left] + right_sum[right];
                    if (odd_mode ? original % 2 == 0 : original == 0) continue;
                    ++candidates;
                    std::array<int, states> column{};
                    int left_code = left, right_code = right;
                    for (int bit = 0; bit < 8; ++bit) {
                        int digit = left_code % 3;
                        column[bit] = digit == 2 ? -1 : digit;
                        left_code /= 3;
                        digit = right_code % 3;
                        column[bit + 8] = digit == 2 ? -1 : digit;
                        right_code /= 3;
                    }
                    U32 aggregate = 0;
                    for (const auto& matrix : coefficients) {
                        U32 coefficient = 0;
                        for (U32 rows = row_tensor; rows; rows &= rows - 1) {
                            int row = __builtin_ctz(rows);
                            for (int index = 0; index < states; ++index)
                                coefficient += U32(column[index]) * matrix[row][index];
                        }
                        aggregate |= coefficient;
                        if (aggregate & mask) break;
                    }
                    if (!(aggregate & mask)) {
                        std::cout << "{\"found\":true,\"precision\":" << precision
                                  << ",\"mode\":\"" << argv[3] << "\",\"row_tensor_mask\":" << row_tensor
                                  << ",\"original_difference\":" << original << ",\"column_difference\":[";
                        for (int index = 0; index < states; ++index) {
                            if (index) std::cout << ',';
                            std::cout << column[index];
                        }
                        std::cout << "],\"row_tensors_tested\":" << index
                                  << ",\"right_half_queries\":" << queries
                                  << ",\"projection_survivors\":" << candidates << "}\n";
                        return 0;
                    }
                }
            }
        }
        std::cout << "{\"found\":false,\"precision\":" << precision
                  << ",\"mode\":\"" << argv[3] << "\",\"row_tensors_tested\":65535"
                  << ",\"right_half_queries\":" << queries
                  << ",\"projection_survivors\":" << candidates << "}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
