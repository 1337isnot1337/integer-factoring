#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

// Complete two-generator regroupings around the four stated pairing targets.
// F4 is F2[w]/(w^2+w+1); both cuts have size 49.
using BinaryMatrix = std::array<std::uint64_t, 49>;
struct BinaryTensor { BinaryMatrix row, column; };
struct FieldRow { std::uint64_t low = 0, high = 0; };
using FieldMatrix = std::array<FieldRow, 49>;
struct FieldTensor { FieldMatrix row, column; };

FieldRow scale(FieldRow value, int scalar) {
    if (scalar == 0) return {};
    if (scalar == 1) return value;
    if (scalar == 2) return {value.high, value.low ^ value.high};
    return {value.low ^ value.high, value.low};
}

int rank(FieldMatrix matrix, int cap = 49) {
    std::array<FieldRow, 49> pivots{};
    int result = 0;
    for (auto value : matrix) {
        while (value.low | value.high) {
            const int pivot = 63 - std::countl_zero(value.low | value.high);
            const int scalar = ((value.low >> pivot) & 1) | (((value.high >> pivot) & 1) << 1);
            if (!(pivots[pivot].low | pivots[pivot].high)) {
                pivots[pivot] = scale(value, scalar == 1 ? 1 : 5 - scalar);
                if (++result > cap) return result;
                break;
            }
            const auto old = scale(pivots[pivot], scalar);
            value.low ^= old.low;
            value.high ^= old.high;
        }
    }
    return result;
}

FieldTensor combine(const BinaryTensor& first, const BinaryTensor& second, int a, int b) {
    FieldTensor result;
    for (int i = 0; i < 49; ++i) {
        result.row[i] = {((a & 1) ? first.row[i] : 0) ^ ((b & 1) ? second.row[i] : 0),
                         ((a & 2) ? first.row[i] : 0) ^ ((b & 2) ? second.row[i] : 0)};
        result.column[i] = {((a & 1) ? first.column[i] : 0) ^ ((b & 1) ? second.column[i] : 0),
                            ((a & 2) ? first.column[i] : 0) ^ ((b & 2) ? second.column[i] : 0)};
    }
    return result;
}

FieldTensor add(FieldTensor first, const BinaryTensor& second) {
    for (int i = 0; i < 49; ++i) {
        first.row[i].low ^= second.row[i];
        first.column[i].low ^= second.column[i];
    }
    return first;
}

bool improves_36(int first_product, int second_product, int fixed_cost) {
    const long long remaining = 36 - fixed_cost;
    const long long difference = remaining * remaining - first_product - second_product;
    return difference > 0 && 4LL * first_product * second_product < difference * difference;
}

int main(int argc, char** argv) {
    if (argc != 2 || (std::string(argv[1]) != "2" && std::string(argv[1]) != "4")) return 2;
    const int field = std::stoi(argv[1]);
    int count;
    if (!(std::cin >> count) || count < 1 || count > 10000) return 2;
    std::array<BinaryTensor, 4> targets;
    const std::array<int, 4> fixed_cost{0, 7, 7, 14};
    for (auto& tensor : targets)
        for (auto* matrix : {&tensor.row, &tensor.column})
            for (auto& value : *matrix) if (!(std::cin >> std::hex >> value)) return 2;
    std::vector<BinaryTensor> candidates(count);
    for (auto& tensor : candidates)
        for (auto* matrix : {&tensor.row, &tensor.column})
            for (auto& value : *matrix) if (!(std::cin >> std::hex >> value)) return 2;
    std::vector<std::array<int, 2>> representatives{{1, 1}};
    if (field == 4) representatives = {{1, 1}, {1, 2}, {2, 1}, {2, 2}, {2, 3}};
    std::uint64_t cases = 0, comparisons = 0, improvements = 0;
    double best = 36;
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < i; ++j) {
            for (const auto scalar : representatives) {
                ++cases;
                const auto tensor = combine(candidates[i], candidates[j], scalar[0], scalar[1]);
                const int row_rank = rank(tensor.row), column_rank = rank(tensor.column);
                const int first_product = row_rank * column_rank;
                for (int target = 0; target < 4; ++target) {
                    const int remaining = 36 - fixed_cost[target];
                    if (first_product >= remaining * remaining) continue;
                    ++comparisons;
                    const auto residual = add(tensor, targets[target]);
                    const int residual_row_rank = rank(residual.row);
                    int cap = 0;
                    while (cap < 49 && improves_36(first_product, residual_row_rank * (cap + 1),
                                                  fixed_cost[target])) ++cap;
                    if (!cap && residual_row_rank) continue;
                    const int residual_column_rank = rank(residual.column, cap);
                    const int second_product = residual_row_rank * residual_column_rank;
                    if (!improves_36(first_product, second_product, fixed_cost[target])) continue;
                    ++improvements;
                    const double cost = fixed_cost[target] + std::sqrt(first_product) + std::sqrt(second_product);
                    if (cost >= best) continue;
                    best = cost;
                    std::cout << "{\"event\":\"improvement\",\"indices\":[" << i << ',' << j
                              << "],\"scalars\":[" << scalar[0] << ',' << scalar[1]
                              << "],\"target\":" << target << ",\"rank_pairs\":[["
                              << row_rank << ',' << column_rank << "],[" << residual_row_rank << ','
                              << residual_column_rank << "]],\"approximate_cost\":" << cost << "}\n";
                }
            }
        }
    }
    std::cout << "{\"event\":\"complete\",\"field\":" << field << ",\"candidates\":" << count
              << ",\"coefficient_cases\":" << cases << ",\"exact_target_comparisons\":" << comparisons
              << ",\"improvements_below_36\":" << improvements << "}\n";
}
