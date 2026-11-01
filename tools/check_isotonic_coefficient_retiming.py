#!/usr/bin/env python3
"""Independent exact certificates for monotone coefficient retiming."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import comb, isqrt
from pathlib import Path
import random

from bell_crossing import transform
from isotonic_coefficient_retiming import (
    optimal_blocks, physical_binary_profiles, rational_schedule, retime_binary_polynomial,
    retime_coefficients,
)


def convolution(a, b, prime):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] = (result[i + j] + x * y) % prime
    return result


def literal_two_copy_ranks(rows, prime):
    result = [0] * (2 * len(rows) - 1)
    for left, right in product(zip(*rows), repeat=2):
        for index, value in enumerate(convolution(left, right, prime)):
            result[index] += bool(value)
    return result


def binary_multiply(a, b):
    result = 0
    while b:
        bit = b & -b
        result ^= a << (bit.bit_length() - 1)
        b ^= bit
    return result


def positive_pressure_difference(reference, candidate):
    """Certify candidate pressure > reference pressure on every 0<r<1."""
    base = Counter(a + b for a in reference[0] for b in reference[1])
    changed = Counter(a + b for a in candidate[0] for b in candidate[1])
    minimum, maximum = min(set(base) | set(changed)), max(set(base) | set(changed))
    polynomial = [changed[i] - base[i] for i in range(minimum, maximum + 1)]
    while polynomial[-1] == 0:
        polynomial.pop()
    while polynomial[0] == 0:
        minimum += 1
        polynomial.pop(0)
    multiplicity = 0
    while sum(polynomial) == 0:
        quotient, running = [], 0
        for coefficient in polynomial[:-1]:
            running += coefficient
            quotient.append(running)
        assert running + polynomial[-1] == 0
        polynomial = quotient
        multiplicity += 1
    degree = len(polynomial) - 1
    bernstein = [sum(polynomial[i] * comb(degree - i, k - i)
                     for i in range(k + 1)) for k in range(degree + 1)]
    assert all(value >= 0 for value in bernstein) and any(bernstein)
    rebuilt = [sum(bernstein[k] * comb(degree - k, j - k) * (-1) ** (j - k)
                   for k in range(j + 1)) for j in range(degree + 1)]
    assert rebuilt == polynomial
    return {"difference_shift": minimum, "one_minus_r_multiplicity": multiplicity,
            "quotient_coefficients": polynomial, "quotient_degree": degree,
            "nonnegative_bernstein_coefficients": len(bernstein),
            "strict_comparison": "candidate > reference for every 0 < r < 1"}


def physical_eight_letter_controls():
    from eight_letter_reflection_retraction import (
        displayed_maps, homogeneous_restriction, shift_polynomial,
    )

    source = homogeneous_restriction(*displayed_maps(True), keep_matrices=True)
    q, heights = 8, source["heights"]
    zero = source["matrices"][0]
    values = [0] * q ** 4
    for x, y, u, v in product(range(q), repeat=4):
        values[x + q * y + q * q * u + q ** 3 * v] = shift_polynomial(
            zero[q * x + y][q * u + v], 16 - heights[x] - heights[y])
    assert max(value.bit_length() - 1 for value in values) == 32
    twice = [2 * h for h in heights]
    empty = [0] * q
    zero_heights = (twice, twice, empty, empty)
    infinity_heights = (empty, empty, twice, twice)
    schedule = list(range(0, 65, 2))
    reference = physical_binary_profiles(values, q, schedule, zero_heights,
                                         infinity_heights, (-32, 32))
    assert reference["zero_orders"] == [2 * value for value in source["zero_orders"]]
    assert reference["infinity_orders"] == [2 * value for value in source["infinity_orders"]]
    baseline = (reference["zero_orders"], reference["infinity_orders"])
    candidates = []
    for index in range(1, 16):
        for delta in (-1, 1):
            modified = schedule[:]
            modified[index] += delta
            modified[32 - index] -= delta
            result = physical_binary_profiles(values, q, modified, zero_heights,
                                               infinity_heights, (-32, 32))
            assert result["total_physical_heights"] == reference["total_physical_heights"]
            result["changed_index"] = index
            result["change"] = delta
            result["all_scale_comparison"] = positive_pressure_difference(
                baseline, (result["zero_orders"], result["infinity_orders"]))
            candidates.append(result)
    return {"scope": "all 30 single reciprocal +/-1 changes from d_i=2i, at the doubled natural physical heights",
            "reference": reference, "candidates": candidates,
            "complete_physical_smith_orders": 31 * 2 * 64,
            "literal_physical_endpoint_entries": 31 * q ** 4,
            "all_scale_strictly_dominated_candidates": len(candidates),
            "arbitrary_height_optimization_claimed": False}


def check():
    rng = random.Random(2026110111)
    counts = {"exact_optimality_certificates": 0,
              "retimed_endpoint_coefficients": 0,
              "literal_two_copy_rank_bounds": 0,
              "smith_counterexample_determinants": 0}
    for degree in range(1, 9):
        for _ in range(24):
            low = [rng.randrange(10) for _ in range(degree + 1)]
            high = [rng.randrange(10) for _ in range(degree + 1)]
            low[0] += 1
            high[-1] += 1
            blocks = optimal_blocks(low, high)
            assert blocks[0][0] == 0 and blocks[-1][1] == degree
            for left, right in zip(blocks, blocks[1:]):
                assert left[1] + 1 == right[0]
                assert left[2] * right[3] <= right[2] * left[3]
            for start, stop, a, b in blocks:
                assert sum(low[start:stop + 1]) == a
                assert sum(high[start:stop + 1]) == b
                assert a > 0 and b > 0
                for end in range(start, stop + 1):
                    assert sum(low[start:end + 1]) * b >= a * sum(high[start:end + 1])
            result = rational_schedule(low, high, accuracy=8)
            scale = Fraction(*result["scale"])
            degrees = result["retimed_degrees"]
            assert degrees[0] == 0
            assert all(a < b for a, b in zip(degrees, degrees[1:]))
            low_sum = sum(Fraction(n) / scale ** j for n, j in zip(low, degrees))
            high_sum = sum(Fraction(n) * scale ** j for n, j in zip(high, degrees))
            bound = low_sum * high_sum
            assert bound == Fraction(*result["pressure"])
            denominator = 1 << 40
            lower_root = sum(Fraction(isqrt(a * b * denominator ** 2), denominator)
                             for _, _, a, b in blocks)
            assert bound <= Fraction(*result["approximation_factor"]) * lower_root ** 2
            assert bound >= lower_root ** 2
            counts["exact_optimality_certificates"] += 1

    # Actual coefficient vectors, both transformed cuts, and literal tensor
    # coefficients verify applicability without assuming reciprocity.
    for prime in (3, 5):
        for degree in (2, 3, 4):
            for _ in range(3):
                rows = [[rng.randrange(prime) for _ in range(4)]
                        for _ in range(degree)]
                rows.append([(1 - sum(row[j] for row in rows)) % prime
                             for j in range(4)])
                if not any(rows[0]):
                    rows[0][0] = 1
                    rows[-1][0] = (rows[-1][0] - 1) % prime
                if not any(rows[-1]):
                    rows[0][0] = (rows[0][0] - 1) % prime
                    rows[-1][0] = 1
                assert any(rows[0]) and any(rows[-1])
                other = [transform(tuple(row), 1, prime) for row in rows]
                low_counts = Counter(next(i for i, row in enumerate(rows) if row[j])
                                     for j in range(4))
                high_counts = Counter(max(i for i, row in enumerate(other) if row[j])
                                      for j in range(4))
                low = [low_counts[i] for i in range(degree + 1)]
                high = [high_counts[i] for i in range(degree + 1)]
                result = rational_schedule(low, high, accuracy=4)
                retimed = retime_coefficients(rows, result["retimed_degrees"])
                retimed_other = retime_coefficients(other, result["retimed_degrees"])
                assert tuple(transform(row, 1, prime) for row in retimed) == retimed_other
                for cut in (retimed, retimed_other):
                    assert all(sum(row[j] for row in cut) % prime == 1 for j in range(4))
                    counts["retimed_endpoint_coefficients"] += 4
                row_ranks = literal_two_copy_ranks(retimed, prime)
                col_ranks = literal_two_copy_ranks(retimed_other, prime)
                bound = Fraction(*result["pressure"])
                for row_rank, col_rank in zip(row_ranks, col_ranks):
                    assert row_rank * col_rank <= bound ** 2
                    counts["literal_two_copy_rank_bounds"] += 1

    # A nonconstant schedule, with a pressure strictly below the trivial 16.
    rows = ((0, 0, 0, 2), (1, 1, 1, 2))
    other = tuple(transform(row, 1, 3) for row in rows)
    assert other == ((2, 1, 1, 1), (2, 0, 0, 0))
    explicit = rational_schedule((1, 3), (3, 1), accuracy=4)
    assert explicit["retimed_degrees"] == [0, 6]
    retimed = retime_coefficients(rows, explicit["retimed_degrees"])
    retimed_other = retime_coefficients(other, explicit["retimed_degrees"])
    bound = Fraction(*explicit["pressure"])
    assert bound < 16
    for row_rank, col_rank in zip(literal_two_copy_ranks(retimed, 3),
                                  literal_two_copy_ranks(retimed_other, 3)):
        assert row_rank * col_rank <= bound ** 2
        counts["literal_two_copy_rank_bounds"] += 1

    # Smith orders are not scalar eigenpolynomial orders. This endpoint-I
    # matrix has original orders (0,3), but nonlinear coefficient retiming
    # gives (0,2), not the falsely inferred (0,4).
    original = ((0b0010, 0b0011), (0b1100, 0b1110))
    schedule = (0, 1, 3, 4)
    retimed = tuple(tuple(retime_binary_polynomial(value, schedule) for value in row)
                    for row in original)
    determinants = []
    orders = []
    for matrix in (original, retimed):
        assert [[value.bit_count() % 2 for value in row] for row in matrix] == [[1, 0], [0, 1]]
        assert min((value & -value).bit_length() - 1
                   for row in matrix for value in row if value) == 0
        determinant = (binary_multiply(matrix[0][0], matrix[1][1])
                       ^ binary_multiply(matrix[0][1], matrix[1][0]))
        determinants.append(determinant)
        orders.append([0, (determinant & -determinant).bit_length() - 1])
        counts["smith_counterexample_determinants"] += 1
    assert determinants == [0b1000, 0b11100]
    assert orders == [[0, 3], [0, 2]]
    physical = physical_eight_letter_controls()
    return {
        "status": "passed", "counts": counts,
        "exact_infimum": "(sum over pooled blocks sqrt(A_block B_block))^2",
        "constructive_factor": "at most (1+1/accuracy)^(D+1)",
        "nonconstant_diagonal_example": explicit,
        "smith_retiming_counterexample": {
            "original_matrix_bitmasks": original,
            "retimed_matrix_bitmasks": retimed,
            "degree_schedule": schedule,
            "determinants_bitmasks": determinants,
            "actual_smith_orders": orders,
            "incorrect_profile_only_prediction": [0, 4],
        },
        "applicability": "constant-basis diagonal cuts; general Smith profiles require recomputation",
        "eight_letter_physical_retimings": physical,
        "factoring_score": "1/2; unchanged",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key != "eight_letter_physical_retimings"} | {
        "eight_letter_physical_retimings": {
            key: value for key, value in result["eight_letter_physical_retimings"].items()
            if key not in ("reference", "candidates")}}, indent=2))


if __name__ == "__main__":
    main()
