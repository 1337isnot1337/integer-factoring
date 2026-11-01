#!/usr/bin/env python3
"""Exact block and gauged-Smith controls for arbitrary common binary codes.

The all-code entropy theorem is proved in the companion note. This replay
uses no floating arithmetic and does not treat finite enumeration as a proof
for larger codes.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
from random import Random

import check_block_smith_crossings as smith
from check_cancellation_module_automata import pmul


ROOT = Path(__file__).resolve().parent.parent


def one_plus_powers(length):
    values = [1]
    for _ in range(length):
        values.append(pmul(values[-1], 3))
    return values


def code_tensor(length, words, second_words=None):
    """Construct literal entries; independently check sum/difference blocks."""
    assert len(set(words)) == len(words) >= 2
    assert all(0 <= word < (1 << length) for word in words)
    other = words if second_words is None else second_words
    assert len(set(other)) == len(other) == len(words)
    assert all(0 <= word < (1 << length) for word in other)
    alphabets = (words, other, other, words)
    powers = one_plus_powers(length)
    q = len(words)
    entries = []
    row_blocks = Counter((x & y, x ^ y) for x, y in product(words, other))
    column_blocks = Counter((x ^ u, x & (x ^ u)) for x, u in product(words, other))
    for xi, yi, ui, vi in product(range(q), repeat=4):
        x, y, u, v = (alphabet[index] for alphabet, index in zip(alphabets, (xi, yi, ui, vi)))
        value = 1
        for bit in range(length):
            xb, yb, ub, vb = ((word >> bit) & 1 for word in (x, y, u, v))
            straight = int(xb == ub and yb == vb)
            crossing = int(xb == vb and yb == ub)
            local = straight | ((straight ^ crossing) << 1)
            value = pmul(value, local)
        assert value.bit_count() % 2 == int(x == v and y == u)

        difference = x ^ y
        if (x & y, difference) == (u & v, u ^ v):
            width = difference.bit_count()
            distance = (x ^ u).bit_count()
            row_prediction = powers[width-distance] << distance
            assert row_prediction.bit_length()-1 == width
        else:
            row_prediction = 0
        assert value == row_prediction

        difference = x ^ u
        if y ^ v == difference and v & difference == x & difference:
            width = difference.bit_count()
            distance = ((x ^ y) & ~difference).bit_count()
            column_prediction = powers[distance] << width
            assert column_prediction & -column_prediction == 1 << width
        else:
            column_prediction = 0
        assert value == column_prediction
        entries.append(((xi, yi, ui, vi), value))

    degree = max((x ^ y).bit_count() for x, y in product(words, other))
    assert max(value.bit_length()-1 for _, value in entries) == degree
    assert sum(row_blocks.values()) == sum(column_blocks.values()) == q*q
    return entries, row_blocks, column_blocks, degree


def endpoint_profile(entries, q, heights, cut, infinity):
    """Return actual Laurent orders, restoring the computational shift."""
    terms = []
    minimum = maximum = None
    for indices, polynomial in entries:
        height = sum(heights[leg][index] for leg, index in enumerate(indices))
        exponents = [height-k if infinity else height+k
                     for k in range(polynomial.bit_length()) if polynomial >> k & 1]
        if exponents:
            minimum = min(exponents) if minimum is None else min(minimum, min(exponents))
            maximum = max(exponents) if maximum is None else max(maximum, max(exponents))
        terms.append((indices, exponents))
    assert minimum is not None and maximum is not None
    matrix = [[0]*(q*q) for _ in range(q*q)]
    for (x, y, u, v), exponents in terms:
        row, column = (q*x+y, q*u+v) if cut == 0 else (q*x+u, q*y+v)
        matrix[row][column] = sum(1 << (exponent-minimum) for exponent in exponents)
    orders, _ = smith.smith_profile(matrix, maximum-minimum)
    assert len(orders) == q*q
    return tuple(order+minimum for order in orders)


def pressure(orders, scale):
    return sum((scale ** -order for order in orders), Fraction(0))


def crossing_block_orders(words, zero_heights, infinity_heights, second_words=None, opposite=False):
    row, column = {}, {}
    other = words if second_words is None else second_words
    row_heights = zero_heights if opposite else infinity_heights
    column_heights = infinity_heights if opposite else zero_heights
    sign = 1 if opposite else -1
    for x, y in product(range(len(words)), repeat=2):
        distance = (words[x] ^ other[y]).bit_count()
        indices = (x, y, y, x)
        row_order = sum(row_heights[leg][index] for leg, index in enumerate(indices)) + sign*distance
        column_order = sum(column_heights[leg][index] for leg, index in enumerate(indices)) - sign*distance
        row_key = (words[x] & other[y], words[x] ^ other[y])
        column_key = (words[x] ^ other[y], words[x] & (words[x] ^ other[y]))
        row[row_key] = min(row.get(row_key, row_order), row_order)
        column[column_key] = min(column.get(column_key, column_order), column_order)
    return tuple(row.values()), tuple(column.values())


def gauge_control(length, words, seed, second_words=None):
    entries, _, _, degree = code_tensor(length, words, second_words)
    q = len(words)
    random = Random(seed)
    heights = [tuple(tuple(random.randrange(-2, 3) for _ in words) for _ in range(4))
               for _ in range(2)]
    zero, infinity = heights
    profiles = {}
    for endpoint, height, reverse in (("zero", zero, False), ("infinity", infinity, True)):
        for cut in (0, 1):
            profiles[endpoint, cut] = endpoint_profile(entries, q, height, cut, reverse)
    row_block, column_block = crossing_block_orders(words, zero, infinity, second_words)
    opposite_row_block, opposite_column_block = crossing_block_orders(
        words, zero, infinity, second_words, opposite=True)
    complement_aligned = second_words is not None and all(
        left ^ right == (1 << length)-1 for left, right in zip(words, second_words))
    energy = [sum(height[leg][index] for height in heights for leg in range(4))
              for index in range(q)]
    comparisons = []
    for root in (Fraction(1), Fraction(6, 5), Fraction(3, 2)):
        # Twelfth powers make the rho=3/4 and rho=1 thresholds rational
        # for arbitrary integer energy, without large gauge exponents.
        scale = root ** 12
        good = pressure(profiles["infinity", 0], scale) * pressure(profiles["zero", 1], scale)
        block = pressure(row_block, scale) * pressure(column_block, scale)
        threshold_three = sum((root ** (-4*value) for value in energy), Fraction(0)) ** 3
        opposite = pressure(profiles["zero", 0], scale) * pressure(profiles["infinity", 1], scale)
        opposite_block = pressure(opposite_row_block, scale) * pressure(opposite_column_block, scale)
        threshold_four = sum((root ** (-3*value) for value in energy), Fraction(0)) ** 4
        assert good >= block >= threshold_three
        assert opposite >= opposite_block >= threshold_three

        # Independently check the four partial-charge factors in the
        # principal-minor majorization argument for the opposite cut.
        if second_words is None or complement_aligned:
            if complement_aligned:
                zero_pairs, infinity_pairs, shift = ((0, 3), (1, 2)), ((0, 1), (2, 3)), length
            else:
                zero_pairs, infinity_pairs, shift = ((0, 2), (1, 3)), ((0, 3), (1, 2)), 0
            opposite_zero = scale ** -shift
            opposite_infinity = scale ** shift
            for first, second in zero_pairs:
                opposite_zero *= sum((scale ** (-zero[first][i]-zero[second][i])
                                      for i in range(q)), Fraction(0))
            for first, second in infinity_pairs:
                opposite_infinity *= sum((scale ** (-infinity[first][i]-infinity[second][i])
                                          for i in range(q)), Fraction(0))
            assert pressure(profiles["zero", 0], scale) >= opposite_zero
            assert pressure(profiles["infinity", 1], scale) >= opposite_infinity
            assert opposite_zero * opposite_infinity >= threshold_four
        comparisons.append({
            "scale": str(scale),
            "good_pressure": str(good),
            "one_factor_per_block_floor": str(block),
            "rho_three_quarters_threshold": str(threshold_three),
            "opposite_pressure": str(opposite),
            "opposite_one_factor_per_block_floor": str(opposite_block),
            "rho_one_threshold_when_applicable": str(threshold_four)
                if second_words is None or complement_aligned else None,
        })
    return {
        "length": length, "words": words, "alphabet": q, "ungauged_degree": degree,
        "paired_second_code": second_words,
        "complement_aligned": complement_aligned,
        "zero_gauges": zero, "infinity_gauges": infinity, "symbol_energy": energy,
        "profiles": {endpoint+"_cut_"+str(cut): sorted(Counter(orders).items())
                     for (endpoint, cut), orders in profiles.items()},
        "comparisons": comparisons,
    }


def run():
    families = tensor_entries = sum_blocks = difference_blocks = 0
    for length in range(1, 4):
        universe = 1 << length
        for selected in range(1, 1 << universe):
            if selected.bit_count() < 2:
                continue
            words = tuple(word for word in range(universe) if selected >> word & 1)
            entries, row, column, _ = code_tensor(length, words)
            families += 1
            tensor_entries += len(entries)
            sum_blocks += len(row)
            difference_blocks += len(column)
    assert families == 259 and tensor_entries == 118188

    cases = (
        (2, (0, 1, 2)),
        (3, (0, 1, 3, 6)),
        (3, (0, 2, 3, 5, 7)),
        (4, (0, 1, 6, 10, 15)),
        (4, (1, 3, 4, 8, 10, 14)),
        (5, (0, 3, 7, 18, 23, 29)),
    )
    controls = [gauge_control(length, words, seed)
                for seed, (length, words) in enumerate(cases, 8130)]
    paired_cases = (
        (3, (0, 1, 3, 6), (2, 4, 5, 7)),
        (4, (0, 1, 6, 10, 15), (2, 3, 7, 9, 14)),
        (4, (1, 3, 4, 8, 10), (14, 12, 11, 7, 5)),
        (5, (0, 3, 7, 18, 23, 29), (8, 30, 12, 3, 21, 1)),
    )
    paired_controls = [gauge_control(length, words, seed, second)
                       for seed, (length, words, second) in enumerate(paired_cases, 9210)]
    return {
        "scope": "arbitrary common binary codes with diagonal physical monomial gauges",
        "proved_good_orientation_rho_floor": "3/4",
        "proved_opposite_orientation_rho_floor": "1",
        "exhaustive_codes_through_length_three": families,
        "literal_tensor_entries_checked": tensor_entries,
        "sum_blocks_checked": sum_blocks,
        "difference_blocks_checked": difference_blocks,
        "nonuniform_gauged_controls": controls,
        "paired_code_controls": paired_controls,
        "complete_laurent_profiles_computed": 4*(len(controls)+len(paired_controls)),
        "exact_pressure_comparison_scales": 3*(len(controls)+len(paired_controls)),
        "entropy_theorem_proved_analytically_in_note": True,
        "floating_arithmetic_used": False,
        "factoring_exponent_improved": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "evidence" / "arbitrary-binary-code-pressure.json")
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
