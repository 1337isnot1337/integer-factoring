#!/usr/bin/env python3
"""Audit factorial-baseline normalization and bilinear state-tree width."""

from __future__ import annotations

import argparse
import json
import random
from fractions import Fraction
from math import factorial, isqrt
from pathlib import Path

from circuit_box_products import factorial_valuation


def valuation_two(value):
    if not value:
        raise ValueError("Zero has infinite valuation")
    value = abs(value)
    return (value & -value).bit_length() - 1


def normalized_capped_row_valuation(number, first_b, count, row, cap):
    """Known factorial baseline plus inverse-window carry bits; no row list."""
    if min(number, first_b, count, row, cap) < 1 or any(x % 2 == 0 for x in (number, first_b, row)):
        raise ValueError("Require positive odd number, row and initial column")
    if cap < count.bit_length():
        raise ValueError("Cap must include the complete factorial baseline")
    baseline = count + factorial_valuation(count)
    numerator = (number - row * first_b) // 2
    bits = []
    for precision in range(1, cap):
        modulus = 1 << precision
        root_residue = numerator * pow(row, -1, modulus) % modulus
        tail = count % modulus
        bits.append(int(root_residue < tail))
    return {"baseline": baseline, "residual_valuation": sum(bits),
            "capped_total_valuation": baseline + sum(bits), "carry_bits": bits}


def literal_capped_row(number, first_b, count, row, cap):
    return sum(min(cap, valuation_two(value)) if value else cap
               for value in (number - row * (first_b + 2 * index) for index in range(count)))


def is_prime(number):
    return number >= 2 and all(number % divisor for divisor in range(2, isqrt(number) + 1))


def next_prime(number):
    while not is_prime(number):
        number += 1
    return number


def residue_extrema(lower, upper, residue, modulus):
    first = lower + (residue - lower) % modulus
    last = upper - 1 - (upper - 1 - residue) % modulus
    assert first <= last
    return first, last


def controls():
    generator = random.Random(20261101)
    valuation_checks = 0
    baseline_identity_checks = 0
    residual_bound_checks = 0
    slope_checks = 0
    rectangle_checks = 0
    for _ in range(400):
        number = 2 * generator.randrange(1, 100) + 1
        row = 2 * generator.randrange(0, 15) + 1
        first_b = 2 * generator.randrange(0, 15) + 1
        count = generator.randrange(1, 13)
        differences = [number - row * (first_b + 2 * j) for j in range(count)]
        height = max(1, max(abs(value) for value in differences).bit_length())
        start_cap = max(count.bit_length(), 2)
        for cap in range(start_cap, max(start_cap, height + 2) + 1):
            computed = normalized_capped_row_valuation(number, first_b, count, row, cap)
            assert computed["capped_total_valuation"] == literal_capped_row(number, first_b, count, row, cap)
            valuation_checks += 1
        z = Fraction(number - row * first_b, 2 * row)
        binomial = Fraction(1)
        product = 1
        for j in range(count):
            binomial *= (z - j) / (j + 1)
            product *= differences[j]
        assert product == (1 << count) * row ** count * factorial(count) * binomial
        baseline_identity_checks += 1
        if product:
            residual = valuation_two(product) - count - factorial_valuation(count)
            assert 0 <= residual < height - 1
            residual_bound_checks += 1
        cap = max(count.bit_length(), height - 1)
        current = normalized_capped_row_valuation(number, first_b, count, row, cap)
        following = normalized_capped_row_valuation(number, first_b, count, row, cap + 1)
        assert following["residual_valuation"] - current["residual_valuation"] == differences.count(0)
        slope_checks += 1
        first_a = 2 * generator.randrange(0, 12) + 1
        rows = generator.randrange(1, 10)
        largest = max(abs(number - a * b)
                      for a in (first_a, first_a + 2 * (rows - 1))
                      for b in (first_b, first_b + 2 * (count - 1)))
        cap = max(count.bit_length(), max(1, largest.bit_length()) - 1)
        slope = 0
        exact = 0
        for i in range(rows):
            a = first_a + 2 * i
            current = normalized_capped_row_valuation(number, first_b, count, a, cap)
            following = normalized_capped_row_valuation(number, first_b, count, a, cap + 1)
            slope += following["residual_valuation"] - current["residual_valuation"]
            exact += sum(a * (first_b + 2 * j) == number for j in range(count))
        assert slope == exact
        rectangle_checks += 1

    count = (1 << 120) + 17
    first_b, row, root_index = 5, 7, (1 << 119) + 9
    number = row * (first_b + 2 * root_index)
    huge_rows = []
    for offset in (0, 2):
        target = number + offset
        maximum = max(abs(target - row * first_b), abs(target - row * (first_b + 2 * (count - 1))))
        cap = max(count.bit_length(), maximum.bit_length() - 1)
        current = normalized_capped_row_valuation(target, first_b, count, row, cap)
        following = normalized_capped_row_valuation(target, first_b, count, row, cap + 1)
        actual = int(offset == 0)
        assert following["residual_valuation"] - current["residual_valuation"] == actual
        huge_rows.append({"number": target, "row": row, "first_b": first_b,
                          "count": count, "cap": cap,
                          "baseline": current["baseline"],
                          "residual_at_cap": current["residual_valuation"],
                          "residual_at_next_cap": following["residual_valuation"],
                          "factor_count": actual})

    width_controls = []
    literal_states = 0
    for exponent in range(4, 16):
        size = 1 << exponent
        first_prime = next_prime(size // 2 + 1)
        second_prime = next_prime(8 * size + 1)
        assert size // 2 < first_prime < size and 8 * size < second_prime < 16 * size
        number = first_prime * second_prime
        assert 4 * size * size < number < 16 * size * size
        intervals = (size, 3 * size), (3 * size, 9 * size)
        assert all(not intervals[0][0] <= divisor < intervals[0][1]
                   for divisor in (1, first_prime, second_prime, number))
        verified = 0
        for depth in range(1, exponent - 2):
            modulus = 1 << depth
            for first_residue in range(1, modulus, 2):
                second_residue = number * pow(first_residue, -1, modulus) % modulus
                low_a, high_a = residue_extrema(*intervals[0], first_residue, modulus)
                low_b, high_b = residue_extrema(*intervals[1], second_residue, modulus)
                assert low_a < high_a and low_b < high_b
                assert low_a * low_b < number < high_a * high_b
                verified += 1
        assert verified == size // 8 - 1
        literal_states += verified
        width_controls.append({"scale": size, "number": number,
                               "number_bits": number.bit_length(),
                               "prime_factors": [first_prime, second_prime],
                               "row_interval": intervals[0], "column_interval": intervals[1],
                               "factor_count_in_rectangle": 0,
                               "last_unpruned_depth": exponent - 3,
                               "last_level_live_states": size // 16,
                               "verified_unpruned_states": verified})
    return {"literal_capped_valuation_checks": valuation_checks,
            "exact_binomial_baseline_checks": baseline_identity_checks,
            "nonzero_residual_height_checks": residual_bound_checks,
            "single_row_factor_slope_checks": slope_checks,
            "rectangle_factor_slope_checks": rectangle_checks,
            "huge_single_row_controls": huge_rows,
            "literal_unpruned_state_checks": literal_states,
            "rough_semiprime_width_controls": width_controls,
            "scope": "Exact factorial-baseline and inverse-window bridge; exponential width for the stated dyadic corner-pruning strategy, not a factoring lower bound."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    receipt = controls()
    payload = json.dumps(receipt, sort_keys=True, indent=2) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
