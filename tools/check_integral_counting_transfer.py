#!/usr/bin/env python3
"""Independent finite checks for integral factor queries and exact count recovery."""

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from math import comb
from pathlib import Path

from integral_schur_queries import (
    schur_mode, translated_angle_mode, translated_power_sum,
)
from check_partition_channels import dual_polynomial, operator, partitions, partition_denominator
from check_integral_schur_crossing import schur


def reference_mode(parameter, output, input_, index=None):
    """Literal cached source recursion and separate Jacobi--Trudi expansions."""
    total = Fraction(0)
    output_character = {mu: coefficient * partition_denominator(mu)
                        for mu, coefficient in schur(output)}
    for mu, first in schur(parameter):
        lowered = ((mu, 1),) if index is None else dual_polynomial(mu)
        for term, translation in lowered:
            angle = 1 if index is None else (
                comb(sum(parameter) + sum(term), index)
                if sum(parameter) + sum(term) >= index else 0)
            for nu, second in schur(input_):
                for rho, character in output_character.items():
                    total += first * translation * angle * second * character * operator(term, rho, nu)
    assert total.denominator == 1
    return total.numerator


def least_factor_from_counts(number, counter):
    """Recover the least divisor greater than one from exact prefix residues."""
    bits = number.bit_length()
    modulus = 1 << (bits + 1)
    prefix = 0
    for remaining in range(bits - 1, -1, -1):
        candidate = 2 * prefix
        residue = counter(number, remaining, candidate, modulus)
        assert 0 <= residue < (1 << bits)
        without_one = residue - int(1 >> remaining == candidate)
        assert without_one >= 0
        prefix = candidate if without_one else candidate + 1
    assert number % prefix == 0 and prefix > 1
    return prefix


def verify():
    translations = 0
    for weight in range(8):
        for mu in partitions(weight):
            original = {nu: value * (-1) ** (weight - len(mu))
                        for nu, value in dual_polynomial(mu)}
            for target_weight in range(weight + 1):
                for nu in partitions(target_weight):
                    assert translated_power_sum(mu, nu) == original.get(nu, 0)
                    translations += 1

    # Caches below belong only to the small finite harness.  The production
    # oracles imported above never retain a table of partition values.
    mode = lru_cache(None)(schur_mode)
    angle = lru_cache(None)(translated_angle_mode)
    factor_checks = 0
    for weight in (1, 2):
        for energy in range(4 * weight + 1):
            for parameter, state, intermediate in product(
                    partitions(weight), partitions(weight), partitions(energy)):
                for output, input_ in ((intermediate, state), (state, intermediate)):
                    assert mode(parameter, output, input_) == reference_mode(parameter, output, input_)
                    factor_checks += 1
                    for index in range(2 * weight + 1):
                        assert angle(parameter, index, output, input_) == reference_mode(
                            parameter, output, input_, index)
                        factor_checks += 1

    piece_entries = crossing_entries = zero_divisor_entries = 0
    moduli = (4, 8, 9, 16, 25)
    for weight in (1, 2):
        basis = partitions(weight)
        for x, y, v, u in product(basis, repeat=4):
            total = 0
            for shift in range(-2 * weight, 2 * weight + 1):
                row_value = column_value = 0
                for index in range(2 * weight + 1):
                    for sigma in partitions(shift + index):
                        row_value += ((-1) ** index * angle(x, index, sigma, y)
                                      * mode(v, u, sigma))
                    for sigma in partitions(2 * weight - shift - index):
                        column_value += ((-1) ** index * angle(x, index, u, sigma)
                                         * mode(v, sigma, y))
                assert row_value == column_value
                for modulus in moduli:
                    assert row_value % modulus == column_value % modulus
                    zero_divisor_entries += 1
                total += row_value
                piece_entries += 1
            assert total == int(x == v and y == u)
            crossing_entries += 1

    prefix_queries = recoveries = 0
    for number in range(2, 4097):
        divisors = [d for d in range(1, number + 1) if number % d == 0]

        def counter(value, remaining, prefix, modulus):
            nonlocal prefix_queries
            assert value == number
            count = sum(d >> remaining == prefix for d in divisors)
            assert count <= number < modulus
            prefix_queries += 1
            return count % modulus

        actual = least_factor_from_counts(number, counter)
        assert actual == divisors[1]
        recoveries += 1

    # A fixed small characteristic can erase a positive count.  The large
    # residue modulus avoids that ambiguity without assuming it is a field.
    erased_count = sum(12 % d == 0 for d in range(1, 13))
    assert erased_count == 6 and erased_count % 2 == 0
    assert 2 * 4 % 8 == 0 and 2 % 8 and 4 % 8
    return {"scope": "integral factor and exact prefix-recovery checks; no full-size runtime experiment",
            "independent_translated_monomial_coefficients": translations,
            "independent_integral_factor_entries": factor_checks,
            "both_operator_order_piece_entries": piece_entries,
            "integer_crossing_entries": crossing_entries,
            "composite_ring_piece_entries": zero_divisor_entries,
            "composite_moduli": list(moduli),
            "least_factor_recoveries": recoveries,
            "exact_prefix_queries": prefix_queries,
            "negative_controls": {"small_characteristic_can_erase_counts": True,
                                  "nontrivial_zero_divisors_allowed": True}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    report = verify()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        Path(args.output).write_text(rendered)
    print(rendered, end="")
