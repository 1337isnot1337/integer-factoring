#!/usr/bin/env python3
"""Exact controls for relative-Frobenius trace by vector carries.

The implementation uses sparse transfer states, with full polynomial and
finite-field enumeration retained only as independent finite oracles.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path
import random

from check_carry_fourier_elimination import BinaryField, geometric
from check_frobenius_rational_trace import field_of_degree


COUNTS = Counter()


def coefficient_vector_carries(field, t, factors, target):
    """Factors are (row, constant, coefficient-list by q-column)."""
    factor_count = len(factors)
    degree_bound = sum(max([0] + [(1 << (row + t * j))
                                  for j, value in enumerate(coefficients) if value])
                       for row, _, coefficients in factors)
    if target > degree_bound:
        return 0
    columns = max(1, (max(1, degree_bound.bit_length()) + t - 1) // t)
    local = []
    zero = (0,) * columns
    for row in range(t):
        values = {zero: 1}
        for factor_row, constant, coefficients in factors:
            if factor_row != row:
                continue
            choices = [(None, constant)] + [(j, value) for j, value in enumerate(coefficients)
                                             if j < columns and value]
            updated = {}
            for exponents, weight in values.items():
                for column, coefficient in choices:
                    if not coefficient:
                        continue
                    output = list(exponents)
                    if column is not None:
                        output[column] += 1
                    output = tuple(output)
                    updated[output] = updated.get(output, 0) ^ field.mul(weight, coefficient)
            values = {key: value for key, value in updated.items() if value}
        local.append(values)
    transition_cache = {}

    def transitions(row, state):
        key = row, state
        if key in transition_cache:
            return transition_cache[key]
        result = []
        for ones, weight in local[row].items():
            totals = tuple(state[j] + ones[j] for j in range(columns))
            if all((totals[j] & 1) == ((target >> (row + t * j)) & 1)
                   for j in range(columns)):
                output = tuple(total >> 1 for total in totals)
                assert all(carry <= factor_count for carry in output)
                result.append((output, weight))
        transition_cache[key] = result
        return result

    answer = 0
    for tail in product(range(factor_count + 1), repeat=columns - 1):
        initial = (0,) + tail
        state = {initial: 1}
        for row in range(t):
            updated = {}
            for incoming, weight in state.items():
                for outgoing, transfer_weight in transitions(row, incoming):
                    updated[outgoing] = (updated.get(outgoing, 0)
                                          ^ field.mul(weight, transfer_weight))
            state = {key: value for key, value in updated.items() if value}
        required_final = tail + (0,)
        answer ^= state.get(required_final, 0)
        COUNTS["initial_carry_vector_controls"] += 1
    COUNTS["computed_row_transfer_states"] += len(transition_cache)
    return answer


def ordinary_product(field, t, factors):
    result = [1]
    for row, constant, coefficients in factors:
        degree = max([0] + [1 << (row + t * j)
                            for j, value in enumerate(coefficients) if value])
        output = [0] * (len(result) + degree)
        for index, weight in enumerate(result):
            output[index] ^= field.mul(weight, constant)
            for column, coefficient in enumerate(coefficients):
                if coefficient:
                    exponent = 1 << (row + t * column)
                    output[index + exponent] ^= field.mul(weight, coefficient)
        result = output
    return result


def value_factor(field, t, factor, point):
    row, constant, coefficients = factor
    value = constant
    for column, coefficient in enumerate(coefficients):
        value ^= field.mul(coefficient, field.power(point, 1 << (row + t * column)))
    return value


def inverse_zero(field, value):
    return field.power(value, (1 << (field.degree + 1)) - 3)


def expand_frobenius_inverses(field, t, dimension, numerators, denominators):
    w = t * dimension
    if field.degree % w:
        raise ValueError("coefficient field does not contain summation field")
    factors = list(numerators)
    for row, constant, coefficients in denominators:
        for power in [0] + list(range(2, field.degree + 1)):
            new_row = (row + power) % t
            new_coefficients = [0] * dimension
            for column, coefficient in enumerate(coefficients):
                new_exponent = (row + power + t * column) % w
                assert new_exponent % t == new_row
                new_coefficients[new_exponent // t] ^= field.power(coefficient, 1 << power)
            factors.append((new_row, field.power(constant, 1 << power), new_coefficients))
    assert len(factors) == len(numerators) + field.degree * len(denominators)
    return factors


def relative_trace(field, t, dimension, scalar, numerators, denominators):
    factors = expand_frobenius_inverses(field, t, dimension, numerators, denominators)
    w, base = t * dimension, 1 << (t * dimension)
    degree_bound = sum(max([0] + [1 << (row + t * j)
                                  for j, value in enumerate(coefficients) if value])
                       for row, _, coefficients in factors)
    total = 0
    for target in range(base - 1, degree_bound + 1, base - 1):
        total ^= coefficient_vector_carries(field, t, factors, target)
        COUNTS["positive_multiplicative_order_coefficient_queries"] += 1
    total = field.mul(scalar, total)
    at_zero = scalar
    for _, constant, _ in numerators:
        at_zero = field.mul(at_zero, constant)
    for _, constant, _ in denominators:
        at_zero = field.mul(at_zero, inverse_zero(field, constant))
    return total, total ^ at_zero


def coefficient_controls():
    randomizer = random.Random(0xB10C5)
    field = field_of_degree(4)
    for t, dimension in [(1, 1), (1, 2), (2, 1), (2, 2), (3, 2)]:
        for _ in range(24):
            factors = [(randomizer.randrange(t), randomizer.randrange(field.size),
                        [randomizer.randrange(field.size) for _ in range(dimension)])
                       for _ in range(randomizer.randrange(1, 5))]
            full = ordinary_product(field, t, factors)
            targets = {0, len(full) - 1, len(full), randomizer.randrange(len(full))}
            for target in targets:
                expected = full[target] if target < len(full) else 0
                assert coefficient_vector_carries(field, t, factors, target) == expected
                COUNTS["exact_vector_coefficient_comparisons"] += 1
    # Four choices of U^8 add to U^32, carrying out of the second q=4 digit.
    factors = [(1, 0, [0, 1])] * 4
    assert coefficient_vector_carries(field, 2, factors, 32) == 1
    assert coefficient_vector_carries(field, 2, factors, 0) == 0
    assert coefficient_vector_carries(field, 2, factors, 33) == 0
    COUNTS["cross_digit_and_final_overflow_controls"] += 3


def trace_controls():
    randomizer = random.Random(0xA771AE)
    model_data = []
    for t, dimension, extension in [(1, 1, 1), (1, 1, 2), (1, 2, 1),
                                     (2, 1, 1), (2, 1, 2), (2, 2, 1)]:
        degree = t * dimension * extension
        field = BinaryField(0b11) if degree == 1 else field_of_degree(degree)
        base = 1 << (t * dimension)
        domain = [point for point in range(field.size) if field.power(point, base) == point]
        assert len(domain) == base
        for value in range(field.size):
            expected = field.inverse(value) if value else 0
            assert inverse_zero(field, value) == expected
            COUNTS["exact_zero_preserving_inverse_controls"] += 1
        for trial in range(40):
            def factor():
                return (randomizer.randrange(t), randomizer.randrange(field.size),
                        [randomizer.randrange(field.size) for _ in range(dimension)])
            numerators = [factor() for _ in range(randomizer.randrange(3))]
            denominators = [factor() for _ in range(randomizer.randrange(2))]
            if trial % 7 == 0:
                denominators = [(0, 0, [1] + [0] * (dimension - 1))]
            scalar = randomizer.randrange(field.size)
            computed = relative_trace(field, t, dimension, scalar, numerators, denominators)
            expanded = expand_frobenius_inverses(field, t, dimension, numerators, denominators)
            complete, punctured = 0, 0
            for point in domain:
                value = scalar
                for numerator in numerators:
                    value = field.mul(value, value_factor(field, t, numerator, point))
                for denominator in denominators:
                    denominator_value = value_factor(field, t, denominator, point)
                    if not denominator_value:
                        COUNTS["zero_denominator_point_controls"] += 1
                    value = field.mul(value, inverse_zero(field, denominator_value))
                polynomial_value = scalar
                for item in expanded:
                    polynomial_value = field.mul(polynomial_value, value_factor(field, t, item, point))
                assert polynomial_value == value
                complete ^= value
                if point:
                    punctured ^= value
                COUNTS["reference_field_point_evaluations"] += 1
            assert computed == (complete, punctured)
            COUNTS["exact_relative_trace_comparisons"] += 1
        for value in range(field.size):
            numerator = 1 ^ field.power(value, base)
            assert geometric(field, value, base) == field.mul(numerator, inverse_zero(field, 1 ^ value))
            COUNTS["all_point_geometric_box_controls"] += 1
        model_data.append({"t": t, "D": dimension, "field_degree": degree,
                           "irreducible_polynomial": bin(field.polynomial)})
    return model_data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    coefficient_controls()
    models = trace_controls()
    result = {
        "claim": "exact relative-Frobenius product traces by vector carries",
        "status": "passed",
        "scope": "The explicit product representation is required; no repeated network closure is claimed.",
        "field_models": models,
        "counts": dict(sorted(COUNTS.items())),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
