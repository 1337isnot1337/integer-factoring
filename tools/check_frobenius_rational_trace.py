#!/usr/bin/env python3
"""Exact replay of polynomial-cost Frobenius rational-star contraction.

Finite-field enumeration is a reference oracle only. The tested algorithm
uses pole traces and binary carries and does not enumerate the summation
field. This replay makes no claim about repeated network elimination.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import random

from check_carry_fourier_elimination import BinaryField, geometric, irreducible


COUNTS = Counter()


def field_of_degree(degree):
    for lower in range(1, 1 << degree, 2):
        polynomial = (1 << degree) | lower
        if irreducible(polynomial):
            return BinaryField(polynomial)
    raise AssertionError("irreducible polynomial not found")


def multiply_series(field, left, right, bound):
    result = [0] * (bound + 1)
    for i, value in enumerate(left):
        if value:
            for j, other in enumerate(right[:bound + 1 - i]):
                result[i + j] ^= field.mul(value, other)
    return result


def coefficient_reference(field, numerators, denominators, target):
    value = [1] + [0] * target
    for a, b, exponent in numerators:
        step = 1 << exponent
        factor = [0] * (target + 1)
        factor[0] = a
        if step <= target:
            factor[step] ^= b
        value = multiply_series(field, value, factor, target)
    for c, exponent in denominators:
        step = 1 << exponent
        factor = [0] * (target + 1)
        power = 1
        for index in range(0, target + 1, step):
            factor[index] = power
            power = field.mul(power, c)
        value = multiply_series(field, value, factor, target)
    return value[target]


def coefficient_carries(field, numerators, denominators, target):
    """Coefficient of product(a+b Z^(2^e))/product(1+c Z^(2^f))."""
    if target < 0:
        return 0
    bits = max(1, target.bit_length())
    scalar = 1
    for a, _, exponent in numerators:
        if exponent >= bits:
            scalar = field.mul(scalar, a)
    state = {0: scalar}
    factor_count = len(numerators) + len(denominators)
    for bit in range(bits):
        choices = [(1, field.power(c, 1 << (bit - exponent)))
                   for c, exponent in denominators if exponent <= bit]
        choices += [(a, b) for a, b, exponent in numerators if exponent == bit]
        local = [1]
        for a, b in choices:
            updated = [0] * (len(local) + 1)
            for index, value in enumerate(local):
                updated[index] ^= field.mul(a, value)
                updated[index + 1] ^= field.mul(b, value)
            local = updated
        updated = {}
        requested_bit = (target >> bit) & 1
        for carry, weight in state.items():
            for ones, local_weight in enumerate(local):
                total = carry + ones
                if total & 1 == requested_bit:
                    next_carry = total >> 1
                    assert next_carry <= factor_count
                    updated[next_carry] = (updated.get(next_carry, 0)
                                           ^ field.mul(weight, local_weight))
        state = updated
        COUNTS["binary_carry_positions"] += 1
    return state.get(0, 0)


def rational_value(field, scalar, numerators, denominators, point):
    numerator, denominator = scalar, 1
    for a, b, exponent in numerators:
        numerator = field.mul(numerator,
                              a ^ field.mul(b, field.power(point, 1 << exponent)))
    for c, d, exponent in denominators:
        denominator = field.mul(denominator,
                                c ^ field.mul(d, field.power(point, 1 << exponent)))
    if not denominator:
        raise ValueError("pole on the evaluation domain")
    return field.mul(numerator, field.inverse(denominator))


def rational_trace(field, w, scalar, numerators, denominators):
    """Return both complete-field and nonzero-field sums, without enumeration."""
    if field.degree % w:
        raise ValueError("the coefficient field must contain the summation field")
    if any(not (0 <= exponent < w) for _, _, exponent in numerators + denominators):
        raise ValueError("initial Frobenius exponents outside the promised range")
    if sum(1 << e for _, _, e in numerators) > sum(1 << e for _, _, e in denominators):
        raise ValueError("positive-degree polynomial part is not supported")
    base = 1 << w
    grouped = {}
    roots = []
    for c, d, exponent in denominators:
        if not d:
            raise ValueError("denominator slope must be nonzero")
        ratio = field.mul(c, field.inverse(d))
        root = field.power(ratio, 1 << ((-exponent) % field.degree))
        assert field.power(root, 1 << exponent) == ratio
        if field.power(root, base) == root:
            raise ValueError("denominator root lies in the summation field")
        roots.append(root)
        grouped[root] = grouped.get(root, 0) + (1 << exponent)
    total = 0
    for root, multiplicity in grouped.items():
        shifted_numerators = [
            (a ^ field.mul(b, field.power(root, 1 << exponent)), b, exponent)
            for a, b, exponent in numerators
        ]
        shifted_denominators = []
        shifted_scalar = scalar
        for (c, d, exponent), other_root in zip(denominators, roots):
            if other_root == root:
                shifted_scalar = field.mul(shifted_scalar, field.inverse(d))
            else:
                constant = c ^ field.mul(d, field.power(root, 1 << exponent))
                assert constant
                inverse = field.inverse(constant)
                shifted_scalar = field.mul(shifted_scalar, inverse)
                shifted_denominators.append((field.mul(d, inverse), exponent))
        delta = field.power(root, base) ^ root
        assert delta
        if multiplicity >= base:
            COUNTS["poles_with_multiplicity_at_least_base"] += 1
        for v in range((multiplicity - 1) // base + 1):
            denominator = shifted_denominators[:]
            term_scalar = shifted_scalar
            for bit in range((v + 1).bit_length()):
                if ((v + 1) >> bit) & 1:
                    constant = field.power(delta, 1 << bit)
                    inverse = field.inverse(constant)
                    term_scalar = field.mul(term_scalar, inverse)
                    denominator.append((inverse, bit))
            target = multiplicity - 1 - base * v
            coefficient = coefficient_carries(field, shifted_numerators,
                                               denominator, target)
            total ^= field.mul(term_scalar, coefficient)
            COUNTS["pole_trace_coefficient_queries"] += 1
            if v:
                COUNTS["higher_Z_to_base_trace_terms"] += 1
    at_zero = rational_value(field, scalar, numerators, denominators, 0)
    return total, total ^ at_zero


def coefficient_controls():
    randomizer = random.Random(0xC0EFF1C1)
    field = field_of_degree(6)
    for trial in range(320):
        numerator_count = randomizer.randrange(8)
        denominator_count = randomizer.randrange(8)
        numerators = [(randomizer.randrange(field.size), randomizer.randrange(field.size),
                       randomizer.randrange(7)) for _ in range(numerator_count)]
        denominators = [(randomizer.randrange(field.size), randomizer.randrange(7))
                         for _ in range(denominator_count)]
        if numerators and trial % 3 == 0:
            _, b, exponent = numerators[0]
            numerators[0] = (0, b, exponent)
            COUNTS["zero_numerator_constant_coefficient_controls"] += 1
        targets = [0, 1, 2, 7, 15, randomizer.randrange(65)]
        for target in targets:
            assert (coefficient_carries(field, numerators, denominators, target)
                    == coefficient_reference(field, numerators, denominators, target))
            COUNTS["exact_coefficient_comparisons"] += 1
    # The coefficient at seven requires a nonzero carry, including a
    # carry across a position with no new numerator factor.
    numerators = [(1, 1, 0)] * 7
    assert coefficient_carries(field, numerators, [], 7) == 1
    assert coefficient_carries(field, numerators, [], 8) == 0
    COUNTS["carry_and_overflow_controls"] += 2


def trace_controls():
    randomizer = random.Random(0x7ACE2026)
    models = []
    for w in range(1, 5):
        field = field_of_degree(2 * w)
        base = 1 << w
        subfield = [x for x in range(field.size) if field.power(x, base) == x]
        outside = [x for x in range(field.size) if field.power(x, base) != x]
        assert len(subfield) == base
        models.append({"w": w, "field_degree": field.degree,
                       "irreducible_polynomial": bin(field.polynomial)})
        for trial in range(100):
            count = randomizer.randrange(1, 8)
            roots = randomizer.choices(outside, k=count)
            if trial % 4 == 0:
                roots = [outside[0]] * count
            exponents = [randomizer.randrange(w) for _ in range(count)]
            if trial % 4 == 0:
                exponents = [w - 1] * count
            denominators = []
            for root, exponent in zip(roots, exponents):
                slope = randomizer.randrange(1, field.size)
                constant = field.mul(slope, field.power(root, 1 << exponent))
                denominators.append((constant, slope, exponent))
            numerator_count = randomizer.randrange(count + 1)
            numerators = [(randomizer.randrange(field.size),
                           randomizer.randrange(field.size), exponents[index])
                          for index in range(numerator_count)]
            if numerators and trial % 5 == 0:
                numerators[0] = denominators[0]
                COUNTS["cancelled_pole_controls"] += 1
            if numerators and trial % 7 == 0:
                _, slope, exponent = numerators[0]
                numerators[0] = (0, slope, exponent)
                COUNTS["zero_numerator_constant_trace_controls"] += 1
            scalar = randomizer.randrange(field.size)
            computed = rational_trace(field, w, scalar, numerators, denominators)
            complete, punctured = 0, 0
            for point in subfield:
                value = rational_value(field, scalar, numerators, denominators, point)
                complete ^= value
                if point:
                    punctured ^= value
                COUNTS["reference_field_point_evaluations"] += 1
            assert computed == (complete, punctured)
            COUNTS["exact_rational_trace_comparisons"] += 1
            if numerator_count < count:
                COUNTS["proper_rational_controls"] += 1
            if len(set(roots)) == len(roots):
                COUNTS["pairwise_distinct_pole_controls"] += 1
        assert rational_trace(field, w, 1, [], []) == (0, 1)
        COUNTS["constant_rational_controls"] += 1
        for root in subfield:
            try:
                rational_trace(field, w, 1, [], [(root, 1, 0)])
            except ValueError as error:
                assert "summation field" in str(error)
                COUNTS["forbidden_domain_pole_rejections"] += 1
            else:
                raise AssertionError("failed to reject a pole in the summation field")
    return models


def kernel_controls():
    for w in range(1, 5):
        field = field_of_degree(2 * w)
        base = 1 << w
        units = [u for u in range(1, field.size) if field.power(u, base) == u]
        circle = [v for v in range(1, field.size) if field.power(v, base + 1) == 1]
        assert len(units) == base - 1 and len(circle) == base + 1
        kernel = {(u, v): geometric(field, field.mul(u, v), base)
                  for u in units for v in circle}
        for u in units:
            total = 0
            for v in circle:
                total ^= kernel[u, v]
                denominator = field.mul(v, 1 ^ field.mul(u, v))
                if denominator:
                    assert kernel[u, v] == field.mul(u ^ v, field.inverse(denominator))
                    COUNTS["exact_Cauchy_box_comparisons"] += 1
                else:
                    assert u == v == 1 and kernel[u, v] == 0
                    COUNTS["removable_pole_controls"] += 1
            assert total == 1
            COUNTS["single_unit_circle_contractions"] += 1
        for v in circle:
            total = 0
            for u in units:
                total ^= kernel[u, v]
            assert total == 1 ^ field.power(field.inverse(v), 2)
            COUNTS["single_subfield_unit_contractions"] += 1
        for left in circle:
            for right in circle:
                total = 0
                for u in units:
                    total ^= field.mul(kernel[u, left], kernel[u, right])
                inverse = field.inverse(field.mul(left, right))
                expected = 1 ^ field.power(inverse, 2)
                if left != right:
                    expected ^= inverse
                assert total == expected
                COUNTS["double_subfield_unit_contractions"] += 1
        for left in units:
            for right in units:
                total = 0
                for v in circle:
                    total ^= field.mul(kernel[left, v], kernel[right, v])
                expected = 1 ^ field.mul(left, right)
                if left == right:
                    expected ^= field.power(left, 2)
                assert total == expected
                COUNTS["double_unit_circle_contractions"] += 1
        randomizer = random.Random(w)
        for _ in range(40):
            gamma = randomizer.randrange(1, field.size)
            u, v = randomizer.choice(units), randomizer.choice(circle)
            exponent = randomizer.randrange(w)
            upower = field.power(u, 1 << exponent)
            value = field.product(gamma, upower, v)
            if value == 1:
                continue
            numerator = 1 ^ field.product(field.power(gamma, base), upower,
                                           field.inverse(v))
            assert (geometric(field, value, base)
                    == field.mul(numerator, field.inverse(1 ^ value)))
            COUNTS["tagged_Frobenius_box_comparisons"] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    coefficient_controls()
    models = trace_controls()
    kernel_controls()
    result = {
        "claim": "exact polynomial-operation contraction of one Frobenius rational star",
        "status": "passed",
        "scope": "No second-variable closure or improved factoring exponent is claimed.",
        "field_models": models,
        "counts": dict(sorted(COUNTS.items())),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
