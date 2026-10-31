#!/usr/bin/env python3
"""Exact controls for read-once homogeneous Boolean filter obstructions.

No floating-point sign decision is used. The uniform results are proved in
research/read-once-boolean-filter-floor.md; this checks independent finite
transcriptions of the probability, derivative, and polynomial identities.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product
import json
from math import prod
from pathlib import Path
import random


@cache
def formulas(variables):
    """All homogeneous read-once XOR/AND polynomials on this variable set.

    Keys are sorted monomial masks. One expression for each polynomial is
    retained; multiplication of disjoint monomials has no cancellations.
    """
    if variables & (variables - 1) == 0:
        return {(variables,): ('leaf', variables.bit_length() - 1)}
    least = variables & -variables
    answer = {}
    split = (variables - 1) & variables
    while split:
        other = variables ^ split
        if other and split & least:
            for left, left_tree in formulas(split).items():
                for right, right_tree in formulas(other).items():
                    joined = tuple(sorted(a | b for a in left for b in right))
                    answer.setdefault(joined, ('and', left_tree, right_tree))
                    if left[0].bit_count() == right[0].bit_count():
                        joined = tuple(sorted(left + right))
                        answer.setdefault(joined, ('xor', left_tree, right_tree))
        split = (split - 1) & variables
    return answer


def truth_value(monomials, state):
    return sum((state & monomial) == monomial for monomial in monomials) % 2


def truth_probability(monomials, probabilities):
    total = F(0)
    derivative = F(0)
    mean = sum(probabilities, F(0))
    for state in range(1 << len(probabilities)):
        if truth_value(monomials, state):
            weight = prod(p if (state >> i) & 1 else 1 - p
                          for i, p in enumerate(probabilities))
            total += weight
            derivative += weight * (state.bit_count() - mean)
    return total, derivative


def tree_probability(tree, probabilities):
    if tree[0] == 'leaf':
        p = probabilities[tree[1]]
        return p, p * (1 - p)
    p, dp = tree_probability(tree[1], probabilities)
    q, dq = tree_probability(tree[2], probabilities)
    if tree[0] == 'and':
        return p * q, dp * q + p * dq
    assert tree[0] == 'xor'
    return p + q - 2 * p * q, dp * (1 - 2 * q) + dq * (1 - 2 * p)


def exact_ratio(p, q, h, a):
    return (1 - q + q * a) * (1 - p + p * h / a)


def polynomial_control(monomials, modes, m, increment):
    degree = monomials[0].bit_count()
    final = m * degree + increment
    original = [1 ^ (((1 << (m * degree)) ^ (1 << final))
                     if truth_value(monomials, state) else 0)
                for state in range(1 << modes)]
    assert all(value.bit_count() % 2 == 1 for value in original)
    values = original[:]
    for axis in range(modes):
        bit = 1 << axis
        for state in range(1 << modes):
            if state & bit:
                source = values[state ^ bit]
                values[state] ^= source ^ (source << m)
    expected = [1 << (m * state.bit_count() +
                      increment * truth_value(monomials, state))
                for state in range(1 << modes)]
    assert values == expected
    # Independent evaluation of the two profile sums is done at x=3/2.
    # The normalizing product is the binary product-measure partition.
    base = F(3, 2)
    energies = tuple(2 * i + 1 for i in range(modes))
    charge = -2
    unshifted = tuple(1 / (1 + base ** e) for e in energies)
    shifted = tuple(1 / (1 + base ** (e + charge)) for e in energies)
    p = truth_probability(monomials, unshifted)[0]
    q = truth_probability(monomials, shifted)[0]
    zero_sum = F(0)
    infinity_sum = F(0)
    for state in range(1 << modes):
        zero_weight = prod(u if (state >> i) & 1 else 1 - u
                           for i, u in enumerate(shifted))
        infinity_weight = prod(u if (state >> i) & 1 else 1 - u
                               for i, u in enumerate(unshifted))
        order = (values[state] & -values[state]).bit_length() - 1
        full_degree = original[state].bit_length() - 1
        zero_sum += zero_weight * base ** (2 * (m * state.bit_count() - order))
        infinity_sum += infinity_weight * base ** (2 * full_degree)
    assert zero_sum * infinity_sum == exact_ratio(
        p, q, base ** (2 * m * degree), base ** (-2 * increment))
    return len(original)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=root / 'evidence' / 'read-once-boolean-filter-floor.json')
    args = parser.parse_args()
    grid = tuple(F(a, b) for b in range(2, 9) for a in range(1, b))
    gate_controls = 0
    for p, q in product(grid, repeat=2):
        w = p + q - 2 * p * q
        lhs = (abs(p * (1 - p) * (1 - 2 * q)) +
               abs(q * (1 - q) * (1 - 2 * p))) / (w * (1 - w))
        a, b = abs(1 - 2 * p), abs(1 - 2 * q)
        assert lhs == (a + b) / (1 + a * b) <= 1
        assert 0 < (1 - p) / (1 - p * q) < 1
        assert 0 < (1 - q) / (1 - p * q) < 1
        gate_controls += 1

    rng = random.Random(0xA11C0DE)
    probability_controls = derivative_controls = odds_controls = ratio_controls = 0
    transform_controls = endpoint_controls = profile_controls = 0
    family_counts = []
    for modes in range(1, 6):
        family = formulas((1 << modes) - 1)
        family_counts.append({'modes': modes, 'distinct_polynomials': len(family)})
        n = 1 << (modes - 1).bit_length()
        for monomials, tree in sorted(family.items()):
            degree = monomials[0].bit_count()
            assert all(m.bit_count() == degree for m in monomials)
            for _ in range(5):
                probabilities = tuple(rng.choice(grid) for _ in range(modes))
                p, dp = truth_probability(monomials, probabilities)
                assert (p, dp) == tree_probability(tree, probabilities)
                assert 0 < p < 1
                assert abs(dp) <= degree * p * (1 - p)
                probability_controls += 1
                derivative_controls += 1
            selected = tuple(sorted(rng.sample(range(n), modes)))
            energies = tuple(2 * n - 1 - 2 * j for j in selected)
            for base in (F(5, 4), F(3, 2), F(2)):
                unshifted = tuple(1 / (1 + base ** e) for e in energies)
                p = truth_probability(monomials, unshifted)[0]
                h = base ** (2 * n * degree)
                for charge in (-2 * n, -n, -1, 0, 1, n, 2 * n):
                    shifted = tuple(1 / (1 + base ** (e + charge)) for e in energies)
                    q = truth_probability(monomials, shifted)[0]
                    assert q * (1 - p) <= h * p * (1 - q)
                    odds_controls += 1
                    for increment in (1, 3):
                        ratio = exact_ratio(p, q, h, base ** (-2 * increment))
                        assert ratio >= 1 - p + p * h > 1
                        ratio_controls += 1
            if modes <= 4:
                for m, increment in ((1, 1), (3, 2), (n, 3)):
                    endpoint_controls += polynomial_control(monomials, modes, m, increment)
                    transform_controls += 1
                    profile_controls += 1

    exterior_controls = parity_controls = 0
    for n in (1, 2, 4, 8):
        for modes in range(1, min(n, 5) + 1):
            for selected in list(combinations(range(n), modes))[:12]:
                monomials = tuple(1 << i for i in range(modes))
                energies = tuple(2 * n - 1 - 2 * j for j in selected)
                for base in (F(17, 16), F(3, 2), F(2)):
                    probabilities = tuple(1 / (1 + base ** e) for e in energies)
                    p = truth_probability(monomials, probabilities)[0]
                    assert 1 - 2 * p == prod(1 - 2 * u for u in probabilities)
                    assert p >= max(probabilities)
                    parity_controls += 1
                    h = base ** (2 * n)
                    assert p * h > 1 / (1 + 1 / base)
                    for charge in (-8 * n, -2 * n - 1, 2 * n + 1, 8 * n):
                        shifted = tuple(1 / (1 + base ** (e + charge)) for e in energies)
                        q = truth_probability(monomials, shifted)[0]
                        for increment in (1, 2 * n):
                            ratio = exact_ratio(p, q, h, base ** (-2 * increment))
                            assert ratio > p * h > 1 / (1 + 1 / base)
                            exterior_controls += 1

    # Relative gains at exterior charges are possible, and must not be
    # confused with a strict full-pressure witness.
    n, base, selected, charge, increment = 4, F(17, 16), (0, 1, 2), -100, 20
    energies = tuple(2 * n - 1 - 2 * j for j in selected)
    monomials = (1, 2, 4)
    p = truth_probability(monomials, tuple(1 / (1 + base ** e) for e in energies))[0]
    q = truth_probability(monomials, tuple(1 / (1 + base ** (e + charge)) for e in energies))[0]
    ratio = exact_ratio(p, q, base ** (2 * n), base ** (-2 * increment))
    assert ratio < 1
    assert ratio > 1 / (1 + 1 / base)
    result = {
        'status': 'All controls passed using exact rational and binary-polynomial arithmetic.',
        'source_baseline': 'c2e951fa1f391edd97ba587d15ec13d3d6c8892d',
        'scope': {
            'uniform_proof': 'Every homogeneous read-once XOR/AND filter is strictly worse than its Gaussian at |J|<=2n; every homogeneous degree-one filter has positive full pressure at all charges.',
            'not_claimed': 'No factoring exponent improvement; no general non-read-once degree bound; no all-charge degree>=2 obstruction.'
        },
        'rational_gate_jacobian_controls': gate_controls,
        'enumerated_homogeneous_read_once_families': family_counts,
        'recursive_vs_truth_probability_controls': probability_controls,
        'rational_common_shift_derivative_controls': derivative_controls,
        'rational_interior_odds_controls': odds_controls,
        'rational_interior_ratio_controls': ratio_controls,
        'complete_polynomial_transform_controls': transform_controls,
        'exact_filter_endpoint_equalities': endpoint_controls,
        'independent_profile_ratio_controls': profile_controls,
        'exact_parity_product_controls': parity_controls,
        'rational_exterior_loss_cap_controls': exterior_controls,
        'strict_exterior_relative_gain': {
            'n': n, 'base_exp_t': [base.numerator, base.denominator],
            'selected': list(selected), 'charge': charge,
            'Boolean_monomials': list(monomials), 'final_degree': n + increment,
            'ratio_numerator': str(ratio.numerator), 'ratio_denominator': str(ratio.denominator),
            'exact_comparison': '1/(1+exp(-t)) < ratio < 1; the exterior Gaussian floor makes the full pressure positive.'
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key != 'strict_exterior_relative_gain'}, indent=2))


if __name__ == '__main__':
    main()
