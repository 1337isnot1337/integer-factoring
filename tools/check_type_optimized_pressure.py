#!/usr/bin/env python3
"""Exact arbitrary-type pressure certificates and physical-gauge controls."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path

import check_block_smith_crossings as smith
from check_schur_gaussian_reduction import seed_data


def profile_at_end(local, heights, cut, infinity, degree):
    changed = []
    for labels, matrix in local.items():
        value = matrix[0][0]
        if not value:
            continue
        height = sum(h[label] for h, label in zip(heights, labels))
        if infinity:
            value = smith.reverse_polynomial(value, degree)
            height -= degree
        order = (value & -value).bit_length() - 1
        changed.append((labels, value >> order, height + order))
    offset = min(order for _, _, order in changed)
    gauged = {labels: [[0]] for labels in local}
    for labels, unit, order in changed:
        gauged[labels] = [[unit << (order - offset)]]
    matrix = smith.flatten(gauged, 4, cut)
    bound = max(value.bit_length() - 1 for row in matrix for value in row)
    orders, precision = smith.smith_profile(matrix, bound)
    checked = 0
    for ell in range(1, max(orders) + 2):
        assert smith.toeplitz_rank(matrix, ell) == sum(max(0, ell - value) for value in orders)
        checked += 1
    return tuple(order + offset for order in orders), checked, precision


def rational_type_certificate(orders_infinity, orders_zero, weights, root, numerator, denominator):
    """Use s=root^numerator and R=numerator/denominator; all arithmetic exact."""
    assert root > 1 and numerator > 0 and denominator > 0
    s = root ** numerator
    raw = sum(s ** (-v) for v in orders_infinity) * sum(s ** (-v) for v in orders_zero)
    tilted = tuple(root ** (-denominator * g) for g in weights)
    total = sum(tilted)
    probabilities = tuple(value / total for value in tilted)
    lhs, rhs = raw ** denominator, total ** numerator
    assert sum(probabilities) == 1 and all(value > 0 for value in probabilities)
    return {'passes': lhs < rhs, 'type': probabilities, 'raw': raw, 's': s,
            'lhs': lhs, 'rhs': rhs, 'ratio': lhs / rhs}


def compositions(total, parts):
    if parts == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def check_type_counts():
    count = 0
    for q in range(2, 6):
        for length in range(q, 11):
            for counts in compositions(length, q):
                dimension = factorial(length)
                denominator = 1
                for value in counts:
                    dimension //= factorial(value)
                    denominator *= value ** value
                entropy_power = Fraction(length ** length, denominator)
                assert dimension <= entropy_power <= (length + 1) ** q * dimension
                count += 1
    return count


def check():
    seed = seed_data()
    local = {tuple(row['indices']): row['matrix'] for row in seed['local']}
    infinity_heights = ((0, 1, 1, 2), (0, 1, 1, 2), (0, 0, 0, 0), (0, 2, 2, 4))
    zero_heights = ((0, 0, 0, 0), (0, 0, 0, 0), (0, 1, 1, 2), (0, -1, -1, -2))
    first, count_a, _ = profile_at_end(local, infinity_heights, 1, True, 8)
    second, count_b, _ = profile_at_end(local, zero_heights, 0, False, 8)
    expected = (0, 1, 1, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 7, 7, 8)
    assert first == second == expected
    weights = tuple(sum(h[a] for h in infinity_heights + zero_heights) for a in range(4))
    assert weights == (0, 4, 4, 8)
    for s in map(Fraction, ('3/2', '2', '5/2', '3')):
        assert sum(s ** (-v) for v in first) == (1 + 1 / s) ** 2 * (1 + s ** -3) ** 2

    # Small-denominator common type used in the displayed all-length bound.
    pi = tuple(Fraction(value, 81) for value in (49, 14, 14, 4))
    assert sum(pi) == 1
    assert sum(p * g for p, g in zip(pi, weights)) == Fraction(16, 9)
    J = Fraction(5, 2) ** 16 * Fraction(7, 5) ** 36 * Fraction(133, 125) ** 36
    Q = Fraction(9 ** 9, 2 ** 2 * 7 ** 7)
    assert J ** 10 < Q ** 61

    passing = rational_type_certificate(first, second, weights, Fraction(33, 32), 31, 10)
    failing = rational_type_certificate(first, second, weights, Fraction(17, 16), 14, 5)
    assert passing['passes']
    assert not failing['passes']
    # A uniform three-power shift changes raw spectra and means together.
    shifted = rational_type_certificate(tuple(v + 3 for v in first), second,
                                        tuple(g + 3 for g in weights), Fraction(33, 32), 31, 10)
    assert shifted['ratio'] == passing['ratio']
    assert shifted['type'] == passing['type']
    # Forgetting the mean shift would fabricate a post-half certificate.
    false_pass = rational_type_certificate(tuple(v + 10 for v in first), second,
                                           weights, Fraction(17, 16), 14, 5)
    assert false_pass['passes']

    def summary(result, root, numerator, denominator):
        return {'root': str(root), 'R': str(Fraction(numerator, denominator)),
                'passes': result['passes'], 'rational_type': list(map(str, result['type'])),
                'lhs_numerator_bits': result['lhs'].numerator.bit_length(),
                'rhs_numerator_bits': result['rhs'].numerator.bit_length()}

    return {
        'status': 'pass', 'infinity_cut': 1, 'zero_cut': 0,
        'infinity_heights': infinity_heights, 'zero_heights': zero_heights,
        'restored_orders_in_both_cuts': first, 'total_letter_heights': weights,
        'independent_toeplitz_ranks': count_a + count_b,
        'type_dimension_controls': check_type_counts(),
        'explicit_type': list(map(str, pi)), 'explicit_s': '5/2',
        'exact_rho_upper': '61/80', 'exact_integer_comparison_J10_lt_Q61': True,
        'rational_variational_control': summary(passing, Fraction(33, 32), 31, 10),
        'post_half_negative_control': summary(failing, Fraction(17, 16), 14, 5),
        'uniform_shift_invariance': True, 'omitted_mean_false_pass_detected': True,
        'scope': 'Exact general type criterion and a weaker-than-half gauged seed bound; no score improvement.',
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = check()
    source_files = (Path(__file__), Path(smith.__file__), Path(__file__).with_name('check_schur_gaussian_reduction.py'))
    report['source_sha256'] = {p.name: sha256(p.read_bytes()).hexdigest() for p in source_files}
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
