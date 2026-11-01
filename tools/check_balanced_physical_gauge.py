#!/usr/bin/env python3
"""Exact balanced-type physical-gauge and shifted-coefficient checks."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path

import check_smith_profile_transfer as basic
import check_auxiliary_gauge_smith as auxiliary


def transform(local, heights, infinity=False):
    sign = -1 if infinity else 1
    dimension = len(next(iter(local.values())))
    return {
        labels: [[{sign*k + sum(height[label] for height, label in zip(heights, labels))
                   for k in auxiliary.exponents(value[a][b])}
                  for b in range(dimension)] for a in range(dimension)]
        for labels, value in local.items()
    }


def mean(heights):
    return Fraction(sum(sum(height) for height in heights), 2)


def balanced_words(length):
    assert length % 2 == 0
    return [word for word in range(1 << length) if word.bit_count() == length//2]


def restricted_coefficients(local, degree, length):
    labels = balanced_words(length)
    dimension = len(next(iter(local.values())))
    size = len(labels)
    rows = [[0]*(size*size) for _ in range(degree*length+1)]
    columns = [[0]*(size*size) for _ in range(degree*length+1)]
    products = {}
    for xi, yi, ui, vi in itertools.product(range(size), repeat=4):
        x, y, u, v = (labels[j] for j in (xi, yi, ui, vi))
        matrix = [[int(a == b) for b in range(dimension)] for a in range(dimension)]
        word = []
        for position in range(length):
            site = tuple(label >> position & 1 for label in (x, y, u, v))
            word.append(site)
            matrix = basic.matrix_product(matrix, local[site])
        polynomial = matrix[0][0]
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        products[tuple(word)] = polynomial
        for k in auxiliary.exponents(polynomial):
            rows[k][xi*size+yi] |= 1 << (ui*size+vi)
            columns[k][xi*size+ui] |= 1 << (yi*size+vi)
    return [list(map(basic.rank, rows)), list(map(basic.rank, columns))], products


def check_family(name, local, degree):
    weight_sets = [((0, 0),)*4,
                   ((0, 1), (0, -1), (0, 0), (0, 0)),
                   ((0, 1), (0, 0), (0, 0), (0, 0)),
                   ((-1, 1), (0, 2), (1, 0), (0, -1)),
                   ((1, 1), (0, 0), (0, 0), (0, 0))]
    profiles = {}
    toeplitz_checks = 0
    for wi, weights in enumerate(weight_sets):
        for infinity, cut, side in itertools.product((False, True), range(2), range(2)):
            changed = transform(local, weights, infinity)
            orders, checks = auxiliary.restored_orders(auxiliary.flatten(changed, cut, side))
            profiles[infinity, cut, side, wi] = [-x for x in orders] if infinity else orders
            toeplitz_checks += checks
    shift_entries = rank_checks = rational_checks = endpoint_entries = 0
    for length in (2, 4):
        ranks, products = restricted_coefficients(local, degree, length)
        endpoint_entries += len(products)
        dimension = len(next(iter(local.values())))
        for wi, weights in enumerate(weight_sets):
            shift = length*mean(weights)
            assert shift.denominator == 1
            for infinity in (False, True):
                changed = transform(local, weights, infinity)
                sign = -1 if infinity else 1
                for word, polynomial in products.items():
                    # Exact multiplication of transformed Laurent matrices;
                    # this does not infer the result from the mean identity.
                    matrix = auxiliary.identity_series(dimension)
                    for site in word:
                        matrix = auxiliary.product_matrix(matrix, changed[site])
                    assert matrix[0][0] == {sign*k + int(shift) for k in auxiliary.exponents(polynomial)}
                    shift_entries += 1
        for first_cut, first_side, second_side, wi, wj in itertools.product(
                range(2), range(2), range(2), range(len(weight_sets)), range(len(weight_sets))):
            degrees = profiles[True, first_cut, first_side, wi]
            orders = profiles[False, 1-first_cut, second_side, wj]
            H_inf, H_zero = mean(weight_sets[wi]), mean(weight_sets[wj])
            high, low = basic.word_counts(degrees, length), basic.word_counts(orders, length)
            for k in range(degree*length+1):
                upper_threshold = k-length*H_inf
                lower_threshold = k+length*H_zero
                r_bound = sum(count*(value-upper_threshold+1) for value, count in high.items() if value >= upper_threshold)
                c_bound = sum(count*(lower_threshold-value+1) for value, count in low.items() if value <= lower_threshold)
                assert ranks[first_cut][k] <= r_bound
                assert ranks[1-first_cut][k] <= c_bound
                rank_checks += 2
                for root in (Fraction(1), Fraction(3, 2), Fraction(2)):
                    # s=root^2 clears all binary mean denominators.
                    s = root**2
                    U = sum(root**int(2*d+2*H_inf) for d in degrees)
                    V = sum(root**int(-2*v+2*H_zero) for v in orders)
                    high_count = sum(count for value, count in high.items() if value >= upper_threshold)
                    low_count = sum(count for value, count in low.items() if value <= lower_threshold)
                    assert high_count <= s**(-k)*U**length
                    assert low_count <= s**k*V**length
                    C_inf = 1+max(Fraction(0), max(degrees)+H_inf)
                    C_zero = 1+max(Fraction(0), degree+H_zero-min(orders))
                    assert ranks[first_cut][k]*ranks[1-first_cut][k] <= C_inf*C_zero*(length+1)**2*(U*V)**length
                    rational_checks += 3
    return {'name': name, 'block_toeplitz_checks': toeplitz_checks,
            'balanced_endpoint_entries': endpoint_entries,
            'shifted_ordered_product_entries': shift_entries,
            'coefficient_rank_checks': rank_checks,
            'exact_rational_bound_checks': rational_checks}


def controls():
    # A uniform one-power shift changes every Smith order; its mean correction
    # exactly cancels the apparent gain in the weighted envelope.
    q = 2; n = q*q; s = Fraction(5)
    old = n*n
    without_mean = n*(n/s)
    with_mean = s*without_mean
    assert with_mean == old and without_mean < old
    # Single-letter weights have different sums on balanced and unbalanced words.
    weights = ((0, 1), (0, 0), (0, 0), (0, 0))
    assert 4*mean(weights) == 2
    assert sum(weights[0][0] for _ in range(4)) == 0
    cardinalities = []
    import math
    for q in (2, 3, 4, 5):
        for multiple in range(1, 9):
            t = q*multiple
            size = math.factorial(t)//math.factorial(multiple)**q
            assert q**t <= size*(t+1)**q and size <= q**t
            cardinalities.append((q, t, size))
    return {'omitted_mean_changes_base_from': old, 'omitted_mean_false_base': str(without_mean),
            'unbalanced_word_shift': 0, 'balanced_word_shift': 2,
            'exact_balanced_cardinality_checks': len(cardinalities)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [('scalar elementary', *basic.elementary()),
             ('ordered auxiliary pairing', *basic.pairing((1, 6, 9), 0)),
             ('noncommuting auxiliary loop', *basic.loop(13, 18))]
    report = {'status': 'passed', 'field': 2,
              'cases': [check_family(name, local, degree) for name, local, degree in cases],
              'negative_controls': controls(),
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'base_smith_sha256': hashlib.sha256(Path(basic.__file__).read_bytes()).hexdigest(),
              'auxiliary_checker_sha256': hashlib.sha256(Path(auxiliary.__file__).read_bytes()).hexdigest(),
              'scope': 'Finite exact checks of the balanced-label transfer; no improving crossing certificate.'}
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
