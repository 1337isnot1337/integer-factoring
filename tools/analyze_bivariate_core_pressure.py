#!/usr/bin/env python3
"""Fit frozen core profiles and check the bivariate valuation lemma.

The profile fit is finite evidence. It is not an all-M profile theorem or
a strict pressure certificate. All algebraic controls use exact integers.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
import math
from pathlib import Path
import random


def add(left, right):
    """Sparse Laurent polynomials over F2, stored as exponent-pair sets."""
    return left ^ right


def multiply(left, right):
    result = set()
    for a, b in left:
        for c, d in right:
            pair = (a+c, b+d)
            if pair in result:
                result.remove(pair)
            else:
                result.add(pair)
    return result


def determinant(matrix):
    n = len(matrix)
    result = set()
    for permutation in permutations(range(n)):
        term = {(0, 0)}
        for row, column in enumerate(permutation):
            term = multiply(term, matrix[row][column])
        result = add(result, term)
    return result


def lex_order(polynomial):
    return min((b, a) for a, b in polynomial)


def specialization_order(polynomial, M):
    coefficients = set()
    for a, b in polynomial:
        exponent = a+b*M
        if exponent in coefficients:
            coefficients.remove(exponent)
        else:
            coefficients.add(exponent)
    return None if not coefficients else min(coefficients)


def width(polynomial):
    return max(a for a, _ in polynomial)-min(a for a, _ in polynomial)


def minors(matrix, size):
    n = len(matrix)
    return [determinant([[matrix[i][j] for j in columns] for i in rows])
            for rows in combinations(range(n), size)
            for columns in combinations(range(n), size)]


def algebra_controls():
    rng = random.Random(718261)
    scalar_count = matrix_minor_count = 0
    for _ in range(100):
        polynomial = {(rng.randrange(-5, 6), rng.randrange(-3, 4))
                      for _ in range(12)}
        slope, offset = lex_order(polynomial)
        for M in (width(polynomial)+1, width(polynomial)+7):
            assert specialization_order(polynomial, M) == slope*M+offset
            scalar_count += 1
    # At a small exceptional M, distinct exponent pairs can cancel.
    cancellation = {(1, 0), (0, 1)}
    assert specialization_order(cancellation, 1) is None
    assert specialization_order(cancellation, 2) == 1
    denominator = {(0, 0), (1, 0), (0, 1)}
    matrices = []
    for _ in range(16):
        matrix = [[{(rng.randrange(-2, 4), rng.randrange(-1, 3))
                    for _ in range(3)} for _ in range(3)] for _ in range(3)]
        all_terms = set().union(*(entry for row in matrix for entry in row))
        wp = width(all_terms)
        by_size = [minors(matrix, k) for k in range(1, 4)]
        rank = max(k for k, entries in enumerate(by_size, 1) if any(entries))
        bound = max(rank*wp, width(denominator))
        pair_q = lex_order(denominator)
        prefixes = [(0, 0)]
        for k, entries in enumerate(by_size[:rank], 1):
            pair = min(lex_order(entry) for entry in entries if entry)
            prefixes.append((pair[0]-k*pair_q[0], pair[1]-k*pair_q[1]))
        pairs = [(b[0]-a[0], b[1]-a[1])
                 for a, b in zip(prefixes, prefixes[1:])]
        for M in (bound+1, bound+11):
            q_order = specialization_order(denominator, M)
            previous = 0
            actual = []
            for k, entries in enumerate(by_size[:rank], 1):
                values = [specialization_order(entry, M)
                          for entry in entries if entry]
                assert all(value is not None for value in values)
                matrix_minor_count += len(values)
                prefix = min(values)-k*q_order
                actual.append(prefix-previous)
                previous = prefix
            assert actual == [a*M+b for a, b in pairs]
        matrices.append({'rank': rank, 'sufficient_M_strictly_above': bound,
                         'smith_affine_pairs': pairs})
    return {'scalar_specializations': scalar_count,
            'matrix_minor_specializations': matrix_minor_count,
            'rational_matrix_cases': matrices,
            'small_M_cancellation_control': {'polynomial': 'z+Z',
                'M1_specialization': 'zero', 'M2_valuation': 1}}


def profile_fit(records, place):
    cases = {entry['M']: entry for entry in records if entry['place'] == place}
    first, second = cases[16]['orders'], cases[32]['orders']
    assert len(first) == len(second)
    pairs = []
    for a, b in zip(first, second):
        assert (b-a) % 16 == 0
        slope = (b-a)//16
        pairs.append((slope, a-16*slope))
    checked = []
    for M, entry in sorted(cases.items()):
        predicted = sorted(a*M+b for a, b in pairs)
        assert predicted == entry['orders']
        checked.append({'M': M, 'orders_checked': len(predicted),
                        'source_sha256': entry['source_sha256']})
    multiplicities = Counter(a for a, _ in pairs)
    smallest = min(multiplicities)
    return {'place': place, 'fit_sizes': [16, 32], 'checked_sizes': checked,
            'affine_pairs': pairs,
            'slope_multiplicities': [[a, multiplicities[a]]
                                     for a in sorted(multiplicities)],
            'smallest_slope': smallest,
            'smallest_slope_offsets': sorted(b for a, b in pairs if a == smallest),
            'smallest_slope_multiplicity': multiplicities[smallest]}


def gaussian_leading(place):
    # Active one-particle orders are a*M+b. The complete profile contains
    # every subset sum. This calculation is independent of the fitted core.
    pairs = ([(0, 1), (0, 3), (0, 5), (2, -5), (2, -3), (2, -1),
              (0, 1), (2, -3), (2, -1)] if place == 'zero' else
             [(-2, 1), (-2, 3), (0, -1), (0, 1), (0, 3), (2, -3), (2, -1)])
    smallest = sum(min(a, 0) for a, _ in pairs)
    multiplicity = 2**sum(a == 0 for a, _ in pairs)
    return {'place': place, 'one_particle_affine_pairs': pairs,
            'smallest_slope': smallest, 'multiplicity': multiplicity}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=root/'evidence/bivariate-core-pressure-inputs.json')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    inputs = json.loads(args.input.read_text())
    fits = [profile_fit(inputs['records'], place) for place in ('zero', 'infinity')]
    gaussian = [gaussian_leading(place) for place in ('zero', 'infinity')]
    numerator = denominator = 1
    for fit, control in zip(fits, gaussian):
        assert fit['smallest_slope'] == control['smallest_slope']
        numerator *= fit['smallest_slope_multiplicity']
        denominator *= control['multiplicity']
    divisor = math.gcd(numerator, denominator)
    numerator //= divisor
    denominator //= divisor
    result = {
        'status': 'Exact finite affine-fit controls; the extrapolated pressure limit is conditional.',
        'profiles': fits, 'gaussian_controls': gaussian,
        'conditional_natural_limit': {
            'physical_core_minimum_terms_cancel': True,
            'leading_tau_coefficient': 0,
            'excess_exponential_ratio': [numerator, denominator],
            'logarithmic_excess_approximation': math.log(numerator/denominator)},
        'algebra_controls': algebra_controls(),
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'tool_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(json.dumps({'fit_counts': {entry['place']: len(entry['affine_pairs'])
                                     for entry in fits},
                      'conditional_limit_ratio': [numerator, denominator],
                      'algebra_controls_passed': True}))


if __name__ == '__main__':
    main()
