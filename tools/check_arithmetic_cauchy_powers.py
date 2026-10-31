#!/usr/bin/env python3
"""Exact controls for every positive arithmetic Cauchy Hadamard power."""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb, factorial, prod
from pathlib import Path
from random import Random

from check_arithmetic_cauchy_factorial import local_smith_orders, valuation


@lru_cache(None)
def stirling(n, k):
    if n == k == 0:
        return 1
    if not 1 <= k <= n:
        return 0
    return stirling(n-1, k-1)+k*stirling(n-1, k)


def rising(m, r):
    return prod(m+j for j in range(r))


def predicted_orders(p, d, m, n):
    e = valuation(d, p)
    return [2*r*e+valuation(factorial(r), p)+valuation(rising(m, r), p)
            for r in range(n)]


def determinant(matrix):
    work = [list(map(Fraction, row)) for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next((row for row in range(column, len(work)) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            result = -result
        value = work[column][column]
        result *= value
        for row in range(column+1, len(work)):
            factor = work[row][column]/value
            for j in range(column+1, len(work)):
                work[row][j] -= factor*work[column][j]
            work[row][column] = 0
    return result


def at_least(value, p, lower):
    return not value or valuation(value, p) >= lower


def difference_controls():
    stirling_count = difference_count = 0
    for p in (2, 3, 5, 7):
        for k in range(35):
            for i in range(k+1):
                value = Fraction(factorial(i)*stirling(k, i), factorial(k))
                assert at_least(value, p, -Fraction(k-i, p-1))
                stirling_count += 1
        for delta in (Fraction(p), Fraction(p*p, 1+p)):
            for i, r in product(range(9), repeat=2):
                value = sum(((-1)**(i-j)*comb(i, j)*(Fraction(j)/(1+delta*j))**r
                             for j in range(i+1)), Fraction())/factorial(r)
                assert at_least(value, p, -Fraction(r-i, p-1))
                difference_count += 1
    return stirling_count, difference_count


def alternant_controls():
    rng = Random(310736)
    records = []
    for p, n in product((2, 3, 5), (2, 3, 4, 5)):
        delta = Fraction(p, 1+p)
        nodes = [Fraction(i)/(1+delta*i) for i in range(n)]
        vandermonde = prod(nodes[j]-nodes[i] for i in range(n) for j in range(i+1, n))
        options = list(combinations(range(1, n+5), n-1))
        rng.shuffle(options)
        for rest in options[:20]:
            exponents = (0,)+rest
            k = sum(r-i for i, r in enumerate(exponents))
            d = sum(valuation(factorial(r), p)-valuation(factorial(i), p)
                    for i, r in enumerate(exponents))
            schur = determinant([[node**r for r in exponents] for node in nodes])/vandermonde
            bound = max(0, d-k//(p-1))
            assert at_least(schur, p, bound)
            records.append({'prime': p, 'size': n, 'exponents': list(exponents),
                            'bound': bound, 'valuation': valuation(schur, p)})
    # Complete symmetric functions give single-row partitions at large n.
    # These force genuinely positive D-B bounds without a large determinant.
    large = []
    for p, power in product((2, 3, 5), (2, 3, 4)):
        n = p**power
        delta = Fraction(p)
        complete = [Fraction(1)]+[Fraction(0)]*4
        for i in range(n):
            node = Fraction(i)/(1+delta*i)
            for degree in range(1, 5):
                complete[degree] += node*complete[degree-1]
        for k in range(1, 5):
            d = valuation(prod(range(n, n+k)), p)
            bound = max(0, d-k//(p-1))
            assert at_least(complete[k], p, bound)
            large.append({'prime': p, 'size': n, 'partition': [k],
                          'bound': bound, 'valuation': valuation(complete[k], p)})
    return records, large


def determinant_controls():
    records = []
    for p, exponent, m, n in product((2, 3, 5), (1, 2), (1, 2, 3, 4, 5), (1, 2, 3, 4)):
        d, t = p**exponent, 1+p
        matrix = [[Fraction(1, (t+d*(i+j))**m) for j in range(n)] for i in range(n)]
        actual = determinant(matrix)
        leading = (Fraction(d**(n*(n-1)), t**(m*n+n*(n-1)))
                   *prod(factorial(r)*rising(m, r) for r in range(n)))
        quotient = actual/leading
        assert valuation(quotient, p) == 0 and at_least(quotient-1, p, 1)
        records.append({'prime': p, 'spacing_valuation': exponent, 'power': m,
                        'size': n, 'determinant_valuation': valuation(actual, p),
                        'normalized_residue': (quotient.numerator*pow(quotient.denominator, -1, p)) % p})
    return records


def smith_controls():
    records = []
    for p, e in product((2, 3, 5, 7), (1, 2)):
        d = p**e
        for m, n, t in product(sorted({1, 2, 3, 4, 5, 8, p**8-1, p**8, p**8+1}),
                               (2, 3, 5, 8, 13), (1, 1+p)):
            expected = predicted_orders(p, d, m, n)
            precision = max(expected)+24
            modulus = p**precision
            matrix = [[pow(pow(t+d*(i+j), -1, modulus), m, modulus)
                       for j in range(n)] for i in range(n)]
            actual = local_smith_orders(matrix, p, precision)
            assert actual == expected
            records.append({'prime': p, 'spacing_valuation': e, 'center': t,
                            'power': m, 'size': n, 'precision': precision,
                            'orders': actual})
    # A negative control outside p|d: the same formula need not hold.
    p, d, m, n = 5, 1, 4, 2
    modulus = p**24
    matrix = [[pow(pow(1+i+j, -1, modulus), m, modulus) for j in range(n)] for i in range(n)]
    actual = local_smith_orders(matrix, p, 24)
    assert actual == [0, 2] and actual != predicted_orders(p, d, m, n)
    return records, {'prime': p, 'spacing': d, 'power': m, 'size': n,
                     'actual_orders': actual, 'incorrect_extension': [0, 0]}


def orthogonal_controls():
    records = []
    for p, m in product((2, 3, 5), (3, 4, 5)):
        size, t, d = 4, 1+p, p
        moments = [Fraction(1, (t+d*k)**m) for k in range(2*size+1)]
        def pairing(a, b):
            return sum((x*y*moments[i+j] for i, x in enumerate(a)
                        for j, y in enumerate(b)), Fraction())
        basis, norms = [], []
        for r in range(size+1):
            polynomial = [Fraction(0)]*r+[Fraction(1)]
            for previous, norm in zip(basis, norms):
                coefficient = pairing(polynomial, previous)/norm
                for j, value in enumerate(previous):
                    polynomial[j] -= coefficient*value
            norm = pairing(polynomial, polynomial)
            assert valuation(polynomial[0], p) == 0
            assert all(at_least(value, p, 0) for value in polynomial)
            assert all(pairing(polynomial, previous) == 0 for previous in basis)
            assert valuation(norm, p) == predicted_orders(p, d, m, size+1)[r]
            basis.append(polynomial)
            norms.append(norm)
        records.append({'prime': p, 'power': m, 'maximum_degree': size,
                        'norm_orders': [valuation(norm, p) for norm in norms],
                        'integral_coefficients': sum(map(len, basis))})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    stirling_count, difference_count = difference_controls()
    alternants, large = alternant_controls()
    determinants = determinant_controls()
    profiles, negative = smith_controls()
    orthogonal = orthogonal_controls()
    gaps = []
    for p, a in product((2, 3, 5, 7), (1, 2, 4, 8)):
        size = min(p**a, 40)+1
        actual = predicted_orders(p, p, p**a, size)
        expected = [0]+[a+2*r+valuation(factorial(r), p)+valuation(factorial(r-1), p)
                        for r in range(1, size)]
        assert actual == expected
        gaps.append({'prime': p, 'power_exponent': a, 'size': size, 'orders': actual})
    root = Path(__file__).resolve().parent.parent
    sources = ('research/arithmetic-cauchy-powers.md', 'tools/check_arithmetic_cauchy_powers.py',
               'tools/check_arithmetic_cauchy_factorial.py')
    result = {'status': 'exact finite controls of an all-size ordinary proof; no new factoring exponent',
              'source_sha256': {path: sha256((root/path).read_bytes()).hexdigest() for path in sources},
              'stirling_denominator_checks': stirling_count,
              'finite_difference_checks': difference_count,
              'alternant_controls': alternants, 'large_node_alternants': large,
              'determinant_congruences': determinants, 'complete_smith_profiles': profiles,
              'integral_orthogonal_bases': orthogonal, 'prime_power_gaps': gaps,
              'spacing_hypothesis_negative_control': negative}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(json.dumps({'stirling': stirling_count, 'finite_differences': difference_count,
                      'alternants': len(alternants)+len(large), 'determinants': len(determinants),
                      'smith_profiles': len(profiles), 'orthogonal_bases': len(orthogonal),
                      'prime_power_gaps': len(gaps), 'negative_control': negative}, indent=2))


if __name__ == '__main__':
    main()
