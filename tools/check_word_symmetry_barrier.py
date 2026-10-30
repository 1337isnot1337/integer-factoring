#!/usr/bin/env python3
"""Exact hook arithmetic and constants for the characteristic-zero word barrier."""

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import factorial, prod
from pathlib import Path


@lru_cache(None)
def partitions(total, length, ceiling):
    if not length:
        return ((),) if total == 0 else ()
    return tuple((first,)+tail for first in range(min(total, ceiling), -1, -1)
                 for tail in partitions(total-first, length-1, first))


def specht_dimension(shape):
    n, q = sum(shape), len(shape)
    numerator = factorial(n)*prod(shape[i]-shape[j]+j-i
                                 for i in range(q) for j in range(i+1, q))
    denominator = prod(factorial(shape[i]+q-i-1) for i in range(q))
    assert numerator % denominator == 0
    return numerator//denominator


def hook_dimension(shape):
    hooks = [shape[i]-j+sum(row > j for row in shape[i+1:])
             for i in range(len(shape)) for j in range(shape[i])]
    numerator, denominator = factorial(sum(shape)), prod(hooks)
    assert numerator % denominator == 0
    return numerator//denominator


def log_bounds(value, terms=32):
    exponent = value.bit_length()-1
    def series(x):
        t = (x-1)/(x+1)
        lower = 2*sum(t**(2*j+1)/Fraction(2*j+1) for j in range(terms))
        error = 2*t**(2*terms+1)/(Fraction(2*terms+1)*(1-t*t))
        return lower, lower+error
    l2, u2 = series(Fraction(2))
    lo, hi = series(Fraction(value, 1 << exponent))
    return exponent*l2+lo, exponent*u2+hi


def exact_rank(matrix):
    a = [[Fraction(v) for v in row] for row in matrix]
    rank = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][j]
        a[rank] = [v/scale for v in a[rank]]
        for i in range(len(a)):
            if i != rank and a[i][j]:
                scale = a[i][j]
                a[i] = [x-scale*y for x, y in zip(a[i], a[rank])]
        rank += 1
    return rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hook_cases = 0
    minima = []
    for q, maximum_k in [(4, 8), (9, 3), (16, 2)]:
        for k in range(1, maximum_k+1):
            n = 2*q*k
            bound_denominator = (factorial((q+1)*k)*factorial(k)**(q-1)
                                 *(n+q)**(q*(q-1)//2))
            smallest = None
            for remainder in partitions(q*k, q, q*k):
                shape = tuple(v+k for v in remainder)
                dimension = specht_dimension(shape)
                assert dimension == hook_dimension(shape)
                assert dimension*bound_denominator >= factorial(n)
                smallest = dimension if smallest is None else min(smallest, dimension)
                hook_cases += 1
            minima.append({'local_matrix_size': q, 'k': k, 'word_length': n,
                           'minimum_specht_dimension_in_checked_region': smallest})
    intervals = []
    for d in range(2, 6):
        q = d*d
        a, b = log_bounds(2*q)
        c, e = log_bounds(q+1)
        f, g = log_bounds(q)
        coefficient = Fraction(q+1, 2*q)
        lower, upper = (a-coefficient*e)/g, (b-coefficient*c)/f
        if d <= 4:
            assert lower > 0 and 2*lower*lower > 1
        else:
            assert 2*upper*upper < 1  # This bound alone does not exclude seed five.
        scale = 10**12
        lo = lower.numerator*scale//lower.denominator
        hi = -((-upper.numerator*scale)//upper.denominator)
        intervals.append({'seed_dimension': d, 'rho_lower_decimal': f'{lo//scale}.{lo%scale:012d}',
                          'rho_upper_decimal': f'{hi//scale}.{hi%scale:012d}',
                          'strictly_above_half_threshold': d <= 4})
    # A small exact Schur-block control: det(A)^2 at local matrix size two.
    # The invariant matrix 2I-A_Johnson+2C on weight-two binary words has rank two.
    matrix = [[0]*16 for _ in range(16)]
    for x in range(16):
        for y in range(16):
            if x.bit_count() == y.bit_count() == 2:
                matrix[x][y] = 2 if (x == y or x ^ y == 15) else -1
    functional = matrix[3][3]-2*matrix[3][5]+matrix[3][12]
    assert functional == 6 and exact_rank(matrix) == specht_dimension((2, 2)) == 2
    result = {'status': 'exact arithmetic for a characteristic-zero invariant-word obstruction',
              'hook_formula_cases': hook_cases, 'checked_region_minima': minima,
              'rho_intervals': intervals, 'small_schur_block_rank': 2,
              'small_functional_value': functional,
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'checked_region_minima'}))


if __name__ == '__main__':
    main()
