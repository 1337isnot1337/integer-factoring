#!/usr/bin/env python3
"""Exact minor, convexity, and entropy-flag controls for Smith pressure."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from random import Random

from weighted_smith_oracle import WeightedSmithOracle, polynomial_multiply


def multiply(first, second):
    n = len(first)
    answer = [[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            for j in range(n):
                answer[i][j] ^= polynomial_multiply(first[i][k], second[k][j])
    return answer


def partition(orders):
    # All controls use half-integral energies: 4^-kappa is rational.
    assert all((2*order).denominator == 1 for order in orders)
    return sum(Fraction(2)**(-int(2*order)) for order in orders)


def profile(oracle, parameters):
    n = oracle.size
    return oracle.profile(parameters[1:n+1], parameters[n+1:], parameters[0])


def flag_prefix(minor, parameters, n):
    return (parameters[0]*minor.valuation
            + sum(parameters[1+i] for i in range(n) if minor.row_mask >> i & 1)
            + sum(parameters[1+n+i] for i in range(n) if minor.column_mask >> i & 1))


def check():
    rng = Random(562901)
    comparisons = flags = prefix_tests = entropy_tests = 0
    for n in (2, 3, 4):
        for instance in range(4):
            lower = [[int(i == j) for j in range(n)] for i in range(n)]
            upper = [[int(i == j) for j in range(n)] for i in range(n)]
            for i in range(n):
                for j in range(i):
                    lower[i][j] = rng.randrange(1, 32)
                    upper[j][i] = rng.randrange(1, 32)
            diagonal = [[(1 << rng.randrange(5)) if i == j else 0
                         for j in range(n)] for i in range(n)]
            matrix = multiply(multiply(lower, diagonal), upper)
            for infinity in (False, True):
                oracle = WeightedSmithOracle(matrix, infinity=infinity)
                fixtures = [[Fraction(rng.randrange(1, 4))]
                            + [Fraction(rng.randrange(-3, 4)) for _ in range(2*n)]
                            for _ in range(5)]
                for left in fixtures:
                    left_orders, left_flag = profile(oracle, left)
                    weights = [Fraction(2)**(-int(2*order)) for order in left_orders]
                    total = sum(weights)
                    probabilities = [weight/total for weight in weights]
                    assert probabilities == sorted(probabilities, reverse=True)
                    for right in fixtures:
                        right_orders, _ = profile(oracle, right)
                        midpoint = [(a+b)/2 for a, b in zip(left, right)]
                        middle_orders, _ = profile(oracle, midpoint)
                        assert partition(middle_orders)**2 <= partition(left_orders)*partition(right_orders)
                        comparisons += 1
                        previous = 0
                        increments = []
                        for k, minor in enumerate(left_flag, 1):
                            value = flag_prefix(minor, right, n)
                            assert sum(right_orders[:k]) <= value
                            increments.append(value-previous)
                            previous = value
                            prefix_tests += 1
                            # Uniform-prefix entropy supports: U^k 4^(S_k) >= k^k.
                            assert partition(right_orders)**k * Fraction(2)**int(2*value) >= k**k
                            entropy_tests += 1
                        assert sum(p*v for p, v in zip(probabilities, right_orders)) <= sum(
                            p*v for p, v in zip(probabilities, increments))
                        flags += 1
    first, second = [[0]*4 for _ in range(4)], [[0]*4 for _ in range(4)]
    for x in range(2):
        for y in range(2):
            for u in range(2):
                for v in range(2):
                    identity = int(x == u and y == v)
                    crossing = int(x == v and y == u)
                    value = identity | ((identity ^ crossing) << 1)
                    first[2*x+y][2*u+v] = value
                    second[2*x+u][2*y+v] = value
                    assert value.bit_count() % 2 == crossing
    row, _ = WeightedSmithOracle(first, infinity=True).profile([0]*4, [0]*4)
    column, _ = WeightedSmithOracle(second).profile([0]*4, [0]*4)
    assert row == (-1, 0, 0, 1) and column == (0, 1, 1, 2)
    probabilities = tuple(map(Fraction, ('9/16', '3/16', '3/16', '1/16')))
    assert sum(probabilities) == 1
    assert sum(p*(a+b) for p, a, b in zip(probabilities, row, column)) == 0
    entropy_product = Fraction(1)
    for p in probabilities:
        entropy_product *= p**(-int(32*p))
    assert entropy_product == Fraction(256, 27)**16
    assert sum(Fraction(3)**(-a) for a in row)*sum(Fraction(3)**(-b) for b in column) == Fraction(256, 27)
    return dict(exact_log_convexity_comparisons=comparisons,
                global_weighted_flag_comparisons=flags, prefix_minor_comparisons=prefix_tests,
                rational_entropy_supports=entropy_tests,
                sharp_elementary_control=dict(row_infinity_orders=list(row), column_zero_orders=list(column),
                                              flag_probabilities=list(map(str, probabilities)),
                                              exact_global_product_floor='256/27', equality_scale=3),
                scope='Square full-rank polynomial matrices; positive valuation scale and arbitrary diagonal physical-energy maps.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
