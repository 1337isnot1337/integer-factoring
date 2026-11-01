#!/usr/bin/env python3
"""Exact controls for the weighted Toeplitz register-profile budget."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb, isqrt
from pathlib import Path
import random

from check_controlled_permutation_smith import binary_rank, phi_rows
from check_multipoint_hermite_crossings import rank as prime_rank
from check_shift_register_copy_wire_barrier import euler_register, encode_digits


def rank(matrix, prime):
    if prime == 2:
        return binary_rank([sum((value % 2) << j for j, value in enumerate(row))
                            for row in matrix])
    return prime_rank(matrix, prime)


def toeplitz(coefficients, length):
    rows, columns = len(coefficients[0]), len(coefficients[0][0])
    return [[coefficients[r - c][i][j] if c <= r else 0
             for c in range(length) for j in range(columns)]
            for r in range(length) for i in range(rows)]


def combination(coefficients, weights, prime):
    rows, columns = len(coefficients[0]), len(coefficients[0][0])
    return [[sum(weight * coefficients[k][i][j] for k, weight in enumerate(weights)) % prime
             for j in range(columns)] for i in range(rows)]


def projected_constant_map(matrix, weights, rows, columns, prime):
    # C selects the constant-input columns; P combines complete output blocks.
    return [[sum(weight * matrix[k * rows + i][j] for k, weight in enumerate(weights)) % prime
             for j in range(columns)] for i in range(rows)]


def known_smith_jets(rows, columns, orders, precision, prime, seed):
    rng = random.Random(seed)
    jets = [[[0] * columns for _ in range(rows)] for _ in range(precision)]
    for i, order in enumerate(orders):
        if order < precision:
            jets[order][i][i] = 1
    # Polynomial elementary shears are unimodular, including those at degree zero.
    for _ in range(16):
        side = rng.randrange(2)
        size = rows if side == 0 else columns
        target, source = rng.sample(range(size), 2)
        degree, coefficient = rng.randrange(precision), rng.randrange(1, prime)
        for k in range(degree, precision):
            if side == 0:
                for j in range(columns):
                    jets[k][target][j] = (jets[k][target][j]
                        + coefficient * jets[k - degree][source][j]) % prime
            else:
                for i in range(rows):
                    jets[k][i][target] = (jets[k][i][target]
                        + coefficient * jets[k - degree][i][source]) % prime
    return jets


def arbitrary_matrix_controls(counts):
    precision = 7
    for prime in (2, 3, 5):
        for rows, columns in ((2, 3), (3, 2), (3, 3), (4, 4)):
            for seed in range(4):
                rng = random.Random(1000 * prime + 100 * rows + 10 * columns + seed)
                orders = sorted(rng.randrange(precision + 2)
                                for _ in range(min(rows, columns)))
                jets = known_smith_jets(rows, columns, orders, precision, prime, seed)
                for length in range(1, precision + 1):
                    matrix = toeplitz(jets, length)
                    expected = sum(max(0, length - order) for order in orders)
                    actual = rank(matrix, prime)
                    assert actual == expected
                    counts['known_smith_toeplitz_ranks'] += 1
                    assert expected <= length * sum(order < length for order in orders)
                    for weights in ([rng.randrange(prime) for _ in range(length)],
                                    [(-1) ** k * comb(length - 1, k) % prime
                                     for k in range(length)]):
                        direct = combination(jets, weights, prime)
                        projected = projected_constant_map(matrix, weights, rows, columns, prime)
                        assert direct == projected
                        assert rank(direct, prime) <= expected
                        counts['coefficient_projection_rank_cases'] += 1


def register_controls(counts):
    records = []
    for radix, digits, prime in ((2, 2, 2), (2, 4, 2), (2, 6, 2),
                                 (3, 2, 3), (3, 4, 3), (4, 2, 2)):
        transition = euler_register(radix, digits)
        partitions = list(combinations(range(digits), digits // 2))
        if digits >= 6 or radix >= 3:
            partitions = [tuple(range(digits // 2))]
        else:
            partitions = [first for first in partitions if 0 in first]
        for first in partitions:
            second = tuple(i for i in range(digits) if i not in first)
            n, q = radix ** digits, radix ** (digits // 2)
            encoding = [encode_digits(s, first, second, radix) for s in range(n)]
            powers, canonical = [], list(range(n))
            for exponent in range(1, digits):
                canonical = [transition[state] for state in canonical]
                power = [0] * n
                for state, target in enumerate(canonical):
                    power[encoding[state]] = encoding[target]
                rows = phi_rows(power, q)
                powers.append([[(row >> column) & 1 for column in range(n)] for row in rows])
            jets = []
            for k in range(digits - 1):
                jets.append([[sum((-1) ** ell * comb(k, ell) * powers[ell][i][j]
                                  for ell in range(k + 1)) % prime
                              for j in range(n)] for i in range(n)])
            profile_ranks = [0]
            for length in range(1, digits):
                weights = [(-1) ** k * comb(length - 1, k) % prime for k in range(length)]
                direct = combination(jets, weights, prime)
                assert direct == powers[length - 1]
                matrix = toeplitz(jets, length)
                projected = projected_constant_map(matrix, weights, n, n, prime)
                assert direct == projected
                budget, power_rank = rank(matrix, prime), rank(direct, prime)
                profile_ranks.append(budget)
                early = profile_ranks[-1] - profile_ranks[-2]
                assert power_rank <= budget <= length * early
                assert length * length <= n
                records.append({'radix': radix, 'digits': digits, 'prime': prime,
                                'first': list(first), 'power': length,
                                'power_rank': power_rank, 'toeplitz_budget': budget,
                                'early_channels': early})
                counts['complete_register_coefficient_cases'] += 1
    return records


def scalar_controls(counts):
    scales = tuple(map(Fraction, ('1', '65/64', '9/8', '3/2', '2', '8')))
    for n in (4, 9, 16, 64, 81, 256):
        for length in range(1, isqrt(n) + 1):
            for order in range(length):
                assert order * (order + 1) <= order * length
                assert order * length * length <= order * n
                for scale in scales:
                    value = scale ** (n - 1 - order) + (n - 1) * scale ** (-1 - order)
                    assert length * value >= n * (length - order)
                    counts['exact_rational_channel_cases'] += 1
    for digits in range(2, 1001):
        assert (digits - 1) ** 2 <= 1 << digits
        counts['register_window_size_controls'] += 1


def witness(digits):
    minimum_sum = isqrt(2 * digits)
    if minimum_sum * minimum_sum < 2 * digits:
        minimum_sum += 1
    length = next(j for j in range(1, digits)
                  if j + (digits + 2 * j - 1) // (2 * j) == minimum_sum)
    exponent = digits // 2 - minimum_sum
    return length, exponent


def cutoff_controls(counts):
    records = []
    for digits in range(380, 576, 2):
        length, exponent = witness(digits)
        power = 12 * exponent - 5 * digits
        assert power >= 0 and (1 << power) >= length ** 12
        records.append({'digits': digits, 'J': length, 'E': exponent,
                        'binary_exponent': power, 'right_side': str(length ** 12)})
    assert len(records) == 98
    assert (records[0]['J'], records[0]['E'], records[0]['binary_exponent']) == (12, 162, 44)
    assert 1 << 44 == 17592186044416
    assert 12 ** 12 == 8916100448256
    for length in range(1, 378):
        exponent = 189 - (378 + 2 * length - 1) // (2 * length) - length
        power = 12 * exponent - 5 * 378
        assert power < 0 or (1 << power) < length ** 12
    counts['finite_cutoff_witnesses'] = len(records)
    counts['last_failure_all_window_checks'] = 377
    tail = []
    for digits in (576, 578, 1024, 1296, 2048, 4096, 16384):
        length = isqrt(digits)
        exponent = digits // 2 - (digits + 2 * length - 1) // (2 * length) - length
        assert (1 << (12 * exponent - 5 * digits)) >= length ** 12
        tail.append({'digits': digits, 'J': length, 'E': exponent})
    assert 2 * 17 ** 2 > 24 ** 2
    return {'finite_witnesses': records, 'analytic_tail_from': 576,
            'tail_controls': tail, 'last_failed_integer_certificate': 378}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    counts = Counter()
    arbitrary_matrix_controls(counts)
    registers = register_controls(counts)
    scalar_controls(counts)
    cutoff = cutoff_controls(counts)
    result = {'statement': 'B(s) >= d^(m+E)/J for every full-cycle reversible shift register.',
              'explicit_even_digit_cutoff': 380, 'bare_rate_lower_bound': '17/24',
              'asymptotic_bare_rate_lower_bound': '3/4', 'counts': dict(counts),
              'register_controls': registers, 'cutoff': cutoff,
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope': 'Bare two-place Smith profiles and balanced coordinate partitions; factoring exponent unchanged.'}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'cutoff': 380, **dict(counts)}, sort_keys=True))


if __name__ == '__main__':
    main()
