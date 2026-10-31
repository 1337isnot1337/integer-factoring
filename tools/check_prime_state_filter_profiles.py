#!/usr/bin/env python3
"""Exact checks of the prime-state finite-difference profile theorem.

The polynomial factor identities are checked before any monomial
specialization. Independent active-matrix Smith calculations include
equal exponents and collisions between different ramification bands.
Full physical n=1 tensors are checked by default; --extended also
constructs both 625-by-625 cuts of one characteristic-five n=2 tensor.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from math import comb, factorial
from pathlib import Path
import time

from check_all_field_pascal_crossings import (
    add, multiply, scale, pascal_kernel, ramify_negative, smith_orders,
    weighted_orders,
)
from check_truncated_symmetric_crossings import direct_kernel, basis
from prime_state_filter_profiles import active_profiles


def polynomial_add(first, second, prime, coefficient=1):
    result = dict(first)
    for degree, value in second.items():
        result[degree] = (result.get(degree, 0) + coefficient * value) % prime
    return {degree: value for degree, value in result.items() if value}


def polynomial_multiply(first, second, prime):
    result = {}
    for degree, a in first.items():
        for other, b in second.items():
            result[degree + other] = (result.get(degree + other, 0) + a * b) % prime
    return {degree: value for degree, value in result.items() if value}


def polynomial_scale(value, coefficient, prime):
    return {degree: a * coefficient % prime for degree, a in value.items()
            if a * coefficient % prime}


def c_power(exponent, prime):
    return {i: (-1)**i * comb(exponent, i) % prime
            for i in range(exponent + 1) if comb(exponent, i) % prime}


def algebra_multiply(first, second, prime):
    result = {}
    for powers, a in first.items():
        for other, b in second.items():
            target = tuple(i + j for i, j in zip(powers, other))
            if max(target) >= prime:
                continue
            result[target] = polynomial_add(
                result.get(target, {}), polynomial_multiply(a, b, prime), prime)
    return {powers: value for powers, value in result.items() if value}


def universal_q(prime):
    # Variable order is V,Y,U,X. The coefficient variable is u; c=1-u.
    return {
        (0, 1, 1, 0): {0: 1},
        (0, 1, 0, 1): {0: prime - 1},
        (1, 0, 1, 0): {0: 1, 1: prime - 1},
        (1, 0, 0, 1): {0: prime - 1, 1: 1},
    }


def exact_factor_control(prime):
    """Compare the displayed L/R factors with the literal coefficient algebra."""
    q = universal_q(prime)
    exponential_vx = {(s, 0, 0, s): {0: pow(factorial(s), -1, prime)}
                      for s in range(prime)}
    powers = [exponential_vx]
    for _ in range(1, prime):
        powers.append(algebra_multiply(powers[-1], q, prime))
    checks = 0
    blocks = []
    for charge in range(1 - prime, prime):
        rows = range(max(0, charge), min(prime, prime + charge))
        modes = range(max(0, -charge), min(prime, prime - charge))
        left_antidiagonal, right_antidiagonal = [], []
        for i, j, k in product(rows, rows, modes):
            actual = {}
            for d in range(prime - k):
                coefficient = (-1)**k * comb(prime - 1 - d, k)
                coefficient *= pow(factorial(d), -1, prime)
                inverse = polynomial_scale(c_power(prime - 1 - d - k, prime),
                                           coefficient, prime)
                term = polynomial_multiply(
                    powers[d].get((j - charge, i, j, i - charge), {}),
                    inverse, prime)
                actual = polynomial_add(actual, term, prime)
            if i + k >= prime or j + k >= prime:
                expected = {}
            else:
                numerator = (-1)**(i - charge) * factorial(i + k) * factorial(j + k)
                denominator = (factorial(i) * factorial(i - charge) * factorial(k)
                               * factorial(k + charge) * factorial(j)
                               * factorial(j - charge))
                expected = polynomial_scale(c_power(prime - 1 - i - k, prime),
                                            numerator * pow(denominator, -1, prime), prime)
            assert actual == expected, (prime, charge, i, j, k, actual, expected)
            checks += 1
        for i in rows:
            k = prime - 1 - i
            left = (-1)**(i - charge) * factorial(prime - 1)
            left *= pow(factorial(i) * factorial(i - charge), -1, prime)
            right = factorial(prime - 1) * pow(
                factorial(k) * factorial(k + charge)
                * factorial(i) * factorial(i - charge), -1, prime)
            assert left % prime and right % prime
            left_antidiagonal.append(left % prime)
            right_antidiagonal.append(right % prime)
        blocks.append({'charge': charge, 'size': len(rows),
                       'beta_indices': list(modes),
                       'left_antidiagonal': left_antidiagonal,
                       'right_antidiagonal': right_antidiagonal})
    return {'prime': prime, 'polynomial_coefficient_identities': checks, 'blocks': blocks}


@lru_cache(maxsize=None)
def universal_alpha_blocks(prime):
    """Build the original active tensor, without the claimed diagonalization."""
    gaussian = {
        (s, r, r, s): {s: pow(factorial(r) * factorial(s), -1, prime)}
        for r, s in product(range(prime), repeat=2)
    }
    powers = [gaussian]
    q = universal_q(prime)
    for _ in range(1, prime):
        powers.append(algebra_multiply(powers[-1], q, prime))
    coefficients = []
    for digit in range(prime):
        tensor = {}
        for j in range(digit, prime):
            scalar = (-1)**(j - digit) * pow(factorial(digit) * factorial(j - digit), -1, prime)
            for powers_key, value in powers[j].items():
                tensor[powers_key] = polynomial_add(
                    tensor.get(powers_key, {}), value, prime, scalar)
        coefficients.append(tensor)
    result = []
    for charge in range(1 - prime, prime):
        rows = [(x, y) for x, y in product(range(prime), repeat=2) if y - x == charge]
        columns = [(v, u) for v, u in product(range(prime), repeat=2) if u - v == charge]
        matrices = []
        for digit in range(prime):
            matrices.append([
                [polynomial_scale(coefficients[digit].get((v, y, u, x), {}),
                                  factorial(x) * factorial(u), prime)
                 for v, u in columns]
                for x, y in rows])
        result.append((charge, tuple(matrices)))
    return tuple(result)


def direct_active_orders(prime, n, ramification, charge, exponents):
    actual = []
    for block_charge, matrices in universal_alpha_blocks(prime):
        size = len(matrices[0])
        sparse = [[{} for _ in range(size)] for _ in range(size)]
        for digit, matrix in enumerate(matrices):
            for i in range(size):
                for j in range(size):
                    for degree, value in matrix[i][j].items():
                        order = 2 * n * ramification * degree - 2 * exponents[digit]
                        sparse[i][j][order] = (sparse[i][j].get(order, 0) + value) % prime
        sparse = [[{degree: value for degree, value in entry.items() if value}
                   for entry in row] for row in sparse]
        minimum = min(degree for row in sparse for entry in row for degree in entry)
        maximum = max(degree for row in sparse for entry in row for degree in entry)
        dense = [[tuple(entry.get(degree, 0) for degree in range(minimum, max(entry) + 1))
                  if entry else () for entry in row] for row in sparse]
        bound = max(512, 1 << (size * (maximum - minimum + 1)).bit_length())
        orders, precision = smith_orders(dense, prime, maximum_precision=bound)
        actual.extend(order + minimum + ramification * (1 - charge) * block_charge
                      for order in orders)
    return tuple(sorted(actual))


def active_controls():
    result = []
    for prime, n, ramification in ((3, 2, 1), (5, 1, 1), (5, 3, 3), (7, 1, 1)):
        patterns = [(0,) * prime,
                    tuple(0 if d == 0 else 1 for d in range(prime)),
                    tuple(0 if d == 0 else -1 for d in range(prime)),
                    tuple((d % 3 - 1) * n * ramification for d in range(prime))]
        if prime == 5:
            patterns.extend(((0, 1, 2, 2, 1), (0, 3, 3, 1, 1)))
        for exponents in patterns:
            charge = (len(result) % 5) - 1
            expected = active_profiles(prime, n, ramification, charge, exponents,
                                       include_spectators=False)
            actual = direct_active_orders(prime, n, ramification, charge, exponents)
            assert actual == expected['zero_active'], (prime, n, ramification, charge, exponents)
            result.append({'prime': prime, 'n': n, 'ramification': ramification,
                           'charge': charge, 'exponents': exponents,
                           'beta_orders': expected['beta_orders'], 'zero_active': actual})
    return result


def tuple_algebra_multiply(first, second, prime):
    result = {}
    for powers, a in first.items():
        for other, b in second.items():
            target = tuple(i + j for i, j in zip(powers, other))
            if max(target) >= prime:
                continue
            result[target] = add(result.get(target, ()), multiply(a, b, prime), prime)
    return {powers: value for powers, value in result.items() if value}


def physical_filtered_tensor(prime, n, ramification, exponents):
    """Construct the original polynomial tensor before any active reduction."""
    one_particle = pascal_kernel(range(n), 2 * n, prime)
    gaussian = direct_kernel(one_particle, prime)
    left, right = {}, {}
    for position, coefficient in ((0, (1, 1)), (n, (1,))):
        powers = [0] * (4 * n)
        powers[position] = 1
        left[tuple(powers)] = coefficient
    for position, coefficient in [(2 * n, (1,))] + [
            (3 * n + b, (0,) * (2 * n - 1 - b) + ((-1)**b % prime,))
            for b in range(n)]:
        powers = [0] * (4 * n)
        powers[position] = 1
        right[tuple(powers)] = coefficient
    q = tuple_algebra_multiply(left, right, prime)
    powers = [gaussian]
    for _ in range(1, prime):
        powers.append(tuple_algebra_multiply(powers[-1], q, prime))
    clearing = max(exponents)
    result = {}
    for j, summand in enumerate(powers):
        multiplier = [0] * (clearing - min(exponents) + 1)
        for digit in range(j + 1):
            shift = clearing - exponents[digit]
            scalar = (-1)**(j - digit + shift) * pow(
                factorial(digit) * factorial(j - digit), -1, prime)
            multiplier[shift] = (multiplier[shift] + scalar) % prime
        for occupations, value in summand.items():
            ramified = [0] * ((len(value) - 1) * ramification + 1)
            for degree, coefficient in enumerate(value):
                ramified[ramification * degree] = coefficient
            result[occupations] = add(result.get(occupations, ()),
                                      multiply(tuple(ramified), tuple(multiplier), prime), prime)
    result = {powers: value for powers, value in result.items() if value}
    endpoint = lambda value: sum(coefficient * (-1)**degree
                                for degree, coefficient in enumerate(value)) % prime
    support = set(result) | set(gaussian)
    assert all(endpoint(result.get(powers, ())) == endpoint(gaussian.get(powers, ()))
               for powers in support)
    return result, clearing, len(support)


def physical_cut_profiles(tensor, prime, n, ramification, charge):
    digits = basis(n, prime)
    lookup = {word: i for i, word in enumerate(digits)}
    dimension = prime**n
    heights = ([0] * n, [2 * i + 1 - charge for i in range(n)], [0] * n,
               [2 * i - 4 * n + 1 + charge for i in range(n)],
               [2 * i - n + 2 for i in range(n)], [-n] * n,
               [2 * i - n + 1 for i in range(n)], [3 * n - 1] * n)
    heights = [[ramification * value for value in row] for row in heights]
    entries = {}
    for powers, value in tensor.items():
        labels = tuple(lookup[powers[i * n:(i + 1) * n]] for i in range(4))
        coefficient = 1
        for digit in powers[2 * n:]:
            coefficient *= factorial(digit)
        entries[labels] = ramify_negative(scale(value, coefficient, prime), prime)
    result = []
    for place, ordering in enumerate(((3, 1, 0, 2), (3, 2, 0, 1))):
        size = dimension**2
        matrix = [[()] * size for _ in range(size)]
        for labels, value in entries.items():
            row = labels[ordering[0]] * dimension + labels[ordering[1]]
            column = labels[ordering[2]] * dimension + labels[ordering[3]]
            matrix[row][column] = value
        rows, columns = [], []
        for first, second in product(range(dimension), repeat=2):
            rows.append(sum(heights[4 * place + ordering[k]][i] * digits[letter][i]
                            for k, letter in enumerate((first, second)) for i in range(n)))
            columns.append(sum(heights[4 * place + ordering[k + 2]][i] * digits[letter][i]
                               for k, letter in enumerate((first, second)) for i in range(n)))
        adjacency = [[] for _ in range(2 * size)]
        for i in range(size):
            for j in range(size):
                if matrix[i][j]:
                    adjacency[i].append(size + j)
                    adjacency[size + j].append(i)
        seen, orders = set(), []
        for start in range(size):
            if start in seen:
                continue
            queue, block_rows, block_columns = [start], [], []
            seen.add(start)
            while queue:
                vertex = queue.pop()
                if vertex < size:
                    block_rows.append(vertex)
                else:
                    block_columns.append(vertex - size)
                for other in adjacency[vertex]:
                    if other not in seen:
                        seen.add(other)
                        queue.append(other)
            assert len(block_rows) == len(block_columns)
            block = [[matrix[i][j] for j in block_columns] for i in block_rows]
            profile, precision = weighted_orders(
                block, [rows[i] for i in block_rows], [columns[j] for j in block_columns],
                prime, infinity=bool(place))
            orders.extend(profile)
        result.append(tuple(sorted(orders)))
    return tuple(result)


def expand_spectators(orders, prime):
    result = [0]
    for order in orders:
        result = [value + digit * order for value in result for digit in range(prime)]
    return result


def physical_control(prime, n, ramification, charge, exponents):
    start = time.monotonic()
    tensor, clearing, endpoints = physical_filtered_tensor(prime, n, ramification, exponents)
    actual = physical_cut_profiles(tensor, prime, n, ramification, charge)
    compact = active_profiles(prime, n, ramification, charge, exponents)
    expected = (
        tuple(sorted(a + b + 2 * clearing for a in compact['zero_active']
                     for b in expand_spectators(compact['zero_spectators'], prime))),
        tuple(sorted(a + b - 2 * clearing for a in compact['infinity_active']
                     for b in expand_spectators(compact['infinity_spectators'], prime))),
    )
    assert actual == expected, (prime, n, ramification, charge, exponents, actual, expected)
    return {'prime': prime, 'n': n, 'ramification': ramification, 'charge': charge,
            'exponents': exponents, 'matrix_size': prime**(2 * n),
            'clearing_power': clearing, 'endpoint_coefficients_checked': endpoints,
            'zero_profile_multiplicities': dict(sorted(Counter(actual[0]).items())),
            'infinity_profile_multiplicities': dict(sorted(Counter(actual[1]).items())),
            'seconds': time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    factors = [exact_factor_control(prime) for prime in (2, 3, 5, 7, 11)]
    active = active_controls()
    physical = [physical_control(5, 1, 1, 0, (0, 1, 2, 2, 1)),
                physical_control(5, 1, 3, 1, (0, 3, 0, -3, 6)),
                physical_control(5, 1, 1, -1, (0, -1, 0, 1, 1))]
    if args.extended:
        physical.append(physical_control(5, 2, 1, 1, (0, 1, 2, 2, 1)))
    result = {
        'status': 'Exact coefficient and Smith-profile controls; no pressure improvement claimed.',
        'polynomial_factorizations': factors,
        'direct_active_profiles': active,
        'original_physical_tensors': physical,
        'extended': args.extended,
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'polynomial_coefficient_identities': sum(
        item['polynomial_coefficient_identities'] for item in factors),
        'active_profiles': len(active), 'physical_cut_profiles': 2 * len(physical),
        'extended': args.extended}))


if __name__ == '__main__':
    main()
