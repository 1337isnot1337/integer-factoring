#!/usr/bin/env python3
"""Exact finite controls for characteristic-p truncated symmetric crossings."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from math import factorial
from pathlib import Path
import random

from check_all_field_pascal_crossings import (
    add, multiply, scale, trim, pascal_kernel, normalized_infinity,
    ramify_negative, smith_orders, weighted_orders,
)


def basis(size, prime):
    return list(product(range(prime), repeat=size))


def induced(matrix, prime):
    """Monomial matrix of linear substitution modulo x_i**p."""
    r = len(matrix)
    words = basis(r, prime)
    index = {word: i for i, word in enumerate(words)}
    result = [[()] * len(words) for _ in words]
    for column, powers in enumerate(words):
        polynomial = {(0,)*r: (1,)}
        for j, power in enumerate(powers):
            for _ in range(power):
                updated = {}
                for exponents, coefficient in polynomial.items():
                    for i in range(r):
                        if exponents[i] == prime-1 or not matrix[i][j]:
                            continue
                        target = list(exponents)
                        target[i] += 1
                        target = tuple(target)
                        value = multiply(coefficient, matrix[i][j], prime)
                        updated[target] = add(updated.get(target, ()), value, prime)
                polynomial = {key: value for key, value in updated.items() if value}
        for exponents, coefficient in polynomial.items():
            result[index[exponents]][column] = coefficient
    return result


def matmul(left, right, prime):
    result = [[()] * len(right[0]) for _ in left]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if not value:
                continue
            for j, other in enumerate(right[k]):
                if other:
                    result[i][j] = add(result[i][j], multiply(value, other, prime), prime)
    return result


def direct_kernel(matrix, prime):
    r = len(matrix)
    polynomial = {(0,)*(2*r): (1,)}
    for i in range(r):
        for j in range(r):
            powers = [(1,)]
            for _ in range(1, prime):
                powers.append(multiply(powers[-1], matrix[i][j], prime))
            updated = {}
            for exponents, coefficient in polynomial.items():
                for k in range(min(prime-exponents[i], prime-exponents[r+j])):
                    target = list(exponents)
                    target[i] += k
                    target[r+j] += k
                    target = tuple(target)
                    value = scale(multiply(coefficient, powers[k], prime),
                                  pow(factorial(k), -1, prime), prime)
                    updated[target] = add(updated.get(target, ()), value, prime)
            polynomial = {key: value for key, value in updated.items() if value}
    return polynomial


def polynomial_value(value, point, prime):
    result = 0
    for coefficient in reversed(value):
        result = (result*point+coefficient) % prime
    return result


def tensor_from_induced(matrix, n, prime):
    d = prime**n
    return {(v, y, u, x): matrix[v*d+y][u*d+x]
            for v, y, u, x in product(range(d), repeat=4)
            if matrix[v*d+y][u*d+x]}


def endpoint_controls(indices, copies, prime):
    n = len(indices)
    ambient = 2*max(indices)+2
    kernel = pascal_kernel(indices, ambient, prime)
    values = [[trim([polynomial_value(value, -1, prime)], prime)
               for value in row] for row in kernel]
    local = induced(values, prime)
    d = prime**n
    digits = basis(n, prime)
    words = list(product(range(d), repeat=copies))
    totals = {word: tuple(sum(digits[letter][i] for letter in word)
                         for i in range(n)) for word in words}
    pairs = [(x, v) for x in words for v in words if totals[x] == totals[v]]
    count = 0
    for x, v in pairs:
        for y, u in pairs:
            actual = 1
            for i in range(copies):
                entry = local[v[i]*d+y[i]][u[i]*d+x[i]]
                actual = actual*(entry[0] if entry else 0) % prime
            assert actual == int(x == v and y == u)
            count += 1
    return {'prime': prime, 'indices': indices, 'copies': copies,
            'independent_type_endpoint_entries': count}


def cut_profile_controls(indices, prime):
    n = len(indices)
    ambient = 2*max(indices)+2
    d = prime**n
    digits = basis(n, prime)
    matrix = induced(pascal_kernel(indices, ambient, prime), prime)
    tensor = {labels: ramify_negative(value, prime)
              for labels, value in tensor_from_induced(matrix, n, prime).items()}
    heights = ([0]*n, [2*a+1 for a in indices], [0]*n,
               [2*a-2*ambient+1 for a in indices],
               [2*a-n+2 for a in indices], [-n]*n,
               [2*a-n+1 for a in indices], [2*ambient-n-1]*n)
    one_infinity, _ = smith_orders(normalized_infinity(indices, prime), prime)
    expected_one = [[2*a+1 for a in indices for _ in range(2)],
                    [2*order+1-2*n for order in one_infinity]]
    cases = []
    for endpoint, order in enumerate(((3, 1, 0, 2), (3, 2, 0, 1))):
        full = [[()] * (d*d) for _ in range(d*d)]
        for labels, value in tensor.items():
            full[labels[order[0]]*d+labels[order[1]]][labels[order[2]]*d+labels[order[3]]] = value
        rows, columns, charges = [], [], []
        for a, b in product(range(d), repeat=2):
            charges.append(sum(digits[a])+(1 if endpoint else -1)*sum(digits[b]))
            rows.append(sum(heights[4*endpoint+order[k]][i]*digits[letter][i]
                            for k, letter in enumerate((a, b)) for i in range(n)))
            columns.append(sum(heights[4*endpoint+order[k+2]][i]*digits[letter][i]
                               for k, letter in enumerate((a, b)) for i in range(n)))
        assert all(charges[i] == charges[j] or not full[i][j]
                   for i in range(d*d) for j in range(d*d))
        actual, precisions = [], []
        for charge in sorted(set(charges)):
            ids = [i for i, value in enumerate(charges) if value == charge]
            block = [[full[i][j] for j in ids] for i in ids]
            orders, precision = weighted_orders(block, [rows[i] for i in ids],
                                                [columns[i] for i in ids], prime,
                                                infinity=bool(endpoint))
            actual.extend(orders)
            precisions.append(precision)
        expected = sorted(sum(a*b for a, b in zip(word, expected_one[endpoint]))
                          for word in basis(2*n, prime))
        assert sorted(actual) == expected
        cases.append({'place': 'infinity' if endpoint else 'zero',
                      'one_particle_orders': expected_one[endpoint],
                      'induced_orders': expected, 'maximum_precision': max(precisions)})
    return {'prime': prime, 'indices': indices, 'cut_profiles': cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = random.Random(880173)
    functor = []
    for prime, size in ((2, 2), (3, 2), (5, 2), (3, 3)):
        count = 0
        for _ in range(6):
            left, right = [[[trim([rng.randrange(prime)], prime)
                             for _ in range(size)] for _ in range(size)] for _ in range(2)]
            assert induced(matmul(left, right, prime), prime) == \
                   matmul(induced(left, prime), induced(right, prime), prime)
            count += (prime**size)**2
        functor.append({'prime': prime, 'one_particle_size': size,
                        'composition_entries_checked': count})
    kernels = []
    for prime in (2, 3, 5):
        matrix = [[(1, 1), (0, 1)], [(1,), (1, 0, 1)]]
        actual, kernel = induced(matrix, prime), direct_kernel(matrix, prime)
        words = basis(2, prime)
        for i, alpha in enumerate(words):
            for j, beta in enumerate(words):
                multiplier = factorial(beta[0])*factorial(beta[1])
                assert actual[i][j] == scale(kernel.get(alpha+beta, ()), multiplier, prime)
        orders, precision = smith_orders(actual, prime)
        one, _ = smith_orders(matrix, prime)
        expected = sorted(sum(a*b for a, b in zip(word, one)) for word in words)
        assert list(orders) == expected
        kernels.append({'prime': prime, 'kernel_entries_checked': len(words)**2,
                        'one_particle_orders': one, 'induced_orders': orders,
                        'precision': precision})
    endpoints = [endpoint_controls(indices, copies, prime)
                 for prime, indices, copies in ((2, (0,), 3), (3, (0,), 3),
                                                (5, (0,), 2), (3, (0, 1), 2),
                                                (3, (0, 2), 2))]
    profiles = [cut_profile_controls(indices, prime)
                for prime, indices in ((2, (0,)), (3, (0,)), (5, (0,)),
                                       (3, (0, 1)), (3, (0, 2)))]
    result = {'status': 'Exact finite controls for truncated symmetric crossings; no pressure improvement claimed.',
              'functoriality': functor, 'kernel_and_smith_controls': kernels,
              'typed_endpoints': endpoints, 'actual_cut_profiles': profiles,
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(json.dumps({'functor_entries': sum(x['composition_entries_checked'] for x in functor),
                      'endpoint_entries': sum(x['independent_type_endpoint_entries'] for x in endpoints),
                      'actual_cut_cases': 2*len(profiles)}))


if __name__ == '__main__':
    main()
