#!/usr/bin/env python3
"""Independent coefficient and binary block-Toeplitz controls for coupled copies."""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from native_matrix_smith_f4 import MatrixOracle
from two_copy_crossing_couplings import CoupledCopies, endpoint_power


def literal_product(first, second):
    result = 0
    for i in range(first.bit_length()):
        for j in range(second.bit_length()):
            if first >> i & 1 and second >> j & 1:
                result ^= 1 << (i+j)
    return result


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def binary_toeplitz_rank(matrix, precision):
    n = len(matrix); rows = []
    for time in range(precision):
        for i in range(n):
            first = second = 0
            for j, (low, high) in enumerate(matrix[i]):
                for degree in range(time+1):
                    a, b = (low >> degree) & 1, (high >> degree) & 1
                    offset = 2*((time-degree)*n+j)
                    # Multiplication by a+b*omega on the ordered basis(1,omega).
                    first ^= a << offset
                    first ^= b << (offset+1)
                    second ^= b << offset
                    second ^= (a ^ b) << (offset+1)
            rows.extend((first, second))
    return binary_rank(rows)


def verify():
    counts = Counter(); cases = []; rng = random.Random(31102699)
    for q, coefficient, polynomial, charged in ((2, 1, 3, False), (2, 2, 5, True), (3, 2, 3, False)):
        values = []
        coordinates = list(product(range(q), repeat=4))
        for index in range(q**4):
            x, y, u, v = ((index//q**j) % q for j in range(4))
            identity, crossing = int(x == u and y == v), int(x == v and y == u)
            values.append((identity | ((identity ^ crossing) << 1), 0))
        zero = (0,)*q
        height = tuple(range(q)) if charged else zero
        heights = (zero, zero, height, height, height, height, zero, zero)
        first = tuple(row for row in coordinates if rng.randrange(3) == 0)
        second = tuple(row for row in coordinates if rng.randrange(3) == 0)
        first_set, second_set = set(first), set(second)
        compiler = CoupledCopies(q, values, 1, heights, precision=128)
        baseline = compiler.baseline_control()
        counts['tensor_product_complete_profile_controls'] += 2
        matrices = compiler.matrices(first, second, polynomial, coefficient)
        profiles = []
        for side, matrix in enumerate(matrices):
            raw = MatrixOracle(matrix, q**4, precision=128).profile()['orders']
            assert all(order == round(order) for order in raw)
            orders = tuple(map(round, raw))
            for precision in range(1, max(orders)+2):
                actual = binary_toeplitz_rank(matrix, precision)
                assert actual == 2*sum(max(0, precision-order) for order in orders)
                counts['independent_binary_toeplitz_ranks'] += 1
            profiles.append([order-2*compiler.scalars[side] for order in orders])
        for bits1 in coordinates:
            x1, y1, u1, v1 = bits1
            p1 = values[x1+q*y1+q*q*u1+q**3*v1][0]
            for bits2 in coordinates:
                x2, y2, u2, v2 = bits2
                p2 = values[x2+q*y2+q*q*u2+q**3*v2][0]
                low, high = literal_product(p1, p2), 0
                if bits1 in first_set and bits2 in second_set:
                    low ^= polynomial if coefficient & 1 else 0
                    high ^= polynomial if coefficient & 2 else 0
                endpoint = (low.bit_count() % 2) | ((high.bit_count() % 2) << 1)
                assert endpoint == int(x1 == v1 and x2 == v2 and y1 == u1 and y2 == u2)
                counts['literal_endpoint_entries'] += 1
                for side, matrix in enumerate(matrices):
                    i1, j1 = (q*x1+u1, q*y1+v1) if side == 0 else (q*x1+y1, q*u1+v1)
                    i2, j2 = (q*x2+u2, q*y2+v2) if side == 0 else (q*x2+y2, q*u2+v2)
                    h = heights[4*side:4*side+4]
                    weight = sum(h[k][bits1[k]]+h[k][bits2[k]] for k in range(4))
                    transformed = []
                    for plane in (low, high):
                        value = 0
                        for power in range(3):
                            if plane >> power & 1:
                                exponent = 2*compiler.scalars[side]+weight+(-power if side == 0 else power)
                                value ^= 1 << exponent
                        transformed.append(value)
                    assert tuple(transformed) == matrix[q*q*i1+i2][q*q*j1+j2]
                    counts['independent_weighted_matrix_entries'] += 1
        for invalid in (1, 1 << 3):
            try:
                compiler.matrices(first, second, invalid, coefficient)
            except ValueError:
                counts['invalid_polynomial_controls'] += 1
            else:
                raise AssertionError('Invalid coupling accepted.')
        cases.append({'source_alphabet': q, 'physical_alphabet': q*q, 'F4_coefficient': coefficient,
                      'coupling_polynomial': polynomial, 'heights': heights,
                      'first_support': first, 'second_support': second,
                      'baseline_profiles': baseline, 'coupled_profiles': profiles})
    for power in range(33):
        expected = 1
        for _ in range(power):
            expected = literal_product(expected, 3)
        assert endpoint_power(power) == expected
        counts['literal_binomial_power_controls'] += 1
    root = Path(__file__).resolve().parent
    return {'scope': 'Exact finite endpoint, coefficient, weighting and complete profile controls; no asymptotic improvement asserted.',
            'checks': dict(counts), 'cases': cases,
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest()
                              for name in ('two_copy_crossing_couplings.py', Path(__file__).name)}}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
