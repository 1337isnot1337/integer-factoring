#!/usr/bin/env python3
"""Exact Gaussian-locus classification for two even-core output parameters."""
from __future__ import annotations
import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from check_even_core_postselection import parameter_basis, gaussian_relation_failures
from check_schur_downset_crossing import binomial_polynomial, sparse_kernel
from weighted_smith_oracle import polynomial_multiply, polynomial_exact_divide

def independent_parameter_basis():
    """Return A,B_y,B_v,C for independent output BCS coefficients."""
    n = 4
    kernel = sparse_kernel(tuple(range(n)), n)
    matrix = [row[:2 * n] for row in kernel[2 * n:]]

    @lru_cache(None)
    def minor(rows, columns):
        if rows.bit_count() != columns.bit_count():
            return 0
        if not rows:
            return 1
        bit = rows & -rows
        i = bit.bit_length() - 1
        answer = 0
        remaining = columns
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            j = cb.bit_length() - 1
            answer ^= polynomial_multiply(matrix[i][j], minor(rows ^ bit, columns ^ cb))
        return answer

    basis = [[0] * 256 for _ in range(4)]
    quartic_rows = (1 << 3) | (1 << 7)
    quartic_columns = 1 | (1 << 4)
    for labels in product(range(4), repeat=4):
        target = sum(a << (2 * leg) for leg, a in enumerate(labels))
        for states in product((0, 1), repeat=4):
            old = [((a & 1) << 1) | ((a >> 1) << 2) | (9 if both else 0)
                   for a, both in zip(labels, states)]
            x, y, u, v = old
            rows, columns = v | (y << n), u | (x << n)
            numerator = polynomial_multiply(5, minor(rows, columns))
            if rows & quartic_rows == quartic_rows and columns & quartic_columns == quartic_columns:
                numerator ^= minor(rows ^ quartic_rows, columns ^ quartic_columns)
            basis[states[1] + 2 * states[3]][target] ^= numerator
    # The binary polynomial 5 is c^2 = 1+w^2.
    return [[polynomial_exact_divide(value, 5) for value in row] for row in basis]


def bivariate_wick(basis, mask):
    a, b, c, d = [1 << i for i in range(8) if mask >> i & 1]
    answer = [[0] * 3 for _ in range(3)]
    for left, right in ((0, mask), (a | b, c | d), (a | c, b | d), (a | d, b | c)):
        for i, j in product(range(4), repeat=2):
            ly, lv = (i & 1) + (j & 1), (i >> 1) + (j >> 1)
            answer[ly][lv] ^= polynomial_multiply(basis[i][left], basis[j][right])
    return answer


def parameter_product(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] ^= polynomial_multiply(a, b)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def formal_axis_relations(constant, linear):
    values = [[a, b] for a, b in zip(constant, linear)]

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return [1]
        bit = mask & -mask
        rest = mask ^ bit
        answer = [0]
        remaining = rest
        while remaining:
            cb = remaining & -remaining
            remaining ^= cb
            term = parameter_product(values[bit | cb], pfaffian(rest ^ cb))
            answer += [0] * max(0, len(term) - len(answer))
            for i, value in enumerate(term):
                answer[i] ^= value
        while len(answer) > 1 and not answer[-1]:
            answer.pop()
        return answer

    checks = 0
    for mask in range(256):
        weight = mask.bit_count()
        if weight % 2 or weight < 4:
            continue
        value = values[mask]
        for _ in range(weight // 2 - 1):
            value = parameter_product(value, values[0])
        while len(value) > 1 and not value[-1]:
            value.pop()
        assert value == pfaffian(mask), mask
        checks += 1
    return checks


def exceptional_gaussian():
    c3, c5 = binomial_polynomial(3), binomial_polynomial(5)
    edges = [(0, 3, c5), (0, 6, 1 << 8), (1, 2, c3), (1, 7, 1 << 8),
             (2, 4, 1), (3, 5, 1), (4, 7, c5), (5, 6, c3)]
    values = [0] * 256
    values[0] = 1
    for i, j, polynomial in edges:
        mask = (1 << i) | (1 << j)
        old = values[:]
        for index, value in enumerate(old):
            if not index & mask:
                values[index | mask] ^= polynomial_multiply(value, polynomial)
    return [value << 8 for value in values], edges


def check(seed_directory=None):
    basis = independent_parameter_basis()
    old = parameter_basis()
    assert basis[0] == old[0] and basis[3] == old[2]
    assert [a ^ b for a, b in zip(basis[1], basis[2])] == old[1]
    degrees = [max(v.bit_length() - 1 for v in row) for row in basis]
    assert degrees == [17, 24, 32, 14]
    for i, row in enumerate(basis):
        for index, value in enumerate(row):
            x, y, u, v = [(index >> (2 * leg)) & 3 for leg in range(4)]
            assert value.bit_count() % 2 == int(i < 3 and x == v and y == u)

    r23, r83 = bivariate_wick(basis, 23), bivariate_wick(basis, 83)
    expected23, expected83 = [[0] * 3 for _ in range(3)], [[0] * 3 for _ in range(3)]
    expected23[1][2] = binomial_polynomial(4) << 10
    expected23[2][2] = binomial_polynomial(4) << 2
    expected83[2][1] = expected83[2][2] = binomial_polynomial(7) << 2
    assert r23 == expected23 and r83 == expected83

    axis_checks = formal_axis_relations(basis[0], basis[1]) + formal_axis_relations(basis[0], basis[2])
    point = [a ^ (by << 8) ^ bv ^ (c << 8) for a, by, bv, c in zip(*basis)]
    gaussian, edges = exceptional_gaussian()
    assert point == gaussian
    assert point[0] == 1 << 8 and not gaussian_relation_failures(point)

    result = {
        'status': 'exact two-parameter endpoint and complete Gaussian-locus classification; no pressure improvement claimed',
        'degrees_A_By_Bv_C': degrees,
        'coefficient_divisions': 1024,
        'endpoint_coefficient_checks': 1024,
        'equal_parameter_recombination_checks': 768,
        'axis_formal_pfaffian_checks': axis_checks,
        'exceptional_point_exterior_coefficients': 256,
        'exceptional_point_pfaffian_checks': 99,
        'exceptional_point_vacuum_polynomial_bits': point[0],
        'exceptional_point_degree': max(v.bit_length() - 1 for v in point),
        'wick_mask_23_coefficients_L_then_M': r23,
        'wick_mask_83_coefficients_L_then_M': r83,
        'basis_sha256': sha256(json.dumps(basis, separators=(',', ':')).encode()).hexdigest(),
        'exceptional_point_sha256': sha256(json.dumps(point, separators=(',', ':')).encode()).hexdigest(),
    }
    if seed_directory:
        seed_directory.mkdir(parents=True, exist_ok=True)
        (seed_directory / 'two-parameter-even-core-basis.json').write_text(json.dumps({'q': 4, 'basis': basis, 'parameter': 'A+L*B_y+M*B_v+L*M*C'}, indent=2) + '\n')
        (seed_directory / 'two-parameter-even-core-gaussian-point.json').write_text(json.dumps({'q': 4, 'values': point, 'strict_endpoint': True, 'vacuum_polynomial_bits': 1 << 8, 'edges': edges}, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('evidence/two-parameter-even-core.json'))
    parser.add_argument('--seed-directory', type=Path)
    args = parser.parse_args()
    result = check(args.seed_directory)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
