#!/usr/bin/env python3
"""Exact four-parameter even-core endpoint and Gaussian-surface identities."""
from __future__ import annotations
import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from check_schur_downset_crossing import binomial_polynomial, sparse_kernel
from weighted_smith_oracle import polynomial_multiply, polynomial_exact_divide
from check_two_parameter_even_core import independent_parameter_basis, exceptional_gaussian

def four_parameter_basis():
    """Return all 16 formal coefficients indexed by occupation bits x,y,u,v."""
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

    basis = [[0] * 256 for _ in range(16)]
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
            basis[sum(a << leg for leg, a in enumerate(states))][target] ^= numerator
    # The binary polynomial 5 is c^2 = 1+w^2.
    return [[polynomial_exact_divide(value, 5) for value in row] for row in basis]


def wick(basis, mask):
    a, b, c, d = [1 << i for i in range(8) if mask >> i & 1]
    result = {}
    for left, right in ((0, mask), (a | b, c | d), (a | c, b | d), (a | d, b | c)):
        for i, j in product(range(16), repeat=2):
            powers = tuple((i >> leg & 1) + (j >> leg & 1) for leg in range(4))
            result[powers] = result.get(powers, 0) ^ polynomial_multiply(basis[i][left], basis[j][right])
    return {powers: value for powers, value in result.items() if value}


def check(seed_directory=None):
    basis = four_parameter_basis()
    for coefficient, row in enumerate(basis):
        for index, value in enumerate(row):
            x, y, u, v = [(index >> (2 * leg)) & 3 for leg in range(4)]
            assert value.bit_count() % 2 == int(coefficient in (0, 6, 9) and x == v and y == u)

    two = [[0] * 256 for _ in range(4)]
    for coefficient, row in enumerate(basis):
        target = (coefficient >> 1 & 1) + 2 * (coefficient >> 3 & 1)
        for i, value in enumerate(row):
            two[target][i] ^= value
    assert two == independent_parameter_basis()

    r23, r83 = wick(basis, 23), wick(basis, 83)
    assert r23 == {(1, 1, 1, 2): binomial_polynomial(4) << 10,
                   (1, 2, 2, 2): binomial_polynomial(4) << 2}
    assert r83 == {(1, 2, 1, 1): binomial_polynomial(7) << 2,
                   (1, 2, 2, 2): binomial_polynomial(7) << 2}

    # Clear r^2 after L=w^8/r and M=1/r. Its sole surviving
    # coefficient in p,r is p*r times the known Gaussian tensor.
    surface = [[[0] * 256 for _ in range(4)] for _ in range(2)]
    for coefficient, row in enumerate(basis):
        exponent = 2 + (coefficient >> 2 & 1) - (coefficient >> 1 & 1) - (coefficient >> 3 & 1)
        for i, value in enumerate(row):
            surface[coefficient & 1][exponent][i] ^= value << (8 * (coefficient >> 1 & 1))
    gaussian, _ = exceptional_gaussian()
    for p, rows in enumerate(surface):
        for r, row in enumerate(rows):
            assert row == (gaussian if (p, r) == (1, 1) else [0] * 256)

    result = {
        'status': 'exact four-parameter endpoint and complete Gaussian-locus classification; no pressure improvement claimed',
        'parameter_order': ['p_x', 'L_y', 'r_u', 'M_v'],
        'coefficient_divisions': 4096,
        'endpoint_coefficient_checks': 4096,
        'nonzero_endpoint_parameter_masks': [0, 6, 9],
        'two_parameter_recombination_coefficients': 1024,
        'cleared_gaussian_surface_coefficients': 2048,
        'wick_mask_23': [{'powers': powers, 'polynomial_bits': value} for powers, value in sorted(r23.items())],
        'wick_mask_83': [{'powers': powers, 'polynomial_bits': value} for powers, value in sorted(r83.items())],
        'basis_sha256': sha256(json.dumps(basis, separators=(',', ':')).encode()).hexdigest(),
    }
    if seed_directory:
        seed_directory.mkdir(parents=True, exist_ok=True)
        (seed_directory / 'four-parameter-even-core-basis.json').write_text(json.dumps({'q': 4, 'basis': basis, 'parameter_bits': 'x,y,u,v'}, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('evidence/four-parameter-even-core.json'))
    parser.add_argument('--seed-directory', type=Path)
    args = parser.parse_args()
    result = check(args.seed_directory)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
