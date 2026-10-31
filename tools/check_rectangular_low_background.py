#!/usr/bin/env python3
"""Exact all-low Gaussian background profiles of the partial rectangular seed.

The ordinary proof identifies both spectra with an interval of length M-2.
This checker audits the explicit Schur complement, top inverse minor, and
small original-tensor profiles. It is not a global pressure obstruction.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from check_rectangular_f4_repair import one_particle, binary_mul, binary_power, build
from native_matrix_smith_f4 import MatrixOracle


def infinity_schur_matrix(M):
    if M < 4 or M & (M-1):
        raise ValueError('Use a power-of-two interval with M>=4.')
    G = one_particle(tuple(range(M)), M)
    degree = 2*M-2
    keep = list(range(1, M-1))+list(range(M+1, 2*M-1))
    rows = [2*a-M+1 for a in range(M)]+[-M]*M
    columns = [2*a-M+1 for a in range(M)]+[M]*M
    matrix, shifts = [], []
    for i in keep:
        line, shift_row = [], []
        for j in keep:
            # The occupied-low block has determinant w^(2M-2). All
            # remaining rows have zero in its first column.
            polynomial = (G[i][j] << degree) ^ binary_mul(
                G[i][M], G[0][j] ^ binary_mul(3, G[M][j]))
            if polynomial:
                d = polynomial.bit_length()-1
                line.append((int(bin(polynomial)[2:][::-1], 2), 0))
                shift_row.append(rows[i]+columns[j]+degree-d)
            else:
                line.append((0, 0))
                shift_row.append(0)
        matrix.append(line)
        shifts.append(shift_row)
    return matrix, shifts


def zero_cross_matrix(M, pascal=False):
    indices = range(1, M-1)
    return [[((((1 << (2*a)) if (a & b == a if pascal else a == b) else 0) ^ 5) << 1, 0)
             for b in indices] for a in indices]


def top_inverse_control(M):
    # Factoring z^(1-2M) out of the top two by two inverse block gives
    # [[u, z^(2M)+q^M], [1,u*sum(q^a*t^(M-1-a))]].
    qM = binary_power(5, M)
    assert qM ^ (1 << (2*M)) == 1
    total = 0
    for a in range(M):
        total ^= binary_power(5, a) << (2*(M-1-a))
    assert total == 1
    assert binary_mul(3, 3) ^ 1 == 4
    return {'M': M, 'inverse_top_unit_block': [[3, 1], [1, 3]],
            'inverse_top_common_shift': 1-2*M,
            'inverse_top_determinant_order': 4-4*M,
            'geometric_polynomial_sum': total}


def original_controls():
    from check_rectangular_defect_coordinates import original_cut
    M, bulk = 4, 4
    seed = build(tuple(range(M)), M, 1, 2, check_types=False)
    A, B = seed['q_x'], seed['q_u']
    records = []
    for place in ('infinity', 'zero'):
        for lam, mu in ((0, 0), (0.25, -0.5), (-0.75, 1.25)):
            matrix, rw, cw = original_cut(seed, place, (lam, mu, [0]*4, [0]*3, [0]*3))
            if place == 'infinity':
                rows = [v+A*y for y in range(bulk) for v in range(bulk)]
                columns = [u+B*x for x in range(bulk) for u in range(bulk)]
                one_particle = [value+lam+mu for value in range(5-2*M, 2*M-4, 2)]
                scalar = 4-4*M+2*(lam+mu)
            else:
                rows = [x+A*y for y in range(bulk) for x in range(bulk)]
                columns = [v+A*u for u in range(bulk) for v in range(bulk)]
                one_particle = [value+sign*(lam+mu) for sign in (-1, 1)
                                for value in range(1, 2*M-4, 2)]
                scalar = 0
            restricted = [[matrix[i][j] for j in columns] for i in rows]
            actual = MatrixOracle(restricted, len(rows)).profile(
                [rw[i] for i in rows], [cw[j] for j in columns], infinity=place == 'infinity')['orders']
            if place == 'infinity':
                actual = [value+seed['denominator_degree'] for value in actual]
            expected = [scalar]
            for value in one_particle:
                expected += [entry+value for entry in expected]
            error = max(abs(a-b) for a, b in zip(sorted(actual), sorted(expected)))
            assert len(actual) == len(expected) == 16 and error == 0
            records.append({'place': place, 'M': M, 'lambda': lam, 'mu': mu,
                            'orders_checked': len(actual), 'maximum_error': error})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records = []
    for M in (4, 8, 16, 32, 64, 128):
        matrix, shifts = infinity_schur_matrix(M)
        precision = max(512, 1 << (4*M).bit_length())
        actual = MatrixOracle(matrix, len(matrix), shifts, precision=precision).profile()['orders']
        expected = list(range(5-2*M, 2*M-4, 2))
        assert actual == expected
        zero = []
        for pascal in (False, True):
            block = zero_cross_matrix(M, pascal)
            orders = MatrixOracle(block, len(block), precision=precision).profile()['orders']
            assert orders == list(range(1, 2*M-4, 2))
            zero.append(orders)
        records.append({'M': M, 'infinity_orders': actual, 'zero_block_orders': zero,
                        'infinity_scalar_order': 4-4*M, 'zero_scalar_order': 0,
                        'top_inverse_identity': top_inverse_control(M)})
    result = {'status': 'exact all-low background equals Gaussian interval M-2 in both cuts',
              'one_particle_controls': records, 'original_tensor_controls': original_controls(),
              'scope': 'partial-repair all-low submatrix and its charge gauges; no unrestricted pressure lower bound'}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
