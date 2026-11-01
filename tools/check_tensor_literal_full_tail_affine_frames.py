#!/usr/bin/env python3
"""Check endpoint-only affine full-tail frames against the original tensor.

The reference contracts the original raw Laurent Pfaffian coefficients.
Binary Toeplitz ranks use ordinary Python bit elimination and do not use
the native Smith implementation or the concentrated operator formulas.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

from check_tensor_literal_full_tail_rank_one import reference_entry
from native_matrix_smith_f4 import MatrixOracle
from tensor_literal_full_tail_affine_frames import affine, profile
from tensor_hot_compounds_f16 import height


def lifted_code(a_mask, b_mask, low_height):
    left, embedding = affine(a_mask), affine(b_mask)
    result = {'modes': 4, 'q': 16, 'heights': [], 'embedding': [], 'retraction': []}
    for virtual in range(16):
        tail, low = virtual & ~1, virtual & 1
        result['heights'].append(height(tail)+low_height*low)
        p, ell = {}, {}
        for occupied in range(2):
            for exponent in range(2):
                if embedding[occupied][low] >> exponent & 1:
                    p[tail | occupied, exponent] = 1
                if left[low][occupied] >> exponent & 1:
                    ell[tail | occupied, exponent] = 1
        result['embedding'].append(p)
        result['retraction'].append(ell)
    return result


def binary_toeplitz_controls(matrix, orders):
    n, maximum = len(matrix), max(orders)+1
    coefficients = [[0]*n for _ in range(maximum)]
    for i, row in enumerate(matrix):
        for j, polynomial in enumerate(row):
            while polynomial:
                bit = polynomial & -polynomial
                exponent = bit.bit_length()-1
                if exponent < maximum:
                    coefficients[exponent][i] |= 1 << j
                polynomial ^= bit
    pivots, checks = {}, []
    for time_index in range(maximum):
        for i in range(n):
            row = 0
            for exponent in range(time_index+1):
                row |= coefficients[exponent][i] << (n*(time_index-exponent))
            while row:
                pivot = row.bit_length()-1
                if pivot in pivots:
                    row ^= pivots[pivot]
                else:
                    pivots[pivot] = row
                    break
        length = time_index+1
        expected = sum(max(0, length-order) for order in orders)
        assert len(pivots) == expected, (length, len(pivots), expected)
        checks.append({'jet_length': length, 'binary_rank': len(pivots)})
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-literal-full-tail-affine-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    entries, order_checks, rank_checks, endpoint_checks, records = 0, 0, 0, 0, []
    cases = [(0, 0, 1), (1, 0, 1), (2, 4, 1), (15, 1, 2)]
    for a_mask, b_mask, low_height in cases:
        code = lifted_code(a_mask, b_mask, low_height)
        cut_records = []
        for cut in ('zero', 'infinity'):
            actual = profile(4, a_mask, b_mask, cut, low_height)
            clearing = actual['clearing_power']
            matrix = [[0]*256 for _ in range(256)]
            for x, y, u, v in product(range(16), repeat=4):
                value = reference_entry(code, x, y, u, v, cut, clearing)
                assert value[1:] == (0, 0, 0)
                polynomial = value[0]
                assert polynomial.bit_count() % 2 == int(x == v and y == u)
                row, column = (16*x+y, 16*v+u) if cut == 'zero' else (16*x+u, 16*v+y)
                matrix[row][column] = polynomial
                entries += 1
                endpoint_checks += 1
            converted = [[(value, 0) for value in row] for row in matrix]
            precision = 128
            while True:
                try:
                    expected = MatrixOracle(converted, 256, precision=precision).profile()['orders']
                    break
                except ArithmeticError:
                    precision *= 2
                    assert precision <= 8192
            assert all(value == round(value) for value in expected)
            expected = sorted(map(round, expected))
            assert expected == actual['cleared_orders']
            assert [value-clearing for value in expected] == actual['orders']
            independent_ranks = binary_toeplitz_controls(matrix, expected)
            order_checks += 256
            rank_checks += len(independent_ranks)
            cut_records.append({'active_profile': actual, 'literal_cleared_orders': expected,
                                'literal_precision': precision, 'binary_toeplitz_checks': independent_ranks})
        records.append({'a_mask': a_mask, 'b_mask': b_mask, 'low_heights': [0, low_height],
                        'cuts': cut_records})
        print(json.dumps({'completed_cases': len(records), 'declared_cases': len(cases),
                          'a_mask': a_mask, 'b_mask': b_mask, 'low_height': low_height,
                          'elapsed_seconds': time.monotonic()-started}), flush=True)
    names = [Path(__file__).name, 'tensor_literal_full_tail_affine_frames.py',
             'tensor_literal_full_tail_core.py', 'check_tensor_literal_full_tail_rank_one.py',
             'tensor_literal_full_tail_rank_one.py', 'check_tensor_literal_hot_pfaffians_f16.py',
             'reflection_dual_rail_gaussian.py', 'tensor_hot_compounds_f16.py',
             'tensor_orthogonal_scaling.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_weighted_smith_f4.py',
             'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'cases': records, 'raw_laurent_literal_entries': entries,
              'complete_smith_order_comparisons': order_checks,
              'independent_binary_toeplitz_rank_checks': rank_checks,
              'endpoint_checks': endpoint_checks,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'seconds': time.monotonic()-started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ('cases', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
