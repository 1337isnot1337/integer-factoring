#!/usr/bin/env python3
"""Independent exact controls for the F16 coherent-compound implementation.

Literal minors use binary integer polynomials and permutation determinants.
Complete Smith orders are checked through incremental binary Toeplitz ranks,
using restriction of scalars from F16 to F2. No native Smith code is used
for those rank calculations.
"""
import argparse
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import time

import coherent_hot_compounds_f4 as f4
import tensor_hot_compounds_f16 as f16


def binary_product(a, b):
    result = 0
    while a:
        bit = a & -a
        result ^= b << (bit.bit_length() - 1)
        a ^= bit
    return result


def binary_determinant(matrix):
    result = 0
    for order in permutations(range(len(matrix))):
        value = 1
        for i, j in enumerate(order):
            value = binary_product(value, matrix[i][j])
        result ^= value
    return result


def literal_entry(modes, row, column):
    planes = [0] * 4
    for rmask, a in row.items():
        rr = [i for i in range(2 * modes) if rmask >> i & 1]
        for cmask, b in column.items():
            cc = [i for i in range(2 * modes) if cmask >> i & 1]
            assert len(rr) == len(cc)
            matrix = [[(3 << ((i % modes) ^ (j % modes))) if i // modes == j // modes
                       else int(i % modes == j % modes) for j in cc] for i in rr]
            exponent = sum(i % modes for i in rr + cc) + len(rr)
            polynomial = binary_determinant(matrix) << exponent
            coefficient = f16.mul(a, b)
            for plane in range(4):
                if coefficient >> plane & 1:
                    planes[plane] ^= polynomial
    return tuple(planes)


def toeplitz_controls(matrix, orders):
    n = len(matrix)
    maximum = max(orders) + 1
    # Multiplication by one scalar as a four-by-four binary matrix.
    scalar_rows = [[sum(((f16.mul(c, 1 << b) >> a) & 1) << b for b in range(4))
                    for a in range(4)] for c in range(16)]
    coefficient_rows = [[0] * (4 * n) for _ in range(maximum)]
    for i, row in enumerate(matrix):
        for j, polynomial in enumerate(row):
            for degree in range(maximum):
                coefficient = sum(((plane >> degree) & 1) << bit
                                  for bit, plane in enumerate(polynomial))
                if coefficient:
                    for bit in range(4):
                        coefficient_rows[degree][4 * i + bit] |= scalar_rows[coefficient][bit] << (4 * j)
    pivots, checks = {}, []
    for time_index in range(maximum):
        for i in range(4 * n):
            row = 0
            for degree in range(time_index + 1):
                row |= coefficient_rows[degree][i] << (4 * n * (time_index - degree))
            while row:
                pivot = row.bit_length() - 1
                if pivot in pivots:
                    row ^= pivots[pivot]
                else:
                    pivots[pivot] = row
                    break
        length = time_index + 1
        expected = 4 * sum(max(0, length - order) for order in orders)
        assert len(pivots) == expected, (n, length, len(pivots), expected)
        checks.append({'jet_length': length, 'binary_rank': len(pivots)})
    return checks


def one_case(embedding, retraction, modes, label):
    heights, charges = f16.validate_maps(embedding, retraction, modes)
    blocks, all_orders, entries, jet_checks = [], {}, 0, 0
    for cut in ('zero', 'infinity'):
        orders = []
        for charge in sorted({a + b for a in charges for b in charges}):
            pairs = [(i, j) for i in range(len(embedding)) for j in range(len(embedding))
                     if charges[i] + charges[j] == charge]
            rows = [f16.tensor_vector(retraction[i], retraction[j] if cut == 'zero' else embedding[j], modes)
                    for i, j in pairs]
            columns = [f16.tensor_vector(embedding[j] if cut == 'zero' else retraction[j], embedding[i], modes)
                       for i, j in pairs]
            native = f16.compound_block(modes, rows, columns, coefficients=True)
            literal = [[literal_entry(modes, row, column) for column in columns] for row in rows]
            assert literal == native['matrix']
            checks = toeplitz_controls(literal, native['orders_t'])
            cap_controls = []
            for cap in (8, 16, 24):
                capped = f16.compound_block(modes, rows, columns, order_cap=cap)
                assert capped['orders_t'] == [min(order, cap) for order in native['orders_t']]
                assert capped['exact_orders_below_cap'] + capped['capped_tail_count'] == len(rows)
                cap_controls.append({'cap': cap, 'tail_count': capped['capped_tail_count']})
            entries += len(rows) ** 2
            jet_checks += len(checks)
            orders.extend(native['orders_t'])
            blocks.append({'cut': cut, 'charge': charge, 'size': len(rows),
                           'native_orders': native['orders_t'], 'toeplitz_checks': checks,
                           'capped_profile_controls': cap_controls})
        all_orders[cut] = sorted(orders)
    assert len(all_orders['zero']) == len(embedding) ** 2
    return {'label': label, 'modes': modes, 'embedding': embedding, 'retraction': retraction,
            'physical_heights': heights, 'q': len(embedding), 'literal_entries': entries,
            'toeplitz_rank_checks': jet_checks, 'complete_orders': all_orders, 'blocks': blocks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-hot-compounds-f16-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    for a in range(16):
        for b in range(16):
            assert f16.library(2).tensor_f16_scalar_control(a, b) == f16.mul(a, b)
    controls = []
    f4_order_comparisons = 0
    for coefficient in (2, 3, 4, 5):
        p = [{0: 1}] + [{1 << i: 1} for i in range(4)] + [{9: 1, 6: coefficient}]
        inv = f16.INV[f16.dot(p[-1], p[-1])]
        left = p[:-1] + [{word: f16.mul(inv, value) for word, value in p[-1].items()}]
        result = one_case(p, left, 4, f'six_letters_coefficient_{coefficient}')
        assert result['complete_orders']['zero'] == result['complete_orders']['infinity']
        if coefficient < 4:
            reference = f4.profiles(p, left, 4)
            assert result['complete_orders']['zero'] == reference['zero_orders_t']
            assert result['complete_orders']['infinity'] == reference['infinity_orders_t']
            f4_order_comparisons += 2 * len(p) ** 2
        controls.append(result)
    code = f16.differential_kernel([{9: 1, 6: 4}], 4, 8)
    assert code['q'] == 8
    left = f16.orthogonal_retraction(code)
    assert left is not None
    controls.append(one_case(code['embedding'], left, 4, 'eight_letter_F16_differential_kernel'))
    hyperoval = [[1, 1, 1, 2, 0, 0], [0, 1, 2, 1, 4, 0], [0, 1, 3, 3, 0, 5]]
    hyperoval_certificate = f16.subset_gram_certificate(hyperoval)
    assert hyperoval_certificate['all_subsets_nondegenerate']
    assert all(record['rank'] == min(3, record['active_mask'].bit_count())
               for record in hyperoval_certificate['subsets'])
    names = [Path(__file__).name, 'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'tensor_orthogonal_scaling.py', 'coherent_hot_compounds_f4.cpp',
             'coherent_hot_compounds_f4.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'scope': __doc__, 'scalar_multiplication_checks': 256,
              'literal_polynomial_entries': sum(c['literal_entries'] for c in controls),
              'independent_toeplitz_rank_checks': sum(c['toeplitz_rank_checks'] for c in controls),
              'complete_smith_orders_checked': sum(2 * c['q'] ** 2 for c in controls),
              'capped_order_comparisons': sum(6 * c['q'] ** 2 for c in controls),
              'embedded_f4_profile_order_comparisons': f4_order_comparisons,
              'controls': controls, 'f16_hyperoval_coefficients': hyperoval,
              'hyperoval_subset_gram_certificate': hyperoval_certificate,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'seconds': time.monotonic() - started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('controls', 'hyperoval_subset_gram_certificate', 'source_sha256')}))


if __name__ == '__main__':
    main()
