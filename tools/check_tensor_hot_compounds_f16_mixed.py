#!/usr/bin/env python3
"""Literal polynomial and independent Toeplitz controls for mixed F16 maps."""
import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import time

import check_tensor_hot_compounds_f16 as controls
import tensor_hot_compounds_f16 as field
import tensor_hot_compounds_f16_mixed as mixed
from tensor_f16_cubic_graph_codes import validate_maps, graph_code, orthogonal_retraction


def literal_entry(modes, row, column):
    result = [0] * 4
    for rmask, a in row.items():
        rr = [i for i in range(2 * modes) if rmask >> i & 1]
        for cmask, b in column.items():
            cc = [i for i in range(2 * modes) if cmask >> i & 1]
            if len(rr) != len(cc):
                continue
            matrix = [[(3 << ((i % modes) ^ (j % modes))) if i // modes == j // modes
                       else int(i % modes == j % modes) for j in cc] for i in rr]
            exponent = sum(i % modes for i in rr + cc) + len(rr)
            polynomial = controls.binary_determinant(matrix) << exponent
            coefficient = field.mul(a, b)
            for bit in range(4):
                if coefficient >> bit & 1:
                    result[bit] ^= polynomial
    return tuple(result)


def check_case(embedding, left, modes, label):
    heights, parities = validate_maps(embedding, left, modes)
    records, entries, jets, cap_orders = [], 0, 0, 0
    for cut in ('zero', 'infinity'):
        for parity in (0, 1):
            pairs = [(i, j) for i in range(len(embedding)) for j in range(len(embedding))
                     if (parities[i] + parities[j]) % 2 == parity]
            rows = [field.tensor_vector(left[i], left[j] if cut == 'zero' else embedding[j], modes)
                    for i, j in pairs]
            columns = [field.tensor_vector(embedding[j] if cut == 'zero' else left[j], embedding[i], modes)
                       for i, j in pairs]
            native = mixed.compound_block(modes, rows, columns, coefficients=True)
            literal = [[literal_entry(modes, row, column) for column in columns] for row in rows]
            assert native['matrix'] == literal
            checks = controls.toeplitz_controls(literal, native['orders_t'])
            caps = []
            for cap in (8, 24, 64):
                capped = mixed.compound_block(modes, rows, columns, order_cap=cap)
                assert capped['orders_t'] == [min(cap, order) for order in native['orders_t']]
                caps.append({'cap': cap, 'tail_count': capped['capped_tail_count'],
                             'precision_attempts': capped['precision_attempts']})
                cap_orders += len(pairs)
            entries += len(pairs) ** 2
            jets += len(checks)
            records.append({'cut': cut, 'parity': parity, 'orders': native['orders_t'],
                            'toeplitz_checks': checks, 'capped_controls': caps})
    if label.endswith('_orthogonal'):
        zero = sorted(order for record in records if record['cut'] == 'zero' for order in record['orders'])
        infinity = sorted(order for record in records if record['cut'] == 'infinity' for order in record['orders'])
        assert zero == infinity
    return {'label': label, 'embedding': embedding, 'retraction': left, 'heights': heights,
            'literal_entries': entries, 'independent_toeplitz_checks': jets,
            'capped_order_comparisons': cap_orders, 'profiles': records}


def determinant(matrix):
    work = [row[:] for row in matrix]
    result = 1
    for j in range(len(work)):
        pivot = next((i for i in range(j, len(work)) if work[i][j]), None)
        if pivot is None:
            return 0
        work[j], work[pivot] = work[pivot], work[j]
        result = field.mul(result, work[j][j])
        inverse = field.INV[work[j][j]]
        for i in range(j + 1, len(work)):
            coefficient = field.mul(work[i][j], inverse)
            if coefficient:
                work[i] = [a ^ field.mul(coefficient, b) for a, b in zip(work[i], work[j])]
    return result


def graph_determinant_control(coefficients):
    code = graph_code(list(coefficients), 49, 8)
    cells = defaultdict(list)
    for height, vector in zip(code['heights'], code['embedding']):
        cells[height].append(vector)
    actual = 1
    for vectors in cells.values():
        gram = [[field.dot(a, b) for b in vectors] for a in vectors]
        actual = field.mul(actual, determinant(gram))
    product = 1
    contains_one = False
    for mask in range(1, 8):
        total = 0
        for i, coefficient in enumerate(coefficients):
            if mask >> i & 1:
                total ^= coefficient
        contains_one |= total == 1
        product = field.mul(product, 1 ^ total)
    expected = product
    for _ in range(3):
        expected = field.mul(expected, expected)
    assert actual == expected
    assert bool(actual) == (not contains_one)
    return {'coefficients': coefficients, 'gram_determinant': actual,
            'product_formula': expected, 'one_in_binary_span': contains_one}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-hot-compounds-f16-mixed-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    cases = []
    for coefficient in (2, 4):
        embedding = [{w: 1} for w in range(8)]
        embedding[-1][16] = coefficient
        coordinate = [{w: 1} for w in range(8)]
        scale = field.INV[field.dot(embedding[-1], embedding[-1])]
        orthogonal = embedding[:-1] + [{w: field.mul(scale, c) for w, c in embedding[-1].items()}]
        cases.append(check_case(embedding, coordinate, 5, f'coefficient_{coefficient}_coordinate'))
        cases.append(check_case(embedding, orthogonal, 5, f'coefficient_{coefficient}_orthogonal'))
    core = graph_code([2, 4, 8], 49, 8)
    core_left = orthogonal_retraction(core)
    assert core['q'] == 128 and core_left is not None
    finite_graph_controls = []
    for cap in (15, 20, 24, 25):
        code = graph_code([2, 4, 8], cap)
        left = orthogonal_retraction(code)
        assert left is not None
        finite_graph_controls.append({'height_cap': cap, 'q': code['q'],
                                      'mixed_vectors': sum(len({w.bit_count() for w in v}) > 1
                                                           for v in code['embedding'])})
    rng = random.Random(20261101)
    coefficient_controls = list(product(range(4), repeat=3))
    coefficient_controls += [tuple(rng.randrange(16) for _ in range(3)) for _ in range(64)]
    coefficient_controls += [(2, 2, 2), (2, 4, 6), (2, 4, 8), (1, 2, 4), (2, 4, 9)]
    determinant_controls = [graph_determinant_control(c) for c in sorted(set(coefficient_controls))]
    names = [Path(__file__).name, 'tensor_hot_compounds_f16_mixed.py',
             'tensor_hot_compounds_f16_mixed.cpp', 'tensor_f16_cubic_graph_codes.py',
             'tensor_hot_compounds_f16_fast.py', 'tensor_hot_compounds_f16_fast.cpp',
             'tensor_hot_compounds_f16.py', 'tensor_hot_compounds_f16.cpp',
             'check_tensor_hot_compounds_f16.py', 'tensor_orthogonal_scaling.py',
             'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'scope': __doc__, 'cases': cases,
              'literal_polynomial_entries': sum(c['literal_entries'] for c in cases),
              'independent_toeplitz_rank_checks': sum(c['independent_toeplitz_checks'] for c in cases),
              'capped_order_comparisons': sum(c['capped_order_comparisons'] for c in cases),
              'full_graph_core': core, 'full_graph_core_retraction': core_left,
              'finite_graph_controls': finite_graph_controls,
              'graph_determinant_controls': determinant_controls,
              'full_core_determinant_comparisons': len(determinant_controls),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'seconds': time.monotonic() - started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('cases', 'full_graph_core', 'full_graph_core_retraction',
                                     'graph_determinant_controls', 'source_sha256')}))


if __name__ == '__main__':
    main()
