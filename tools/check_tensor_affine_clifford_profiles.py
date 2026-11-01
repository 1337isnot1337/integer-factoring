#!/usr/bin/env python3
"""Exact controls for the auxiliary-mode Clifford profile theorem."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import tensor_affine_clifford_profiles as clifford
from check_tensor_literal_full_tail_rank_one import reference_entry
from native_matrix_smith_f4 import MatrixOracle
from tensor_hot_compounds_f16 import height, mul


def permutation_controls():
    states, operator_checks = 0, 0
    for n in range(1, 7):
        U = lambda state: state ^ ((1 << n) if (state & ((1 << n)-1)).bit_count() % 2 else 0)
        for state in range(1 << (n+1)):
            assert U(U(state)) == state
            states += 1
            for mode in range(n):
                for creation in (False, True):
                    source = U(state)
                    actual = None if bool(source >> mode & 1) == creation else U(source ^ (1 << mode))
                    expected = None if bool(state >> mode & 1) == creation else state ^ (1 << mode) ^ (1 << n)
                    assert actual == expected
                    operator_checks += 1
    return {'permutation_states': states, 'exact_CAR_state_comparisons': operator_checks}


def formal_local_formula():
    # Degree-two monomials in a,b,c,d have distinct exponents under
    # a=t, b=t^3, c=t^9, d=t^27. Thus this verifies the universal formula.
    a, b, c, d = [(1 << (3**j), 0) for j in range(4)]
    X, adj = [[a, b], [c, d]], [[d, b], [c, a]]
    U = lambda state: state ^ ((state & 1) << 1)
    large = [[X[U(i) & 1][U(j) & 1] if U(i) >> 1 == U(j) >> 1 else clifford.ZERO
              for j in range(4)] for i in range(4)]
    inverse_numerator = [[adj[U(i) & 1][U(j) & 1] if U(i) >> 1 == U(j) >> 1 else clifford.ZERO
                          for j in range(4)] for i in range(4)]
    operators = []
    for mode, creation in ((0, True), (0, False), (1, True), (1, False)):
        operator = [[clifford.ZERO]*4 for _ in range(4)]
        for state in range(4):
            if bool(state >> mode & 1) != creation:
                operator[state ^ (1 << mode)][state] = clifford.ONE
        operators.append(operator)
    expected, determinant = clifford.local_numerator(X)
    count = 0
    for column, operator in enumerate(operators):
        conjugated = clifford.matrix_product(clifford.matrix_product(large, operator), inverse_numerator)
        reconstructed = [[clifford.ZERO]*4 for _ in range(4)]
        for coefficient, basis in zip([row[column] for row in expected], operators):
            for i in range(4):
                for j in range(4):
                    reconstructed[i][j] = clifford.add(reconstructed[i][j],
                                                       clifford.multiply(coefficient, basis[i][j]))
        assert conjugated == reconstructed
        count += 16
    product_matrix = clifford.matrix_product(large, inverse_numerator)
    assert product_matrix == [[determinant if i == j else clifford.ZERO for j in range(4)] for i in range(4)]
    return {'universal_degree_two_conjugation_entries': count, 'adjugate_entries': 16}


def lifted_code(L, P, low_height):
    code = {'modes': 4, 'q': 16, 'heights': [], 'embedding': [], 'retraction': []}
    for word in range(16):
        tail, low = word & ~1, word & 1
        code['heights'].append(height(tail)+low_height*low)
        vectors = []
        for source, transpose in ((P, False), (L, True)):
            vector = {}
            for occupied in range(2):
                planes = source[low][occupied] if transpose else source[occupied][low]
                for exponent in range(max(p.bit_length() for p in planes)):
                    coefficient = ((planes[0] >> exponent) & 1) | (((planes[1] >> exponent) & 1) << 1)
                    if coefficient:
                        vector[tail | occupied, exponent] = coefficient
            vectors.append(vector)
        code['embedding'].append(vectors[0])
        code['retraction'].append(vectors[1])
    return code


def toeplitz_prefix(matrix, orders, depth=8):
    n = len(matrix)
    scalar_rows = [[sum(((mul(c, 1 << b) >> a) & 1) << b for b in range(2))
                    for a in range(2)] for c in range(4)]
    coefficients = [[0]*(2*n) for _ in range(depth)]
    for i, row in enumerate(matrix):
        for j, planes in enumerate(row):
            for exponent in range(depth):
                scalar = ((planes[0] >> exponent) & 1) | (((planes[1] >> exponent) & 1) << 1)
                for bit in range(2):
                    coefficients[exponent][2*i+bit] |= scalar_rows[scalar][bit] << (2*j)
    pivots, checks = {}, []
    for time_index in range(depth):
        for i in range(2*n):
            row = 0
            for exponent in range(time_index+1):
                row |= coefficients[exponent][i] << (2*n*(time_index-exponent))
            while row:
                pivot = row.bit_length()-1
                if pivot in pivots:
                    row ^= pivots[pivot]
                else:
                    pivots[pivot] = row
                    break
        length = time_index+1
        assert len(pivots) == 2*sum(max(0, length-order) for order in orders)
        checks.append({'jet_length': length, 'binary_rank': len(pivots)})
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('evidence/tensor-affine-clifford-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    algebra = {**permutation_controls(), **formal_local_formula()}
    evidence_path = Path('evidence/tensor-literal-full-tail-affine-m16-census.json')
    census = json.loads(evidence_path.read_text())
    assert census['completed'] and len(census['records']) == 256
    comparisons, census_records = 0, []
    for record in census['records']:
        a, b = [sum(((record[key] >> i) & 1) << (2*i) for i in range(4))
                for key in ('a_mask', 'b_mask')]
        profiles = []
        for expected in record['profiles']:
            actual = clifford.profile(record['M'], clifford.affine(a, 1), clifford.affine(b, 1),
                                       expected['cut'], record['low_heights'][1])
            assert actual['orders'] == expected['orders']
            assert 128*actual['mean_order_twice'] == sum(actual['orders'])
            comparisons += 256
            profiles.append({k: v for k, v in actual.items() if k != 'orders'})
        census_records.append({'binary_a_mask': record['a_mask'], 'binary_b_mask': record['b_mask'],
                               'profiles': profiles})
    print(json.dumps({'binary_census_order_comparisons': comparisons,
                      'seconds': time.monotonic()-started}), flush=True)
    raw_entries, raw_order_comparisons, rank_checks, raw_records = 0, 0, 0, []
    for a, b, degree, low_height in ((8, 32, 3, 2), (85, 12, 2, 1)):
        L, P = clifford.affine(a, degree), clifford.affine(b, degree)
        code, cuts = lifted_code(L, P, low_height), []
        for cut in ('zero', 'infinity'):
            actual = clifford.profile(4, L, P, cut, low_height)
            clearing = actual['clearing_power']
            matrix = [[clifford.ZERO]*256 for _ in range(256)]
            for x, y, u, v in product(range(16), repeat=4):
                value = reference_entry(code, x, y, u, v, cut, clearing)
                assert value[2:] == (0, 0)
                endpoint = (value[0].bit_count() % 2) | ((value[1].bit_count() % 2) << 1)
                assert endpoint == int(x == v and y == u)
                row, column = (16*x+y, 16*v+u) if cut == 'zero' else (16*x+u, 16*v+y)
                matrix[row][column] = value[:2]
                raw_entries += 1
            precision = 128
            while True:
                try:
                    values = MatrixOracle(matrix, 256, precision=precision).profile()['orders']
                    break
                except ArithmeticError:
                    precision *= 2
                    assert precision <= 8192
            assert all(value == round(value) for value in values)
            cleared_orders = sorted(map(round, values))
            assert actual['orders'] == [value-clearing for value in cleared_orders]
            independent_ranks = toeplitz_prefix(matrix, cleared_orders)
            raw_order_comparisons += 256
            rank_checks += len(independent_ranks)
            cuts.append({'clifford_profile': actual, 'literal_cleared_orders': cleared_orders,
                         'literal_precision': precision, 'independent_toeplitz_checks': independent_ranks})
        raw_records.append({'a_mask': a, 'b_mask': b, 'degree': degree,
                            'low_height': low_height, 'cuts': cuts})
        print(json.dumps({'completed_original_f4_cases': len(raw_records),
                          'seconds': time.monotonic()-started}), flush=True)
    names = [Path(__file__).name, 'tensor_affine_clifford_profiles.py',
             'tensor_literal_full_tail_core.py', 'check_tensor_literal_full_tail_rank_one.py',
             'tensor_literal_full_tail_rank_one.py', 'tensor_literal_rank_one_cell_search.py',
             'check_tensor_literal_hot_pfaffians_f16.py', 'check_tensor_hot_compounds_f16.py',
             'reflection_dual_rail_gaussian.py', 'tensor_hot_compounds_f16.py',
             'tensor_orthogonal_scaling.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_weighted_smith_f4.py',
             'native_rectangular_smith_f4.cpp']
    report = {'passed': True, 'algebra_controls': algebra,
              'binary_census_complete_order_comparisons': comparisons,
              'original_f4_raw_laurent_entries': raw_entries,
              'original_f4_complete_order_comparisons': raw_order_comparisons,
              'independent_f4_toeplitz_rank_checks': rank_checks,
              'binary_census_reconstruction': census_records, 'original_f4_cases': raw_records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names},
              'evidence_sha256': {str(evidence_path): sha256(evidence_path.read_bytes()).hexdigest()},
              'seconds': time.monotonic()-started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, separators=(',', ':'))+'\n')
    print(json.dumps({k: v for k, v in report.items() if k not in
                      ('binary_census_reconstruction', 'original_f4_cases', 'source_sha256', 'evidence_sha256')}),
          flush=True)


if __name__ == '__main__':
    main()
