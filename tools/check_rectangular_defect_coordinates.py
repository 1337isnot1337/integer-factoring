#!/usr/bin/env python3
"""Independent M4 entry and arbitrary-core-weight controls (NumPy, C++17)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import random

import numpy as np

from check_rectangular_defect_core_basis import core_basis_data, data_for
from check_rectangular_defect_fock import (
    TABLE, mul, power, eval_binary, matrix_product, original_matrix,
    infinity_transform, zero_transform,
)
from check_rectangular_f4_repair import build
from native_matrix_smith_f4 import MatrixOracle
from rectangular_defect_coordinates import coordinate_matrix


def point_car(matrix, image, z, transpose=False):
    result = np.zeros_like(matrix)
    for creation, label in ((True, 'creation'), (False, 'annihilation')):
        for mode, (numerator, denominator) in enumerate(image[label]):
            if not numerator:
                continue
            coefficient = mul(eval_binary(numerator, z), power(eval_binary(denominator, z), -1))
            actual_creation = creation != transpose
            bit = 1 << mode
            sources = [i for i in range(len(matrix)) if bool(i & bit) != actual_creation]
            targets = [i ^ bit for i in sources]
            result[targets] ^= TABLE[matrix[sources], coefficient]
    return result


def point_basis(basis, z, side):
    n, names = basis['matrix_order'], basis['core_names']
    images = basis['regular_images']
    P = np.eye(n, dtype=np.uint8)
    for name in names:
        occupied = point_car(point_car(P, images[name+'_annihilate'], z), images[name+'_create'], z)
        P ^= occupied
    blocks = []
    for mask in range(1 << len(names)):
        value = P.copy()
        for i, name in enumerate(names):
            if mask >> i & 1:
                if side == 'left':
                    value = point_car(value.T, images[name+'_annihilate'], z, transpose=True).T
                else:
                    value = point_car(value, images[name+'_create'], z)
        blocks.append(value[basis['vacuum_rows']] if side == 'left' else value[:, basis['vacuum_columns']])
    return np.concatenate(blocks, axis=0 if side == 'left' else 1)


def evaluate_polynomial_matrix(matrix, shifts, metadata, z):
    assert metadata['exact_polynomial_matrix']
    inverse = power(eval_binary(int(metadata['unit_denominator_hex'], 16), z), -1)
    result = np.zeros((len(matrix), len(matrix)), dtype=np.uint8)
    for i, row in enumerate(matrix):
        for j, (lo, hi) in enumerate(row):
            value = eval_binary(lo, z)^mul(6, eval_binary(hi, z))
            result[i, j] = mul(mul(value, power(z, shifts[i][j])), inverse)
    return result


def point_control(data, compiled, z, seed=None):
    place = data['place']
    K = original_matrix(z, place, seed)
    K = infinity_transform(K, z, data.get('full_core', False)) if place == 'infinity' else zero_transform(K, z)
    masks = [sum(((mask >> j)&1) << index for j, index in enumerate(data['active_indices']))
             for mask in range(1 << len(data['active_indices']))]
    active = K[np.ix_(masks, masks)]
    left, right = [core_basis_data(data, side) for side in ('left', 'right')]
    transformed = matrix_product(matrix_product(point_basis(left, z, 'left'), active), point_basis(right, z, 'right'))
    infinity = place == 'infinity'
    full_alphabet = compiled[2].get('full_alphabet_gamma') is not None
    rows = [i for i in range(len(active))
            if full_alphabet or ((i//left['multiplicity'] >> 2)&3 if infinity else (i//left['multiplicity'])&3) != 3]
    columns = [j for j in range(len(active))
               if full_alphabet or ((j//right['multiplicity'])&3 if infinity else (j//right['multiplicity'] >> 2)&3) != 3]
    dropped_rows = sorted(set(range(len(active)))-set(rows))
    dropped_columns = sorted(set(range(len(active)))-set(columns))
    assert not np.any(transformed[dropped_rows])
    assert not np.any(transformed[:, dropped_columns])
    expected = transformed[np.ix_(rows, columns)]
    actual = evaluate_polynomial_matrix(*compiled, z)
    mismatch = np.argwhere(actual != expected)
    assert not len(mismatch), (place, z, mismatch[:8].tolist())
    return {'place': place, 'M': data['M'], 'z': z,
            'entries_checked': len(expected)**2,
            'zero_rows_checked': len(dropped_rows), 'zero_columns_checked': len(dropped_columns)}


def physical_masks(M, alphabet):
    bulk = 1 << (M-2)
    core = (1, 0, 1 << (M-1), (1 << (M-1))|1)
    return [((letter % bulk) << 1)|core[letter//bulk] for letter in range(alphabet)]


def coordinate_weights(metadata, gauges):
    lam, mu, v, y, u = gauges[:5]
    x = gauges[5] if len(gauges) > 5 else [0]*4
    inverse_core = (1, 0, 2, 3)
    rw, cw = [], []
    if metadata['place'] == 'infinity':
        for state, charge in zip(metadata['row_core_states'], metadata['row_charges']):
            rw.append(lam*charge+v[inverse_core[state&3]]+y[inverse_core[(state >> 2)&3]])
        for state, charge in zip(metadata['column_core_states'], metadata['column_charges']):
            cw.append(mu*charge+u[inverse_core[state&3]]+x[inverse_core[(state >> 2)&3]])
    else:
        for state, charge in zip(metadata['row_core_states'], metadata['row_charges']):
            rw.append(lam*charge+y[inverse_core[state&3]]+x[inverse_core[(state >> 2)&3]])
        for state, charge in zip(metadata['column_core_states'], metadata['column_charges']):
            cw.append(mu*charge+v[inverse_core[state&3]]+u[inverse_core[(state >> 2)&3]])
    return rw, cw


def original_cut(seed, place, gauges):
    M = seed['ambient'] if 'ambient' in seed else 4
    lam, mu, vcore, ycore, ucore = gauges[:5]
    xcore = gauges[5] if len(gauges) > 5 else [0]*4
    A, B, bulk = seed['q_x'], seed['q_u'], 1 << (M-2)
    alphabets = (A, B, B, A)
    masks = [physical_masks(M, alphabet) for alphabet in alphabets]
    if place == 'infinity':
        mode_weights = ([M+mu]*M, [-M+lam]*M,
                        [2*a-M+1+mu for a in range(M)], [2*a-M+1+lam for a in range(M)])
    else:
        mode_weights = ([-2*M+1+2*a+lam for a in range(M)], [1+2*a-lam for a in range(M)],
                        [-mu]*M, [mu]*M)
    core_weights = (xcore, ycore, ucore, vcore)
    weights = [[sum(mode_weights[leg][a] for a in range(M) if mask >> a & 1)+core_weights[leg][letter//bulk]
                for letter, mask in enumerate(masks[leg])] for leg in range(4)]
    n = A*B
    matrix = [[(0, 0)]*n for _ in range(n)]
    for x in range(A):
        for y in range(B):
            for u in range(B):
                for v in range(A):
                    row, column = ((v+A*y, u+B*x) if place == 'infinity' else (x+A*y, v+A*u))
                    matrix[row][column] = seed['values'][x+A*(y+B*(u+B*v))]
    if place == 'infinity':
        rw = [weights[3][v]+weights[1][y] for y in range(B) for v in range(A)]
        cw = [weights[2][u]+weights[0][x] for x in range(A) for u in range(B)]
    else:
        rw = [weights[0][x]+weights[1][y] for y in range(B) for x in range(A)]
        cw = [weights[3][v]+weights[2][u] for u in range(B) for v in range(A)]
    return matrix, rw, cw


def weighted_controls(data, compiled, count, seed=None):
    if seed is None:
        seed = build(tuple(range(4)), 4, 1, 2, check_types=False)
    matrix, shifts, metadata = compiled
    reduced = MatrixOracle(matrix, len(matrix), shifts)
    original, _, _ = original_cut(seed, data['place'], (0, 0, [0]*4, [0]*3, [0]*3))
    dense = MatrixOracle(original, len(original))
    rng = random.Random(18721 if data['place'] == 'infinity' else 9713)
    records = []
    for case in range(count):
        full = data.get('full_core', False)
        number = 16 if full else 12
        values = [rng.randrange(-8, 9)/4 for _ in range(number)] if case else [0]*number
        gauges = (values[0], values[1], values[2:6], values[6:9], values[9:12])
        if full:
            gauges += (values[12:16],)
        rw, cw = coordinate_weights(metadata, gauges)
        expected = reduced.profile(rw, cw)['orders']
        inactive = [i for i in range(2*data['M']) if i not in data['active_indices']]
        for index, order in zip(inactive, data['spectator_orders']):
            sign = 1 if data['place'] == 'infinity' or index < data['M'] else -1
            spectator = order+sign*(gauges[0]+gauges[1])
            expected = sorted(expected+[value+spectator for value in expected])
        _, rw, cw = original_cut(seed, data['place'], gauges)
        actual = dense.profile(rw, cw, infinity=data['place'] == 'infinity')['orders']
        if data['place'] == 'infinity':
            actual = [value+seed['denominator_degree'] for value in actual]
        error = max(abs(a-b) for a, b in zip(sorted(actual), sorted(expected)))
        assert len(actual) == len(expected) and error < 1e-9, (data['place'], case, error)
        records.append({'case': case, 'orders_checked': len(actual), 'maximum_error': error,
                        'lambda': gauges[0], 'mu': gauges[1], 'v': gauges[2], 'y': gauges[3], 'u': gauges[4]})
        if full:
            records[-1]['x'] = gauges[5]
    return {'place': data['place'], 'M': 4, 'cases': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--weighted-cases', type=int, default=5)
    args = parser.parse_args()
    points, weighted = [], []
    for place in ('infinity', 'zero'):
        data = data_for(place, 4)
        compiled = coordinate_matrix(data)
        points.extend(point_control(data, compiled, z) for z in (2, 3, 4, 5))
        weighted.append(weighted_controls(data, compiled, args.weighted_cases))
    result = {'status': 'independent original-tensor coordinate and real-core-gauge controls',
              'field': {'modulus': 'X^4+X+1', 'omega_image': 6},
              'entrywise_cases': points, 'weighted_profiles': weighted}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
