#!/usr/bin/env python3
"""Independent selected F4 minors, complete profiles and symbol derivatives."""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from f4_pascal_flags import change_flag
from f4_pascal_coordinate_sectors import (F4SectorPascalLetterPressure,
                                         coordinate_exterior_tensor)
from pascal_coordinate_sectors import SectorPascalLetterPressure
from permuted_pascal_kernels import permuted_pascal_kernel
from weighted_smith_oracle_f4 import WeightedSmithOracle, polynomial_multiply
from windowed_f4_smith import pack


def independent_entries(seed, alphabet):
    n = seed['n']
    matrix = [[pack(v) if isinstance(v, (tuple, list)) else pack((v, 0)) for v in row]
              for row in seed['matrix']]

    @lru_cache(None)
    def minor(row_mask, column_mask):
        rows = [i for i in range(2*n) if row_mask >> i & 1]
        columns = [i for i in range(2*n) if column_mask >> i & 1]
        if len(rows) != len(columns):
            return 0
        states = {0: 1}
        for row in rows:
            following = {}
            for used, value in states.items():
                for i, column in enumerate(columns):
                    if not used >> i & 1 and matrix[row][column]:
                        key = used | (1 << i)
                        following[key] = following.get(key, 0) ^ polynomial_multiply(value, matrix[row][column])
            states = following
        return states.get((1 << len(rows))-1, 0)
    return {labels: minor(labels[2] | (labels[3] << n), labels[0] | (labels[1] << n))
            for labels in product(alphabet, repeat=4)}


def endpoint_expected(seed, labels):
    b, x, v, y = labels
    if b != y or x.bit_count() != v.bit_count():
        return 0
    columns = [j for j in range(seed['n']) if x >> j & 1]
    basis = {}
    for i in range(seed['n']):
        if not v >> i & 1:
            continue
        row = sum(seed['pascal'][i][j] << a for a, j in enumerate(columns))
        while row:
            pivot = row.bit_length()-1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return int(len(basis) == len(columns))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng, records = np.random.default_rng(0xF4C005), []
    entries_checked, orders_checked, binary_orders = 0, 0, 0
    for trial, (images, alphabet) in enumerate((([0,1], [0,1,2,3]),
                                               ([0,2,1,3], [0,1,3,7,15]),
                                               ([0,1,2,3], [0,1,2,3,5,7,11,15]))):
        original = permuted_pascal_kernel(images)
        flag = [[int(i == j) for j in range(len(images))] for i in range(len(images))]
        flag[0][1] = 2
        if len(images) > 2:
            flag[2][1], flag[2][3] = 3, 1
        seed = change_flag(original, flag)
        tensor = coordinate_exterior_tensor(seed, alphabet)
        independent = independent_entries(seed, alphabet)
        for labels, value in independent.items():
            pair = tensor['values'].get(labels, (0,0))
            assert pack(pair) == value
            endpoint = pair[0].bit_count() % 2 | ((pair[1].bit_count() % 2) << 1)
            assert endpoint == endpoint_expected(seed, labels)
        entries_checked += len(independent)
        model = F4SectorPascalLetterPressure(seed, alphabet)
        parameters = model.natural(1.)
        parameters[:-1] += rng.integers(-2, 3, len(parameters)-1)
        model.evaluate(parameters)
        literal_profiles = [[int(v) for v in cut['orders']]
                            for cut in model.last['complete_original_profiles']]
        for place, (_, legs) in enumerate(model.cuts):
            q, index = model.q, {v: i for i, v in enumerate(alphabet)}
            matrix = [[0]*(q*q) for _ in range(q*q)]
            for labels, value in independent.items():
                i = index[labels[legs[0]]]+q*index[labels[legs[1]]]
                j = index[labels[legs[2]]]+q*index[labels[legs[3]]]
                matrix[i][j] = value
            heights = parameters[:-1].reshape(2,4,q).astype(int)
            row = [int(heights[place,legs[0],a]+heights[place,legs[1],b]) for b in range(q) for a in range(q)]
            col = [int(heights[place,legs[2],a]+heights[place,legs[3],b]) for b in range(q) for a in range(q)]
            exact, _ = WeightedSmithOracle(matrix, infinity=bool(place)).profile(row, col)
            assert list(exact) == model.last['complete_original_profiles'][place]['orders']
            orders_checked += len(exact)
        binary = F4SectorPascalLetterPressure(original, alphabet)
        old = SectorPascalLetterPressure(original, alphabet)
        binary.evaluate(parameters)
        old.evaluate(parameters)
        assert binary.last['complete_original_profiles'] == old.last['complete_original_profiles']
        binary_orders += 2*len(alphabet)**2
        control = model.fitted_natural_start()
        control[:-1] += rng.normal(size=len(control)-1)*.001
        _, derivative = model.evaluate(control, True)
        gradients = []
        for coordinate in (0, 1, model.q, 4*model.q, 8*model.q):
            step = 1e-7/(seed['outer']*seed['n']) if coordinate == 8*model.q else 1e-7
            delta = np.zeros_like(control)
            delta[coordinate] = step
            difference = (model.evaluate(control+delta)-model.evaluate(control-delta))/(2*step)
            assert abs(difference-derivative[coordinate]) < 1e-4
            gradients.append({'coordinate': coordinate, 'analytic': float(derivative[coordinate]),
                              'finite_difference': float(difference)})
        records.append({'images': images, 'alphabet': alphabet, 'flag': flag,
                        'original_entries': [[list(a), list(v)] for a,v in sorted(tensor['values'].items())],
                        'parameters': parameters.tolist(), 'construction': model.construction,
                        'integer_profiles': literal_profiles,
                        'gradient_controls': gradients})
    root = Path(__file__).resolve().parents[1]
    names = ['tools/check_f4_pascal_coordinate_sectors.py', 'tools/f4_pascal_coordinate_sectors.py',
             'tools/f4_pascal_flags.py', 'tools/pascal_coordinate_sectors.py',
             'tools/pascal_letter_gauges.py', 'tools/permuted_pascal_kernels.py',
             'tools/cube_f4_bivariate_smith.py', 'tools/weighted_smith_oracle_f4.py',
             'tools/native_matrix_smith_f4.py', 'tools/windowed_f4_smith.py']
    report = {'status': 'pass', 'independent_original_entries_and_endpoints': entries_checked,
              'independent_nonbinary_complete_orders': orders_checked,
              'binary_specialization_complete_orders': binary_orders, 'controls': records,
              'source_sha256': {n:sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope':'Complete selected original-polynomial profiles over F4, with arbitrary '
              'symbol gauges and both physical costs. No stronger contraction or factoring score.'}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('status','independent_original_entries_and_endpoints',
        'independent_nonbinary_complete_orders','binary_specialization_complete_orders')}))


if __name__ == '__main__':
    main()
