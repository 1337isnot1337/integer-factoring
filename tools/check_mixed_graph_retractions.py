#!/usr/bin/env python3
"""Independent coefficient, endpoint, profile and pressure controls for split graphs."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np

from check_schur_gaussian_reduction import polynomial_multiply as binary_multiply
from mixed_graph_retractions import (WORDS, FixedProfileSupport, ambient_source,
    deflated_kernel, frobenius, graph_split, natural, sparse_neighborhood, tensor)
from physical_split_crossings import flatten
from scaled_asymmetric_crossing_pressure import ScaledAsymmetricPressure
from search_f4_recentered_crossings import five_letter_source
from weighted_smith_oracle_f4 import WeightedSmithOracle, polynomial_multiply, scalar_multiply
from windowed_f4_smith import pack


def determinant(matrix):
    """Permutation subset DP in characteristic two, independent of Pfaffian recursion."""
    n = len(matrix)
    dp = [0] * (1 << n)
    dp[0] = 1
    for mask in range(1 << n):
        i = mask.bit_count()
        if i == n or not dp[mask]:
            continue
        for j in range(n):
            if not mask >> j & 1 and matrix[i][j]:
                dp[mask | 1 << j] ^= binary_multiply(dp[mask], matrix[i][j])
    return dp[-1]


def direct_contraction(source, coefficients):
    """Literal sum over all four ambient indices, with packed F4 arithmetic."""
    lx, px = graph_split(coefficients[:5], coefficients[5:10])
    ly, py = graph_split(coefficients[10:15], coefficients[15:])
    maps = (lx, ly, py, px)
    result = []
    for v, u, y, x in product(range(5), repeat=4):
        terms = [[(j, c) for j, c in enumerate(matrix[i]) if c]
                 for matrix, i in zip(maps, (x, y, u, v))]
        value = 0
        for term in product(*terms):
            indices, scalars = zip(*term)
            coefficient = 1
            for c in scalars:
                coefficient = polynomial_multiply(coefficient, c)
            index = sum(j * 6**k for k, j in enumerate(indices))
            value ^= scalar_multiply(pack(source['values'][index]), coefficient)
        result.append(value)
    return result


def independent_profiles(seed, parameters):
    answer = []
    for side, cut in enumerate((1, 0)):
        h = [[0] + [int(v) for v in parameters[side*16+leg*4:side*16+(leg+1)*4]]
             for leg in range(4)]
        matrix, rows, columns = flatten(seed, cut)
        row_weights = [h[0][a] + h[2 if cut else 1][b] for a, b in rows]
        col_weights = [h[1 if cut else 2][a] + h[3][b] for a, b in columns]
        oracle = WeightedSmithOracle([[pack(v) for v in row] for row in matrix], infinity=bool(cut))
        orders, _ = oracle.profile(row_weights, col_weights, int(parameters[-1]))
        answer.extend(orders)
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    source, kernel = ambient_source(), deflated_kernel()
    old = five_letter_source()[0]
    assert len(source['values']) == 1296
    reference = tensor((0,) * 20)
    assert tuple(reference['values']) == old
    determinants = 0
    for index, (value, high) in enumerate(source['values']):
        coordinates = tuple(index // 6**leg % 6 for leg in range(4))
        mask = sum(WORDS[c] << (4*leg) for leg, c in enumerate(coordinates))
        selected = [i for i in range(16) if mask >> i & 1]
        k = len(selected)
        if not k:
            assert value == 1 << 8 and high == 0
            continue
        matrix = [[sum(1 << (p+8) for p in kernel[i][j]) for j in selected] for i in selected]
        square = binary_multiply(value, value)
        assert determinant(matrix) == (square << (8*k-16) if k >= 2 else 0)
        if k == 1:
            assert square == 0
        determinants += 1
    reps = sparse_neighborhood()
    rng = np.random.default_rng(60520261101)
    examples = [(0,) * 20,
                (2, 0, 0, 0, 0, 0, 3, 0, 0, 0) + (0,) * 10,
                tuple(int(x) for x in rng.integers(0, 4, 20))]
    profiles, derivatives, fixed_controls, contractions, conjugates = [], [], [], 0, 0
    for example_index, coefficients in enumerate(examples):
        seed = tensor(coefficients)
        assert [pack(v) for v in seed['values']] == direct_contraction(source, coefficients)
        contractions += 625
        conjugate = tensor(frobenius(coefficients))
        assert conjugate['values'] == [(low ^ high, high) for low, high in seed['values']]
        conjugates += 625
        support = ScaledAsymmetricPressure(seed)
        unit = natural()
        trial = unit.astype(int)
        trial[:-1] += rng.integers(-3, 4, len(trial)-1)
        native = support.base.oracle.prefixes(trial, 1)['orders']
        independent = independent_profiles(seed, trial)
        assert native == independent
        profiles.append({'coefficients': coefficients, 'parameters': trial.tolist(),
                         'independent_orders': independent})
        # Native fixed heights and the closed unit-profile formula must agree.
        unit_profiles = support.base.oracle.prefixes(unit, 1)['orders']
        fixed = FixedProfileSupport(np.asarray(unit_profiles).reshape(2, 25),
                                    (0, 10, 6, 2, 8), (0, 10, 6, 2, 8))
        ordinary = support.base.oracle.types(unit*.21, 1)[0]
        balanced = fixed.evaluate((math.log(.21/math.sqrt(2)), 0.), False)
        assert abs(balanced-math.sqrt(2)*ordinary) < 1e-9
        for logt, logtan in ((math.log(.14), -.7), (math.log(.43), 0.), (math.log(.2), .9)):
            variables = np.r_[unit[:-1], logt, logtan]
            rate, derivative, _ = support.evaluate(variables)
            rate2, derivative2 = fixed.evaluate((logt, logtan))
            assert abs(rate-rate2) < 1e-10
            assert np.max(np.abs(derivative[-2:]-derivative2)) < 1e-9
            fixed_controls.append({'example': example_index, 'log_scale': logt,
                                   'log_tangent': logtan, 'rate': rate})
        point = np.r_[unit[:-1] + rng.normal(0, .07, 32), math.log(.14), .2]
        _, gradient, _ = support.evaluate(point)
        for coordinate in range(34):
            delta = np.zeros(34)
            delta[coordinate] = 1e-6
            finite = (support.evaluate(point+delta, False)[0] - support.evaluate(point-delta, False)[0])/2e-6
            error = abs(float(finite-gradient[coordinate]))
            assert error < 2e-5, (example_index, coordinate, finite, gradient[coordinate])
            derivatives.append({'example': example_index, 'coordinate': coordinate, 'error': error})
    names = ['mixed_graph_retractions.py', 'check_mixed_graph_retractions.py',
             'physical_split_crossings.py', 'reflection_dual_rail_gaussian.py',
             'scaled_asymmetric_crossing_pressure.py', 'asymmetric_crossing_pressure.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
             'weighted_smith_oracle_f4.py']
    report = {'status': 'pass', 'ambient_endpoints': 1296, 'reference_entries': 625,
              'ambient_degree': source['degree'], 'independent_determinants': determinants,
              'direct_packed_contraction_entries': contractions, 'frobenius_entries': conjugates,
              'independent_untruncated_orders': sum(len(p['independent_orders']) for p in profiles),
              'profiles': profiles, 'derivatives': derivatives, 'balanced_normalizations': 3,
              'fixed_profile_controls': fixed_controls, 'sparse_representatives': len(reps),
              'seconds': time.monotonic()-started,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names},
              'scope': 'Exact finite controls. Floating pressures and finite differences are numerical checks, '
                       'not a strict contraction, a global minimum or a factoring improvement.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ('status', 'ambient_degree', 'independent_determinants',
          'direct_packed_contraction_entries', 'independent_untruncated_orders', 'seconds')}))


if __name__ == '__main__':
    main()
