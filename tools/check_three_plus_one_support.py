#!/usr/bin/env python3
"""Diagnostic for a 3+1/1+3 alternative leg cut.

The q4 endpoint tensor is cut as (x,y,u)|(v) and (x)|(y,u,v).  Both cut
matrices are rectangular 64x4 and have rank four.  The existing asymmetric
support expression is evaluated on the four Smith channels at the zero and
infinity ends, with zero physical letter heights.  The q4 rectangular 2+2
cuts are replayed as a calibration.

This is not a crossing criterion: under these alternative cuts the endpoint
tensor itself has rank q on both sides, so the trivial one-piece seed has
cost q < q**sqrt(2).  The low value is an artifact of changing the standard
two flattenings, not evidence of a transfer or a complete-bound improvement.
"""

from __future__ import annotations

import argparse
import json
import math
from hashlib import sha256
from pathlib import Path

import numpy as np

from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import tensor

ROOT = Path(__file__).resolve().parents[1]


def cut_matrices(values, q):
    arr = np.empty(len(values), dtype=object)
    arr[:] = values
    a = arr.reshape(q, q, q, q)
    return {
        'rect0': [[a[x, y, u, v] for u in range(q) for v in range(q)]
                  for x in range(q) for y in range(q)],
        'rect1': [[a[x, y, u, v] for y in range(q) for v in range(q)]
                  for x in range(q) for u in range(q)],
        'three_one': [[a[x, y, u, v] for v in range(q)]
                      for x in range(q) for y in range(q) for u in range(q)],
        'one_three': [[a[x, y, u, v] for y in range(q) for u in range(q) for v in range(q)]
                      for x in range(q)],
    }


def profile(matrix, rank, scale):
    n = max(len(matrix), len(matrix[0]))
    square = [[(0, 0)] * n for _ in range(n)]
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            square[i][j] = value
    oracle = MatrixOracle(square, rank, [[0] * n for _ in range(n)])
    zero = oracle.profile([0.] * n, [0.] * n, scale, False)
    infinity = oracle.profile([0.] * n, [0.] * n, scale, True)
    return infinity['orders'], zero['orders']


def endpoint_matrices(q):
    one, zero = (1, 0), (0, 0)
    rect0 = [[one if x == v and y == u else zero
              for u in range(q) for v in range(q)]
             for x in range(q) for y in range(q)]
    rect1 = [[one if x == v and y == u else zero
              for y in range(q) for v in range(q)]
             for x in range(q) for u in range(q)]
    one_three = [[one if y == u and x == v else zero
                  for y in range(q) for u in range(q) for v in range(q)]
                 for x in range(q)]
    return rect0, rect1, one_three


def support(orders_infinity, orders_zero, a, b, q):
    gx = np.zeros(q)
    gy = np.zeros(q)
    raw = (a * float(np.log(np.sum(np.exp(-np.asarray(orders_infinity) / a))))
           + b * float(np.log(np.sum(np.exp(-np.asarray(orders_zero) / b)))))
    constant = raw + min(gx) + min(gy)

    def equation(rate):
        return constant - rate * (
            float(np.log(np.sum(np.exp(-(gx - min(gx)) / rate))))
            + float(np.log(np.sum(np.exp(-(gy - min(gy)) / rate)))))

    low, high = 1e-12, 8.0
    if equation(low) <= 0:
        raise ArithmeticError('support equation has no positive root')
    while equation(high) > 0:
        high *= 2
    for _ in range(100):
        middle = (low + high) / 2
        if equation(middle) > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    seed = tensor(4)
    q = seed['q_x']
    mats = cut_matrices(seed['values'], q)
    cases = []
    for scale in (1e-4, 1e-3, 1e-2):
        rect_inf, rect_zero = profile(mats['rect0'], q * q, scale)
        three_inf, three_zero = profile(mats['three_one'], q, scale)
        # The second cut is the other pairing; the calibration keeps the
        # ordinary 2+2 support, while the candidate uses 1+3.
        rect2_inf, rect2_zero = profile(mats['rect1'], q * q, scale)
        one_inf, one_zero = profile(mats['one_three'], q, scale)
        for name, a in (('balanced', 1 / math.sqrt(2)), ('asymmetric', 0.2)):
            b = math.sqrt(1 - a * a)
            rect_value = support(rect_inf, rect2_zero, a, b, q)
            candidate = support(three_inf, one_zero, a, b, q)
            cases.append({
                'scale': scale, 'angle_control': name,
                'a': a, 'b': b,
                'rectangular_support': rect_value,
                'three_plus_one_support': candidate,
                'rectangular_orders_infinity_4': rect_inf[:4],
                'rectangular_orders_zero_4': rect2_zero[:4],
                'three_plus_one_orders_infinity': three_inf,
                'three_plus_one_orders_zero': one_zero,
            })
            if name == 'balanced':
                assert rect_value > 1.4
                assert candidate < 0.85
            else:
                assert rect_value > 1.1
                assert candidate < 0.75
    endpoint = endpoint_matrices(q)
    endpoint_rect0_inf, _ = profile(endpoint[0], q * q, 1e-4)
    endpoint_rect1_inf, endpoint_rect1_zero = profile(endpoint[1], q * q, 1e-4)
    _, endpoint_one_three_zero = profile(endpoint[2], q, 1e-4)
    endpoint_control = []
    for angle in (0.05, 0.2, 0.5, 0.95):
        b = math.sqrt(1 - angle * angle)
        endpoint_control.append({
            'a': angle, 'b': b,
            'gaussian_std0_13': support(rect_inf, one_zero, angle, b, q),
            'endpoint_std0_13': support(endpoint_rect0_inf, endpoint_one_three_zero, angle, b, q),
            'gaussian_std1_13': support(rect2_inf, one_zero, angle, b, q),
            'endpoint_std1_13': support(endpoint_rect1_inf, endpoint_one_three_zero, angle, b, q),
            'gaussian_31_std1': support(three_inf, rect2_zero, angle, b, q),
            'endpoint_31_std1': support(three_inf, endpoint_rect1_zero, angle, b, q),
            'gaussian_minus_endpoint_std0': (
                support(rect_inf, one_zero, angle, b, q)
                - support(endpoint_rect0_inf, endpoint_one_three_zero, angle, b, q)),
        })
        assert endpoint_control[-1]['gaussian_minus_endpoint_std0'] >= -1e-6

    zero_scale = {}
    for a in (1 / math.sqrt(2), 0.2):
        b = math.sqrt(1 - a * a)
        zero_scale[str(a)] = {'rectangular': a + b, 'three_plus_one': (a + b) / 2}
    payload = {
        'scope': ('Diagnostic for the q4 3+1/1+3 alternative leg-cut pair; not '
                  'a crossing criterion because the one-piece endpoint has cost '
                  'q < q**sqrt(2) for this pair; mixed standard/3+1 pairs are '
                  'untested; no transfer or complete-bound claim'),
        'alphabet': q,
        'cut_ranks': {'rectangular_2_plus_2': q * q, 'three_plus_one': q},
        'zero_scale_identity': zero_scale,
        'cases': cases,
        'small_scale_candidate_minimum': min(c['three_plus_one_support'] for c in cases),
        'mixed_endpoint_control': endpoint_control,
        'trivial_endpoint': {
            'three_plus_one_rank': q, 'one_three_rank': q,
            'rank_sum_cost': q, 'threshold': q ** math.sqrt(2),
            'passes_alternative_criterion': q < q ** math.sqrt(2),
        },
    }
    raw = json.dumps(payload, separators=(',', ':'), sort_keys=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(raw + '\n')
    print(json.dumps({
        'cases': len(cases),
        'small_scale_candidate_minimum': payload['small_scale_candidate_minimum'],
        'sha256': sha256(raw.encode()).hexdigest(),
    }, sort_keys=True))


if __name__ == '__main__':
    main()
