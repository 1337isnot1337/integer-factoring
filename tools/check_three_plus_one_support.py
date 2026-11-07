#!/usr/bin/env python3
"""Numerical support diagnostic for a 3+1/1+3 two-cut profile.

The q4 endpoint tensor is cut as (x,y,u)|(v) and (x)|(y,u,v).  Both cut
matrices are rectangular 64x4 and have rank four.  The existing asymmetric
support equation is applied to the four Smith channels at the zero and
infinity ends, with zero physical letter heights.  The q4 rectangular
2+2 cuts are replayed as a calibration: their support is above one, while
the 3+1 profile is below one at small intrinsic scale.  This is a candidate
diagnostic, not a transfer theorem or a complete-bound claim.
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
    zero_scale = {}
    for a in (1 / math.sqrt(2), 0.2):
        b = math.sqrt(1 - a * a)
        zero_scale[str(a)] = {'rectangular': a + b, 'three_plus_one': (a + b) / 2}
    payload = {
        'scope': ('Numerical support diagnostic for the q4 3+1/1+3 cut; '
                  'no transfer theorem, endpoint invoice or complete-bound claim'),
        'alphabet': q,
        'cut_ranks': {'rectangular_2_plus_2': q * q, 'three_plus_one': q},
        'zero_scale_identity': zero_scale,
        'cases': cases,
        'small_scale_candidate_minimum': min(c['three_plus_one_support'] for c in cases),
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
