#!/usr/bin/env python3
"""Map the numerical negative band of the H=39 Jensen-tail test."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np

Z_GRID = np.linspace(0.998, 0.99999, 300)
A_RANGE = range(50, 261)
B_RANGE = range(70, 371)
RHO_FLOOR = Fraction(705, 1000)
RHO_CEILING = Fraction(7072, 10000)
NEGATIVE_MARGIN = -1e-9
EXACT_ANCHORS = ((70, 99, 899, 900), (113, 160, 4997, 5000))


def q_coefficients(cap: int) -> Counter[int]:
    q = Counter({0: 1})
    for odd in range(1, cap + 1, 2):
        for height, count in list(q.items()):
            if height + odd <= cap:
                q[height + odd] += count
    return q


def exact_excludes(orders: Counter[int], q: Counter[int], a: int, b: int, z: Fraction) -> bool:
    left = sum(count * z ** (a * order) for order, count in orders.items())
    right = sum(count * z ** (b * height) for height, count in q.items())
    return left ** b > right ** (2 * a)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    core = json.loads(Path('evidence/height24-core-low-orders.json').read_text())
    assert core['completed'] is True and len(core['words']) == 81
    profile = core['profile']
    retained = profile['cuts'][0]['retained_orders_t']
    assert retained == profile['cuts'][1]['retained_orders_t']
    assert len(retained) == 3359 and sum(retained) == 173663
    assert all(cut['tail_count'] == 3202 for cut in profile['cuts'])
    orders = Counter(retained)
    orders[82] += 3202
    assert sum(orders.values()) == 6561

    q = q_coefficients(39)
    assert sum(q.values()) == 450
    ords = np.array(sorted(orders), dtype=float)
    cnts = np.array([orders[int(order)] for order in ords], dtype=float)
    heights = np.array(sorted(q), dtype=float)
    counts = np.array([q[int(height)] for height in heights], dtype=float)
    logz = np.log(Z_GRID)

    rows = []
    with np.errstate(over='ignore', invalid='ignore'):
        for a in A_RANGE:
            for b in B_RANGE:
                rho = Fraction(a, b)
                if rho < RHO_FLOOR or rho > RHO_CEILING:
                    continue
                if 2 * a * a >= b * b:
                    continue
                left = np.sum(
                    cnts[None, :] * np.exp(a * ords[None, :] * logz[:, None]),
                    axis=1,
                )
                right = np.sum(
                    counts[None, :] * np.exp(b * heights[None, :] * logz[:, None]),
                    axis=1,
                )
                log_ratio = b * np.log(left) - 2 * a * np.log(right)
                negative = np.where(log_ratio < NEGATIVE_MARGIN)[0]
                if not len(negative):
                    continue
                rows.append({
                    'a': a,
                    'b': b,
                    'rho': [a, b],
                    'rho_decimal': float(rho),
                    'z_interval_numeric': [float(Z_GRID[negative[0]]), float(Z_GRID[negative[-1]])],
                    'minimum_log_ratio_numeric': float(log_ratio[negative].min()),
                    'negative_grid_points': int(len(negative)),
                })

    assert rows, 'no negative points found'
    rows.sort(key=lambda row: row['rho_decimal'])
    best = rows[0]
    exact = []
    for a, b, z_num, z_den in EXACT_ANCHORS:
        z = Fraction(z_num, z_den)
        assert 0 < z < 1
        assert 2 * a * a < b * b
        excludes = exact_excludes(orders, q, a, b, z)
        assert excludes is False
        exact.append({
            'a': a,
            'b': b,
            'rho': [a, b],
            'rho_decimal': float(Fraction(a, b)),
            'z': [z_num, z_den],
            'test_excludes': excludes,
        })
    assert any(row['a'] == 113 and row['b'] == 160 for row in exact)

    receipt = {
        'completed': True,
        'maximum_height': 39,
        'ambient_height_inventory_q': sum(q.values()),
        'grid': {
            'a': [min(A_RANGE), max(A_RANGE)],
            'b': [min(B_RANGE), max(B_RANGE)],
            'rho_floor': [RHO_FLOOR.numerator, RHO_FLOOR.denominator],
            'rho_ceiling': [RHO_CEILING.numerator, RHO_CEILING.denominator],
            'z': [float(Z_GRID[0]), float(Z_GRID[-1])],
            'z_points': len(Z_GRID),
            'negative_margin': NEGATIVE_MARGIN,
        },
        'negative_pair_count': len(rows),
        'numeric_min_rho': best,
        'numeric_max_rho': rows[-1],
        'exact_anchors': exact,
        'scope': (
            'A numerical map of the H=39 Jensen-tail necessary condition plus two '
            'exact rational anchor points. The map is not a certificate of the '
            'numerical band boundaries, and the anchors are not a split or a score.'
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({
        'completed': True,
        'negative_pairs': len(rows),
        'numeric_min_rho': best['rho_decimal'],
        'exact_anchors': [[row['rho'][0], row['rho'][1]] for row in exact],
    }))


if __name__ == '__main__':
    main()
