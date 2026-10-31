#!/usr/bin/env python3
"""Exhaust a finite F4 polynomial neighborhood of a compressed binary seed.

Each bit direction adds either 1 or omega times w^j(1+w) to one
tensor entry. Thus every candidate has exactly the strict endpoint.
The fixed height ray is screened by complete integer Smith profiles;
scale and free-gauge refinements are floating exploration only.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

from native_rectangular_smith_f4 import RectangularOracle
from screen_full_alphabet_monomial_pressure import optimize
from search_binary_affine_crossings import critical, rate


BASE = (14, 30, 30, 60, 15, 15, 31, 30, 15, 31, 15, 30, 15, 15, 15, 14)
HEIGHTS = (0, 0, -1, -1, -1, -1, 0, 0)


def candidate(directions, degree=5):
    values = [[value, 0] for value in BASE]
    for direction in directions:
        entry, remainder = divmod(direction, 2 * degree)
        power, plane = divmod(remainder, 2)
        values[entry][plane] ^= 3 << power
    return {'values': values, 'q_x': 2, 'q_u': 2, 'dimensions': [2] * 4,
            'denominator_degree': 0, 'strict_endpoint': True}


def main():
    # SciPy is optional for the exact checkers that import this module.
    from scipy.optimize import minimize_scalar
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--weight', type=int, choices=(1, 2, 3, 4), default=3)
    parser.add_argument('--degree', type=int, default=5)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--top', type=int, default=100)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert 5 <= args.degree <= 12
    start = time.monotonic()
    deadline = start + args.seconds
    directions = 32 * args.degree
    total = math.comb(directions, args.weight)
    scales = (.1, .2, .29277969091981376, .4, .6, 1.)
    energies = (-2, -2)
    classes = {}
    completed = 0
    unresolved = []
    result = {'status': 'Finite exact polynomial neighborhood at one height ray; exploratory rate refinements.',
              'base_coefficients': BASE, 'directions': directions, 'hamming_weight': args.weight,
              'maximum_degree': args.degree, 'height_ray': HEIGHTS, 'scales': scales,
              'planned_candidates': total, 'optimized': [],
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    def save():
        result.update(completed_candidates=completed, complete=completed == total and not unresolved,
                      unresolved=unresolved, profile_classes=len(classes), seconds=time.monotonic() - start)
        result['classes'] = sorted(classes.values(), key=lambda row: row['critical_grid'])
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(result, indent=2) + '\n')
        temporary.replace(args.output)
    for moves in itertools.combinations(range(directions), args.weight):
        if time.monotonic() >= deadline:
            break
        seed = candidate(moves, args.degree)
        oracle = RectangularOracle(seed['values'], 2, 2, endpoint_check=True)
        try:
            _, _, orders = oracle.evaluate(list(HEIGHTS) + [1.], (.5, .5), (.5, .5), 1, True)
            assert all(order == int(order) for order in orders)
            profiles = (tuple(map(int, orders[:4])), tuple(map(int, orders[4:])))
        except ArithmeticError as error:
            unresolved.append({'directions': moves, 'reason': str(error)})
            completed += 1
            continue
        if profiles not in classes:
            pressures = [critical(profiles, energies, scale) for scale in scales]
            classes[profiles] = {'directions': moves, 'profiles': profiles, 'count': 0,
                                  'critical_grid': min(pressures),
                                  'best_grid_scale': scales[min(range(len(scales)), key=pressures.__getitem__)]}
        classes[profiles]['count'] += 1
        completed += 1
        if completed % 10000 == 0:
            save()
            print(json.dumps({'completed': completed, 'classes': len(classes),
                              'best_critical': result['classes'][0]['critical_grid'],
                              'seconds': result['seconds']}), flush=True)
    save()
    rng = np.random.default_rng(310772 + args.weight)
    for row in result['classes'][:args.top]:
        if time.monotonic() >= deadline:
            break
        fixed = minimize_scalar(lambda x: rate(row['profiles'], energies, math.exp(x)),
                                bounds=(-6., 2.), method='bounded')
        scale = math.exp(fixed.x)
        outcome = optimize(candidate(row['directions'], args.degree), cut=1, rng=rng,
                           restarts=5, maxfun=600, height_bound=48,
                           deadline=min(deadline, time.monotonic() + 30),
                           starts=[[height * scale for height in HEIGHTS] + [scale]])
        result['optimized'].append({'directions': row['directions'], 'fixed_ray_rho': float(fixed.fun),
                                    'optimization': outcome})
        save()


if __name__ == '__main__':
    main()
