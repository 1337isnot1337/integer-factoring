#!/usr/bin/env python3
"""Enumerate binary affine polynomial crossings with a fixed strict endpoint.

Every 16-bit table A gives T(w)=A+w(A+E), where E is the crossing
permutation tensor. The complete finite coefficient family is enumerated;
height rays and subsequent floating optimization are explicitly limited.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from native_rectangular_smith_f4 import RectangularOracle
from screen_full_alphabet_monomial_pressure import optimize


def seed(code):
    values = []
    for packed in range(16):
        x, y, u, v = ((packed >> leg) & 1 for leg in range(4))
        endpoint = int(x == v and y == u)
        a = (code >> packed) & 1
        values.append((a | ((a ^ endpoint) << 1), 0))
    return {'values': values, 'q_x': 2, 'q_u': 2, 'dimensions': [2] * 4,
            'denominator_degree': 0, 'strict_endpoint': True, 'code': code}


def log_partition(values):
    values = list(values)
    ground = min(values)
    values.remove(ground)
    return -ground + math.log1p(math.fsum(math.exp(ground - value) for value in values))


def physical_energies(heights):
    return heights[0] + heights[3] + heights[4] + heights[7], heights[1] + heights[2] + heights[5] + heights[6]


def critical(orders, energies, scale):
    R = math.sqrt(2)
    return sum(log_partition([scale * value for value in profile]) for profile in orders) - R * sum(
        log_partition((0, scale * energy / R)) for energy in energies)


def rate(orders, energies, scale):
    # SciPy is optional for the exact checkers that import this module.
    from scipy.optimize import brentq
    ground = sum(min(profile) for profile in orders)
    physical_ground = sum(min(0, value) for value in energies)
    profiles = [[value - min(profile) for value in profile] for profile in orders]
    physical = [abs(value) for value in energies]
    numerator = (physical_ground - ground) * scale + sum(
        log_partition([scale * value for value in profile]) for profile in profiles)
    def equation(R):
        return numerator - R * sum(math.log1p(math.exp(-scale * energy / R)) for energy in physical)
    if equation(1e-9) <= 0:
        return 2.
    upper = 4.
    while equation(upper) > 0 and upper < 1024:
        upper *= 2
    if upper >= 1024:
        return 2.
    return brentq(equation, 1e-9, upper) / 2


def main():
    # SciPy is optional for the exact checkers that import this module.
    from scipy.optimize import minimize_scalar
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--count', type=int, default=65536)
    parser.add_argument('--random-rays', type=int, default=12)
    parser.add_argument('--top', type=int, default=80)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--seed', type=int, default=310771)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert 0 <= args.start < args.start + args.count <= 65536
    rng = np.random.default_rng(args.seed)
    rays = {(0,) * 8, (0, 0, -1, -1, -1, -1, 0, 0)}
    for J in (-1, 0, 1):
        for K in (-1, 0, 1):
            rays.add((1, -1 + K, 0, K, J, -J, 1, -1))
    for _ in range(args.random_rays):
        rays.add(tuple(map(int, rng.integers(-1, 2, 8))))
    rays = sorted(rays)
    scales = (.08, .16, .3, .6, 1.2, 2.4)
    start = time.monotonic()
    deadline = start + args.seconds
    classes = {}
    completed = queries = 0
    failures = []
    result = {'status': 'Complete requested coefficient range; finite height rays and floating optimization only.',
              'family': 'A+w(A+E), A in F2^16', 'start_code': args.start, 'planned_count': args.count,
              'height_rays': rays, 'scales': scales, 'random_seed': args.seed, 'optimized': [],
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    def save():
        result.update(completed_candidates=completed, complete=completed == args.count and not failures,
                      exact_profile_queries=queries, unresolved=failures, profile_classes=len(classes),
                      seconds=time.monotonic() - start)
        result['best_classes'] = sorted(classes.values(), key=lambda row: row['critical_grid'])[:400]
        temporary = args.output.with_suffix('.tmp')
        temporary.write_text(json.dumps(result, indent=2) + '\n')
        temporary.replace(args.output)
    for code in range(args.start, args.start + args.count):
        if time.monotonic() >= deadline:
            break
        candidate = seed(code)
        oracle = RectangularOracle(candidate['values'], 2, 2, endpoint_check=True)
        for ray in rays:
            try:
                _, _, values = oracle.evaluate(list(ray) + [1.], (.5, .5), (.5, .5), 1, True)
                assert all(value == int(value) for value in values)
                orders = (tuple(map(int, values[:4])), tuple(map(int, values[4:])))
                queries += 1
            except ArithmeticError as error:
                failures.append({'code': code, 'ray': ray, 'error': str(error)})
                continue
            energies = physical_energies(ray)
            key = orders, energies
            if key not in classes:
                pressures = [critical(orders, energies, scale) for scale in scales]
                i = min(range(len(scales)), key=pressures.__getitem__)
                classes[key] = {'code': code, 'ray': ray, 'profiles': orders, 'energies': energies,
                                'critical_grid': pressures[i], 'scale': scales[i], 'count': 0}
            classes[key]['count'] += 1
        completed += 1
        if completed % 1024 == 0:
            save()
            print(json.dumps({'completed': completed, 'classes': len(classes), 'queries': queries,
                              'best_critical': result['best_classes'][0]['critical_grid'],
                              'seconds': result['seconds']}), flush=True)
    save()
    seen = set()
    for row in result['best_classes']:
        if len(result['optimized']) >= args.top or time.monotonic() >= deadline:
            break
        if row['code'] in seen:
            continue
        seen.add(row['code'])
        fixed = minimize_scalar(lambda x: rate(row['profiles'], row['energies'], math.exp(x)),
                                bounds=(-6., 2.), method='bounded')
        scale = math.exp(fixed.x)
        outcome = optimize(seed(row['code']), cut=1, rng=rng, restarts=4, maxfun=500,
                           height_bound=32, deadline=min(deadline, time.monotonic() + 35),
                           starts=[[value * scale for value in row['ray']] + [scale]])
        result['optimized'].append({'code': row['code'], 'fixed_ray_rho': float(fixed.fun),
                                    'fixed_ray_scale': scale, 'optimization': outcome})
        save()


if __name__ == '__main__':
    main()
