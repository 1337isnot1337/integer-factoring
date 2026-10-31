#!/usr/bin/env python3
"""Finite coefficient and full-physical-gauge searches of closed auxiliary blocks.

Bare profiles have exact integer orders; numerical minimization and all
physical-gauge results are exploratory and need an independent replay.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path
import random
import time

import numpy as np
from scipy.optimize import minimize

from closed_auxiliary_blocks import product_masks, tensor
from native_rectangular_smith_f4 import RectangularOracle
from search_rectangular_f4_pressure import paired_value_gradient


def bare_minimum(orders, q):
    first, second = (np.asarray(side, dtype=float) for side in orders)
    def evaluate(t):
        result = derivative = 0.
        for side in (first, second):
            weights = -t*side
            offset = float(max(weights))
            probabilities = np.exp(weights-offset)
            total = float(sum(probabilities))
            result += offset+math.log(total)
            derivative -= float(np.dot(side, probabilities))/total
        return result, derivative
    if evaluate(0.)[1] >= 0:
        t = 0.
    else:
        low, high = 0., 1.
        while evaluate(high)[1] < -1e-12 and high < 256:
            high *= 2
        for _ in range(70):
            mid = (low+high)/2
            if evaluate(mid)[1] < 0:
                low = mid
            else:
                high = mid
        t = (low+high)/2
    return {'rho': evaluate(t)[0]/(4*math.log(q)), 'log_scale': t}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--block-length', type=int, default=3)
    parser.add_argument('--physical-words', help='Comma-separated binary word integers; default full alphabet.')
    parser.add_argument('--random-pairs', type=int, default=256)
    parser.add_argument('--catalogue-limit', type=int)
    parser.add_argument('--weighted-classes', type=int, default=40)
    parser.add_argument('--restarts', type=int, default=24)
    parser.add_argument('--maxfun', type=int, default=1200)
    parser.add_argument('--seconds', type=float, default=3600.)
    parser.add_argument('--random-seed', type=int, default=31102673)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    labels = None if args.physical_words is None else tuple(map(int, args.physical_words.split(',')))
    q = (1 << args.block_length) if labels is None else len(labels)
    if q > 32:
        parser.error('This finite driver limits the physical alphabet to 32.')
    start = time.monotonic(); deadline = start+args.seconds
    randomizer = random.Random(args.random_seed)
    generator = np.random.default_rng(args.random_seed)
    pairs = [(0, 0)]+list(combinations_with_replacement(product_masks(), 2))
    pairs += [(randomizer.randrange(65536), randomizer.randrange(65536)) for _ in range(args.random_pairs)]
    if args.catalogue_limit is not None:
        shuffled = pairs[1:]
        randomizer.shuffle(shuffled)
        pairs = (pairs[:1]+shuffled)[:args.catalogue_limit]
    script_root = Path(__file__).resolve().parent
    source_names = (Path(__file__).name, 'closed_auxiliary_blocks.py',
                    'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp',
                    'search_rectangular_f4_pressure.py')
    report = {'scope': 'Finite closed auxiliary-block construction and floating gauge search; no new factoring bound.',
              'arguments': vars(args) | {'output': str(args.output)},
              'source_sha256': {name: sha256((script_root/name).read_bytes()).hexdigest()
                               for name in source_names},
              'catalogue_size': len(pairs), 'seeds': [], 'profile_classes': [],
              'weighted_records': [], 'completed_catalogue': False,
              'completed_weighted_schedule': False}
    classes = {}; errors = Counter(); queries = 0
    best = None

    def save():
        report.update(elapsed_seconds=time.monotonic()-start,
                      native_profile_queries=queries, rejected_evaluations=dict(errors),
                      best_weighted=best)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(report, indent=2)+'\n')
        temporary.replace(args.output)

    class Budget(Exception):
        pass

    try:
        for index, (b, c) in enumerate(pairs):
            if time.monotonic() >= deadline:
                raise Budget()
            data = tensor(args.block_length, b, c, labels)
            oracle = RectangularOracle(data['values'], q, q)
            profiles = []; costs = []
            for cut in (0, 1):
                _, _, orders = oracle.evaluate([0.]*(oracle.parameters-1)+[1.],
                                                [1/q]*q, [1/q]*q, cut, True)
                queries += 1
                assert all(abs(value-round(value)) < 1e-9 for value in orders)
                sides = [tuple(round(v) for v in orders[:oracle.n]),
                         tuple(round(v) for v in orders[oracle.n:])]
                profiles.append(sides)
                costs.append(bare_minimum(sides, q))
            signature = tuple(tuple(side) for profile in profiles for side in profile)
            item = {'index': index, 'mask_b': b, 'mask_c': c,
                    'bare_rho': min(cost['rho'] for cost in costs), 'orientations': costs}
            report['seeds'].append(item)
            if signature not in classes:
                classes[signature] = len(report['profile_classes'])
                report['profile_classes'].append({**item, 'profiles': profiles})
            if index % 25 == 0:
                save()
                print(json.dumps({'phase': 'bare', 'seed': index, 'classes': len(classes),
                                  'best': min(r['bare_rho'] for r in report['profile_classes']),
                                  'seconds': time.monotonic()-start}), flush=True)
        report['completed_catalogue'] = True
        selected = sorted(report['profile_classes'], key=lambda r: r['bare_rho'])[:args.weighted_classes]
        for class_index, item in enumerate(selected):
            data = tensor(args.block_length, item['mask_b'], item['mask_c'], labels)
            oracle = RectangularOracle(data['values'], q, q)
            local_best = None; history = []
            for cut in (0, 1):
                def objective(normalized):
                    nonlocal queries, local_best, best
                    if time.monotonic() >= deadline:
                        raise Budget()
                    scale = normalized[-1]
                    physical = np.r_[normalized[:-1]*scale, scale]
                    queries += 1
                    try:
                        rho, derivative, px, py = paired_value_gradient(oracle, physical, cut)
                    except (ArithmeticError, ValueError, OverflowError) as error:
                        errors[type(error).__name__+': '+str(error)] += 1
                        return 100., np.zeros(oracle.parameters)
                    gradient = np.r_[derivative[:-1]*scale,
                                     derivative[-1]+np.dot(derivative[:-1], normalized[:-1])]
                    if local_best is None or rho < local_best['rho_floating']-1e-10:
                        local_best = {'rho_floating': float(rho), 'cut': cut,
                                      'normalized_parameters': normalized.tolist(),
                                      'parameters': physical.tolist(), 'type_x_v': px,
                                      'type_y_u': py, 'query': queries}
                        if best is None or rho < best['rho_floating']-1e-10:
                            _, _, orders = oracle.evaluate(physical, px, py, cut, True)
                            best = {**local_best, **data,
                                    'scaled_profiles': [orders[:oracle.n], orders[oracle.n:]],
                                    'coefficient_sha256': sha256(json.dumps(data['values'], separators=(',', ':')).encode()).hexdigest()}
                            save()
                    return rho, gradient
                for restart in range(args.restarts):
                    if restart == 0:
                        point = np.r_[np.zeros(oracle.parameters-1),
                                      max(.01, min(8., item['orientations'][cut]['log_scale']))]
                    elif restart % 5 == 4:
                        point = np.r_[generator.normal(0, 3., oracle.parameters-1), .4]
                    else:
                        point = np.array(local_best['normalized_parameters'])
                        point[:-1] += generator.normal(0, (.03, .15, .6, 2.)[restart % 4], oracle.parameters-1)
                        point[-1] *= math.exp(generator.normal(0, .2))
                    point[:-1] = np.clip(point[:-1], -12., 12.)
                    point[-1] = min(8., max(.01, point[-1]))
                    result = minimize(objective, point, method='L-BFGS-B', jac=True,
                                      bounds=[(-12., 12.)]*(oracle.parameters-1)+[(.01, 8.)],
                                      options={'maxfun': args.maxfun, 'maxiter': args.maxfun,
                                               'ftol': 1e-12, 'gtol': 1e-8, 'maxls': 30})
                    history.append({'cut': cut, 'restart': restart, 'evaluations': int(result.nfev),
                                    'success': bool(result.success), 'message': str(result.message)})
                    report['weighted_in_progress'] = {'class_index': class_index, **item,
                                                       'best': local_best, 'restarts': history}
                    save()
            report['weighted_records'].append({**item, **local_best, 'restarts': history})
            report.pop('weighted_in_progress', None)
            save()
            print(json.dumps({'phase': 'weighted', 'class_index': class_index,
                              'rho': local_best['rho_floating'], 'best': best['rho_floating'],
                              'seconds': time.monotonic()-start}), flush=True)
        report['completed_weighted_schedule'] = True
    except Budget:
        report['time_budget_reached'] = True
    save()


if __name__ == '__main__':
    main()
