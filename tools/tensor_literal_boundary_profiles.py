#!/usr/bin/env python3
"""Literal Pfaffian profiles of homogeneous low-cell deformations.

Every other height cell is retained, including cells above the altered
one. Deletion closure is neither assumed nor used by this evaluator.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import logsumexp

import tensor_hot_compounds_f16 as field
import tensor_literal_hot_pfaffians_f16 as literal
from tensor_hot_boundary_retractions import boundary_code


def pressure_fit(profile):
    zero = np.asarray(profile['zero_orders_t'], dtype=float)
    infinity = np.asarray(profile['infinity_orders_t'], dtype=float)
    heights = np.asarray(profile['heights'], dtype=float)
    queries = 0

    def objective(log_scale):
        nonlocal queries
        queries += 1
        scale = float(np.exp(log_scale))
        target = float(logsumexp(-scale * zero) + logsumexp(-scale * infinity))
        equation = lambda rho: target - 4 * rho * float(logsumexp(-scale * heights / rho))
        return float(brentq(equation, 1e-7, 2, xtol=1e-13))

    grid = np.linspace(np.log(1e-4), np.log(2), 81)
    values = [objective(x) for x in grid]
    candidates = [{'rho': values[i], 'scale_t': float(np.exp(grid[i])), 'source': 'grid'}
                  for i in (0, len(grid)-1)]
    for i in range(1, len(grid)-1):
        if values[i] <= min(values[i-1], values[i+1]):
            fit = minimize_scalar(objective, bounds=(grid[i-1], grid[i+1]), method='bounded',
                                  options={'xatol': 1e-10, 'maxiter': 80})
            candidates.append({'rho': float(fit.fun), 'scale_t': float(np.exp(fit.x)),
                               'source': 'refined local minimum', 'completed': bool(fit.success)})
    best = min(candidates, key=lambda x: x['rho'])
    return {**best, 'queries': queries, 'grid_points': len(grid), 'scale_range': [1e-4, 2],
            'local_candidates': candidates, 'global_minimum_certified': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--height', type=int, default=20)
    parser.add_argument('--boundary-height', type=int, default=15)
    parser.add_argument('--coefficients', type=int, nargs='+', default=[1, 1, 0])
    parser.add_argument('--retraction', choices=('orthogonal', 'coordinate'), default='orthogonal')
    parser.add_argument('--order-cap', type=int, default=64)
    parser.add_argument('--cache-budget', type=int, default=8_000_000)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    modes = (args.height + 1) // 2
    code = boundary_code(args.height, args.boundary_height, args.coefficients, modes=modes,
                         orthogonal=args.retraction == 'orthogonal')
    p, left, heights, q = code['embedding'], code['retraction'], code['heights'], code['q']
    parities = Counter(h % 2 for h in heights)
    sizes = [parities[0]**2 + parities[1]**2, 2*parities[0]*parities[1]]
    assert max(sizes) <= 4608
    names = [Path(__file__).name, 'tensor_literal_hot_pfaffians_f16.py',
             'tensor_literal_hot_pfaffians_f16.cpp', 'tensor_hot_boundary_retractions.py',
             'tensor_hot_compounds_f16_fast.cpp', 'tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16.cpp', 'tensor_orthogonal_scaling.py',
             'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'code': code, 'order_cap': args.order_cap,
              'cache_budget': args.cache_budget, 'completed': False, 'in_progress': None,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started = time.monotonic()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        report['elapsed_seconds'] = time.monotonic() - started
        temporary = args.output.with_suffix(args.output.suffix + '.tmp')
        temporary.write_text(json.dumps(report, separators=(',', ':')) + '\n')
        temporary.replace(args.output)

    print(json.dumps({'starting': args.coefficients, 'height': args.height, 'q': q,
                      'parity_blocks': sizes, 'retraction_kind': args.retraction,
                      'deletion_closed': {k: v['deletion_closed'] for k, v in code['deletion_diagnostics'].items()}}), flush=True)
    cuts = []
    for cut in (('zero',) if args.retraction == 'orthogonal' else ('zero', 'infinity')):
        orders, blocks = [], []
        for parity in (0, 1):
            pairs = [(i, j) for i in range(q) for j in range(q) if (heights[i]+heights[j]) % 2 == parity]
            report['in_progress'] = {'cut': cut, 'parity': parity, 'block_size': len(pairs)}
            save()
            rows = [field.tensor_vector(left[i], left[j] if cut == 'zero' else p[j], modes) for i, j in pairs]
            columns = [field.tensor_vector(p[j] if cut == 'zero' else left[j], p[i], modes) for i, j in pairs]
            before = time.monotonic()
            result = literal.compound_block(modes, rows, columns, order_cap=args.order_cap,
                                            cache_budget=args.cache_budget)
            orders.extend(result['orders_t'])
            blocks.append({key: value for key, value in result.items() if key != 'orders_t'})
            blocks[-1].update(parity=parity, seconds=time.monotonic()-before)
        assert len(orders) == q*q
        cuts.append({'cut': cut, 'orders': sorted(orders), 'blocks': blocks})
    zero = cuts[0]['orders']
    infinity = zero[:] if args.retraction == 'orthogonal' else cuts[1]['orders']
    data = {'q': q, 'modes': modes, 'heights': heights, 'zero_orders_t': zero,
            'infinity_orders_t': infinity, 'order_cap': args.order_cap,
            'computed_cuts': [x['cut'] for x in cuts],
            'infinity_profile_source': ('literal orthogonal role-swap identity; not separately computed'
                                       if args.retraction == 'orthogonal' else 'directly computed'),
            'blocks': [{'cut': x['cut'], 'blocks': x['blocks']} for x in cuts],
            'seconds': time.monotonic()-started}
    lower = {**data, 'zero_orders_t': [x for x in zero if x < args.order_cap],
             'infinity_orders_t': [x for x in infinity if x < args.order_cap]}
    fits = {'upper_counting_bound_fit': pressure_fit(data), 'lower_counting_bound_fit': pressure_fit(lower)}
    data['zero_counts'] = sorted(Counter(data.pop('zero_orders_t')).items())
    data['infinity_counts'] = sorted(Counter(data.pop('infinity_orders_t')).items())
    report.update(completed=True, in_progress=None, profile=data, pressure_fits=fits)
    save()
    print(json.dumps({'completed': True, 'q': q, 'fits': fits, 'seconds': data['seconds']}), flush=True)


if __name__ == '__main__':
    main()
