#!/usr/bin/env python3
"""Finite greedy row-target paths for strict identity-endpoint qutrit operators."""
from __future__ import annotations

import argparse
from itertools import product
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

from prime_deficit_operator_synthesis import coefficient_system, synthesize
from prime_deficit_identity_operators import zero_profiles, infinity_profiles, Pressure
from coupled_identity_operator_filters import identity_kernel
from search_qutrit_deficit_operators import encode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', type=int, default=27)
    parser.add_argument('--selected', default='0,13')
    parser.add_argument('--powers', choices=('unit', 'parent'), default='unit')
    parser.add_argument('--infinity-flag', choices=('identity', 'pascal'), default='identity')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    selected = tuple(int(value) for value in args.selected.split(','))
    if len(selected) != 2:
        raise ValueError('This finite driver selects two qutrit directions.')
    height = 1 if args.powers == 'unit' else args.parent
    if args.infinity_flag != 'identity' and args.powers != 'parent':
        raise ValueError('The nontrivial flag screen uses the parent pole for its proved starting kernel.')
    exponents = (-1, 0, 1) if height == 1 else (-height, -1, 0, 1, height)
    system = coefficient_system(3, args.parent, selected, 1, -height, height, exponents)
    data, n = system['data'], args.parent-1
    zero0, infinity0 = data['zero_orders'], data['infinity_orders']
    basis = None
    if args.infinity_flag == 'pascal':
        basis = [[int(i == j) for j in range(9)] for i in range(9)]
        middle = [i for i, state in enumerate(data['states']) if sum(state) == 2]
        for i, row in enumerate(((1, 1, 1), (0, 1, 2), (0, 0, 1))):
            for j, value in enumerate(row):
                basis[middle[i]][middle[j]] = value
        universal_zero = tuple(min(min(value) for source in data['sources'].values()
                                   for (row, column), value in source.items() if row == i)-height
                               for i in range(len(zero0)))
    orders = (tuple(sorted(range(len(zero0)), key=lambda i: (zero0[i], i))),
              tuple(sorted(range(len(zero0)), key=lambda i: (-zero0[i], i))))
    grid = tuple(product(math.log(2*n)/n*np.array((.25, .5, .75, 1., 1.5, 2.)),
                         sorted(set((-2*n, -n, -4, -2, -1, 0, 1, 2, 4, n/2, n, 2*n)))))
    scales, charges = np.array([x[0] for x in grid]), np.array([x[1] for x in grid])
    profiles, kernels, cache = {}, set(), {}
    calls = feasible = completed = 0
    start = time.monotonic()

    def solve(zero, infinity):
        nonlocal calls, feasible
        key = tuple(zero), tuple(infinity)
        if key not in cache:
            calls += 1
            value = synthesize(3, args.parent, selected, 1, -height, height, *key,
                               exponents=exponents, infinity_basis=basis)
            cache[key] = value
            feasible += value is not None
        return cache[key]

    def record(witness):
        if witness is None:
            return
        kernel = witness['kernel']
        key = tuple((a, b, tuple(sorted(value.items()))) for (a, b), value in sorted(kernel.items()))
        if key in kernels:
            return
        kernels.add(key)
        options = (3, args.parent, 1, 0, 2, 1)
        zero = zero_profiles(*options, kernel, selected_indices=selected)
        infinity = infinity_profiles(*options, kernel, selected_indices=selected)
        profile_key = tuple((b['charge'], b['orders']) for b in zero), infinity
        if profile_key in profiles:
            profiles[profile_key]['witness_count'] += 1
            return
        pressure = Pressure(*options, zero, infinity, selected_indices=selected)
        values, relative = pressure.grid(scales, charges), pressure.grid(scales, charges, True)
        i, j = int(np.argmin(values)), int(np.argmin(relative))
        lower = min(power for value in kernel.values() for power in value)
        upper = max(power for value in kernel.values() for power in value)
        profiles[profile_key] = {'witness_count': 1, 'kernel': encode(kernel),
            'zero_blocks': zero, 'infinity_active': infinity,
            'zero_targets': witness['zero_targets'], 'infinity_targets': witness['infinity_targets'],
            'affine_rank': witness['affine_rank'], 'constraints': witness['constraints'],
            'critical_grid': {'pressure': float(values[i]), 'scale': float(scales[i]), 'charge': float(charges[i])},
            'relative_grid': {'change': float(relative[j]), 'scale': float(scales[j]), 'charge': float(charges[j])},
            'common_Laurent_clearing': max(0, -lower),
            'full_polynomial_degree_upper_bound': 8*n*n+8*n+upper+max(0, -lower)}

    if basis is None:
        record(solve(zero0, infinity0))
    else:
        record({'kernel': identity_kernel(3, 2), 'zero_targets': zero0,
                'infinity_targets': tuple(min(infinity0[j] for j in range(9) if basis[i][j]) for i in range(9)),
                'affine_rank': None, 'constraints': 0})
    gains = (1,) if height == 1 or basis is not None else (1, height)
    selected_rows = range(len(infinity0)) if basis is None else (-1, *middle)
    for gain, index, priority in product(gains, selected_rows, orders):
        zero = [value-gain for value in zero0] if basis is None else list(universal_zero)
        infinity = list(infinity0)
        if index >= 0:
            infinity[index] += 2*gain
        witness = solve(zero, infinity)
        assert witness is not None
        for row in priority:
            lo, hi = zero[row], zero0[row]+1
            while lo < hi:
                mid = (lo+hi+1)//2
                proposal = list(zero)
                proposal[row] = mid
                value = solve(proposal, infinity)
                if value is None:
                    hi = mid-1
                else:
                    lo, witness = mid, value
            zero[row] = lo
        record(witness)
        completed += 1
        print(json.dumps({'completed_paths': completed, 'affine_queries': calls,
                          'distinct_kernels': len(kernels), 'exact_profiles': len(profiles),
                          'seconds': time.monotonic()-start}), flush=True)
    # A complementary direction grants an infinity loss and raises one zero row.
    for index, priority in product(range(len(zero0)) if basis is None else (), (tuple(range(len(infinity0))), tuple(reversed(range(len(infinity0)))))):
        zero = list(zero0)
        zero[index] += 1
        infinity = [value-2 for value in infinity0]
        witness = solve(zero, infinity)
        assert witness is not None
        for row in priority:
            proposal = list(infinity)
            proposal[row] = infinity0[row]
            value = solve(zero, proposal)
            if value is not None:
                infinity[row], witness = proposal[row], value
        record(witness)
        completed += 1
    records = sorted(profiles.values(), key=lambda row: row['critical_grid']['pressure'])
    for row in records[:20]:
        pressure = Pressure(3, args.parent, 1, 0, 2, 1, row['zero_blocks'], row['infinity_active'], selected_indices=selected)
        point = row['critical_grid']
        fit = minimize(lambda value: float(pressure.grid([math.exp(float(value[0]))], [float(value[1])])[0]),
                       [math.log(point['scale']), point['charge']], method='Nelder-Mead',
                       bounds=[(-math.log(n)-4, math.log(1.5)), (-2*n, 2*n)],
                       options={'maxiter': 200, 'fatol': 1e-10})
        row['refined'] = {'pressure': float(fit.fun), 'scale': math.exp(float(fit.x[0])),
                          'charge': float(fit.x[1]), 'success': bool(fit.success)}
    result = {'status': 'Completed specified greedy row-target paths over F3; exact affine witnesses and profiles, not exhaustive over polynomial kernels or Smith lattices.',
              'prime': 3, 'parent': args.parent, 'n': n, 'selected': selected,
              'coefficient_exponents': exponents, 'coefficient_variables': len(system['variables']),
              'infinity_row_flag': args.infinity_flag, 'infinity_basis': basis,
              'source_known_below': data['known_below'], 'completed_paths': completed,
              'affine_queries': calls, 'feasible_queries': feasible, 'distinct_witness_kernels': len(kernels),
              'exact_profile_classes': len(records), 'grid_size': len(grid),
              'best_critical_pressure': min(row.get('refined', row['critical_grid'])['pressure'] for row in records),
              'best_relative_grid_change': min(row['relative_grid']['change'] for row in records),
              'records': records, 'seconds': time.monotonic()-start}
    args.output.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'records'}), flush=True)


if __name__ == '__main__':
    main()
