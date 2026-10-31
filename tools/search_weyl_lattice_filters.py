#!/usr/bin/env python3
"""Finite polynomial-filter search through exact reduced row-lattice bounds.

Every reported bound is feasible over F4 and has an explicit endpoint
repair. A returned polynomial witness is also evaluated with its exact
profiles. Neither a finite catalogue nor a non-improving sufficient bound
excludes all polynomials with the same degree caps.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import accumulate, permutations, product
import json
from pathlib import Path
from time import monotonic

from f4_repaired_identity_filters import reduced_zero_profiles
from search_f4_repaired_identity_filters import pressure_grid
from weyl_lattice_filter_synthesis import (
    bound_blocks, coefficient_tables, json_ready, lattice_basis,
    polynomial_order, synthesize,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, choices=(2, 3), default=2)
    parser.add_argument('--geometry', default='mixed',
                        help='A proved geometry, including split:1:2 for three directions.')
    parser.add_argument('--degree-grid', help='Comma-separated nonnegative caps for the first three states.')
    parser.add_argument('--top-degree', type=int, default=0)
    parser.add_argument('--common-offset', type=int, default=0)
    parser.add_argument('--all-degrees', action='store_true')
    parser.add_argument('--start', choices=('minimum', 'gaussian'), default='minimum')
    parser.add_argument('--keep-witnesses', action='store_true')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--priority-limit', type=int)
    parser.add_argument('--degree-family', choices=('occupation', 'split-cold-hot'))
    parser.add_argument('--hot-grid', default='0,1', help='Degree caps for cold-occupied split states.')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    n, ramification, geometry = args.n, args.ramification, args.geometry
    modes, size = args.modes, 2 ** args.modes
    data = lattice_basis(n, ramification, geometry, modes)
    grid = (tuple(map(int, args.degree_grid.split(','))) if args.degree_grid else
            tuple(sorted({0, ramification, n * ramification // 2, n * ramification,
                          n * ramification + 1, 2 * n * ramification,
                          2 * n * ramification + 1, 4 * n * ramification})))
    if args.common_offset < 0 or args.top_degree < 0 or any(value < 0 for value in grid):
        raise ValueError('Use nonnegative polynomial degree caps and offset.')
    if modes == 3 and args.all_degrees:
        raise ValueError('The three-mode driver uses occupation-class degree caps; omit --all-degrees.')
    if modes == 2:
        degrees = tuple(product(grid, repeat=4)) if args.all_degrees else tuple(
            values + (args.top_degree,) for values in product(grid, repeat=3))
    elif (args.degree_family or ('split-cold-hot' if geometry == 'split:1:2' else 'occupation')) == 'occupation':
        degrees = tuple(tuple((values + (args.top_degree,))[sum(state)] for state in data['states'])
                        for values in product(grid, repeat=3))
    else:
        if geometry != 'split:1:2':
            raise ValueError('The cold/hot degree family uses split:1:2.')
        hot_grid = tuple(map(int, args.hot_grid.split(',')))
        degrees = set()
        for base, slope, hot in product(grid, (-1, 0, 1), product(hot_grid, repeat=3)):
            cold_empty = tuple(base + slope * sum(state[1:]) for state in data['states'][:4])
            if min(cold_empty) >= 0:
                degrees.add(cold_empty + hot + (args.top_degree,))
        degrees = tuple(sorted(degrees))
    degrees = tuple(tuple(value + args.common_offset for value in bounds) for bounds in degrees)
    if args.limit:
        degrees = degrees[:args.limit]
    if modes == 2:
        priorities = tuple(permutations(range(size)))
    else:
        orders = (tuple(range(size)), tuple(reversed(range(size))),
                  (0, 1, 4, 2, 5, 3, 6, 7), (7, 6, 3, 5, 2, 4, 1, 0))
        priorities = tuple(dict.fromkeys(order[shift:] + order[:shift]
                                        for order in orders for shift in range(size)))
    if args.priority_limit:
        priorities = priorities[:args.priority_limit]
    scales = sorted({value / n for value in (.1, .3, 1, 2, 3, 5, 8)} | {.125, .25, .5, 1., 2., 4.})
    charges = sorted({1 - 2 * n, -n, -3, -1, 0, 1, 2, 3, n / 2, n, 2 * n})
    evaluate = pressure_grid(n, ramification, geometry, scales, charges, modes)
    minimum_rows = tuple(min(polynomial_order(poly) for basis in data['coefficients']
                             for (row, column), poly in basis.items() if row == i)
                         for i in range(size))
    rows, seen = [], set()
    calls, started = 0, monotonic()
    result = {
        'status': 'Finite exact-feasibility polynomial search; floating pressure is exploratory.',
        'n': n, 'ramification': ramification, 'geometry': geometry, 'modes': modes,
        'coefficient_field': 'F4', 'degree_grid': grid, 'common_offset': args.common_offset,
        'all_state_degree_caps_varied': args.all_degrees, 'top_degree_cap': args.top_degree,
        'degree_cap_family': ('individual states' if modes == 2 else
                             (args.degree_family or ('split-cold-hot' if geometry == 'split:1:2'
                                                     else 'total occupation classes'))),
        'split_hot_degree_grid': args.hot_grid if geometry == 'split:1:2' else None,
        'degree_vectors': len(degrees), 'row_priority_orders': priorities,
        'starting_rows': args.start, 'gaussian_row_orders': data['gaussian_orders'],
        'minimum_row_orders': minimum_rows, 'scale_grid': scales, 'charge_grid': charges,
        'gain_detection_tolerance': 1e-9, 'rows': rows,
    }
    for index, bounds in enumerate(degrees):
        variables, endpoints, equations = coefficient_tables(n, ramification, geometry, bounds, modes)
        maxima = tuple(max((power for power, value in row), default=0) + 1 for row in equations)
        for priority in priorities:
            target = list(minimum_rows if args.start == 'minimum' else
                          tuple(value + min(bounds) for value in data['gaussian_orders']))
            assert synthesize(n, ramification, geometry, bounds, target, construct=False, modes=modes)
            for row in priority:
                low, high = target[row], maxima[row]
                while high - low > 1:
                    middle = (low + high) // 2
                    target[row] = middle
                    calls += 1
                    if synthesize(n, ramification, geometry, bounds, target, construct=False, modes=modes):
                        low = middle
                    else:
                        high = middle
                target[row] = low
            key = bounds, tuple(target)
            if key in seen:
                continue
            seen.add(key)
            certified = bound_blocks(n, ramification, geometry, target, modes)
            critical_bound, relative_bound = evaluate(tuple(-value for value in bounds), certified)
            witness = synthesize(n, ramification, geometry, bounds, target, modes=modes)
            assert witness
            actual = reduced_zero_profiles(n, ramification, modes, witness['alpha'], geometry)
            for block, bound in zip(actual, certified):
                assert block['charge'] == bound['charge']
                assert all(a >= b for a, b in zip(accumulate(block['orders']), accumulate(bound['orders'])))
            critical, relative = evaluate(tuple(-value for value in witness['actual_degrees']), actual)
            portable = json_ready(witness)
            row = {'degree_bounds': bounds, 'row_targets': tuple(target),
                   'critical_bound_grid': critical_bound, 'relative_bound_grid': relative_bound,
                   'actual_degrees': witness['actual_degrees'], 'endpoint_values': witness['endpoint_values'],
                   'endpoint_image_dimension': witness['endpoint_image_dimension'],
                   'actual_critical_grid': critical, 'actual_relative_grid': relative,
                   'exact_zero_blocks_at_J1': [{'charge': block['charge'], 'orders': block['orders']}
                                               for block in actual],
                   'witness_polynomial_sha256': sha256(json.dumps(
                       portable['alpha'], separators=(',', ':')).encode()).hexdigest()}
            if args.keep_witnesses or critical[0] < -1e-9 or relative[0] < -1e-9:
                row['witness'] = portable
            rows.append(row)
            if relative[0] < -1e-9 or critical[0] < -1e-9:
                print('CANDIDATE_GAIN', row, flush=True)
        if (index + 1) % 20 == 0:
            print('PROGRESS', index + 1, len(degrees), calls, len(rows), monotonic() - started, flush=True)
    result['exact_feasibility_queries'] = calls
    result['distinct_degree_row_pairs'] = len(rows)
    result['minimum_sampled_critical_bound'] = min(row['critical_bound_grid'][0] for row in rows)
    result['minimum_sampled_actual_critical_pressure'] = min(row['actual_critical_grid'][0] for row in rows)
    result['minimum_sampled_actual_relative_log_cost'] = min(row['actual_relative_grid'][0] for row in rows)
    result['improving_sufficient_bound_found'] = any(row['relative_bound_grid'][0] < -1e-9 for row in rows)
    result['same_gauge_gain_in_returned_witness'] = any(row['actual_relative_grid'][0] < -1e-9 for row in rows)
    result['negative_critical_pressure_in_returned_witness'] = any(row['actual_critical_grid'][0] < -1e-9 for row in rows)
    result['elapsed_seconds'] = monotonic() - started
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print('DONE', len(rows), calls, result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
