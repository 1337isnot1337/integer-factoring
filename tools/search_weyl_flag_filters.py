#!/usr/bin/env python3
"""Finite polynomial search using constant F4 flags and repaired endpoints."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
from time import monotonic

from f16_diagonal_filter_profiles import zero_profiles, verify_flag_witness
from search_f4_repaired_identity_filters import pressure_grid
from weyl_lattice_filter_synthesis import bound_blocks, json_ready
from weyl_flag_filter_synthesis import coefficient_tables, flag_basis, full_flags, identity_flag, synthesize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, choices=(2, 3), default=2)
    parser.add_argument('--geometry', default='mixed')
    parser.add_argument('--degree-grid', default='0,1,32,33,64,65')
    parser.add_argument('--hot-grid', default='0,1')
    parser.add_argument('--common-offset', type=int, default=0)
    parser.add_argument('--top-degree', type=int, default=0)
    parser.add_argument('--flag-family', choices=('all', 'nonbinary', 'one-sector', 'paired'), default='nonbinary')
    parser.add_argument('--flag-limit', type=int)
    parser.add_argument('--degree-limit', type=int)
    parser.add_argument('--priority-limit', type=int)
    parser.add_argument('--field-size', type=int, choices=(4, 16), default=16)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    n, P, modes, g = args.n, args.ramification, args.modes, args.geometry
    states = tuple(product(range(2), repeat=modes))
    grid = tuple(map(int, args.degree_grid.split(',')))
    if modes == 2:
        degrees = tuple(values + (args.top_degree,) for values in product(grid, repeat=3))
        flags = tuple((0, index, 0) for index, flag in enumerate(full_flags(2))
                      if args.flag_family != 'nonbinary' or any(value > 1 for row in flag for value in row))
        priorities = tuple(permutations(range(4)))
    else:
        if g == 'split:1:2':
            hot = tuple(map(int, args.hot_grid.split(',')))
            degrees = tuple(sorted({(base,) * 4 + tail + (args.top_degree,)
                                    for base in grid for tail in product(hot, repeat=3)}))
        else:
            degrees = tuple(tuple((values + (args.top_degree,))[sum(state)] for state in states)
                            for values in product(grid, repeat=3))
        identity = identity_flag(3)
        if args.flag_family == 'paired':
            flags = tuple((0, first, second, 0) for first, second in product(range(105), repeat=2))
        else:
            flags = tuple(dict.fromkeys([(0, index, identity, 0) for index in range(105)] +
                                        [(0, identity, index, 0) for index in range(105)]))
            if args.flag_family == 'nonbinary':
                flags = tuple(value for value in flags if any(
                    entry > 1 for index in value[1:3] for row in full_flags(3)[index] for entry in row))
        orders = (tuple(range(8)), tuple(reversed(range(8))),
                  (0, 1, 4, 2, 5, 3, 6, 7), (7, 6, 3, 5, 2, 4, 1, 0))
        priorities = tuple(dict.fromkeys(order[shift:] + order[:shift]
                                        for order in orders for shift in range(8)))
    degrees = tuple(tuple(value + args.common_offset for value in caps) for caps in degrees)
    if args.flag_limit:
        flags = flags[:args.flag_limit]
    if args.degree_limit:
        degrees = degrees[:args.degree_limit]
    if args.priority_limit:
        priorities = priorities[:args.priority_limit]
    scales = sorted({value / n for value in (.1, .3, 1, 2, 3, 5, 8)} | {.125, .25, .5, 1., 2., 4.})
    charges = sorted({1 - 2 * n, -n, -3, -1, 0, 1, 2, 3, n / 2, n, 2 * n})
    evaluate = pressure_grid(n, P, g, scales, charges, modes)
    rows, seen, actual_seen = [], set(), set()
    calls = 0
    started = monotonic()
    for flag_index, flag in enumerate(flags):
        data = flag_basis(n, P, g, modes, flag)
        for degree_index, caps in enumerate(degrees):
            _, _, equations = coefficient_tables(n, P, g, modes, flag, caps)
            maxima = tuple(max((power for power, value in row), default=0) + 1 for row in equations)
            for priority in priorities:
                target = list(data['minimum_rows'])
                assert synthesize(n, P, g, caps, target, flag, modes, False, args.field_size)
                for index in priority:
                    low, high = target[index], maxima[index]
                    while high - low > 1:
                        middle = (low + high) // 2
                        target[index] = middle
                        calls += 1
                        if synthesize(n, P, g, caps, target, flag, modes, False, args.field_size):
                            low = middle
                        else:
                            high = middle
                    target[index] = low
                key = flag, caps, tuple(target)
                if key in seen:
                    continue
                seen.add(key)
                bounds = bound_blocks(n, P, g, target, modes)
                cb, rb = evaluate(tuple(-value for value in caps), bounds)
                witness = synthesize(n, P, g, caps, target, flag, modes, True, args.field_size)
                witness['exact_verification'] = verify_flag_witness(witness)
                actual = zero_profiles(n, P, modes, witness['alpha'], g)
                critical, relative = evaluate(tuple(-value for value in witness['actual_degrees']), actual)
                profile_key = witness['actual_degrees'], tuple((b['charge'], b['orders']) for b in actual)
                actual_seen.add(profile_key)
                portable = json_ready(witness)
                row = {'flags': flag, 'degree_bounds': caps, 'row_targets': tuple(target),
                       'critical_bound_grid': cb, 'relative_bound_grid': rb,
                       'actual_degrees': witness['actual_degrees'], 'endpoint_values': witness['endpoint_values'],
                       'endpoint_image_dimension': witness['endpoint_image_dimension'],
                       'actual_critical_grid': critical, 'actual_relative_grid': relative,
                       'exact_zero_blocks_at_J1': actual,
                       'witness_polynomial_sha256': sha256(json.dumps(portable['alpha'], separators=(',', ':')).encode()).hexdigest()}
                if critical[0] < -1e-9 or relative[0] < -1e-9:
                    row['witness'] = portable
                    print('CANDIDATE_GAIN', json.dumps(row), flush=True)
                rows.append(row)
            if (degree_index + 1) % 20 == 0:
                print('PROGRESS', flag_index + 1, len(flags), degree_index + 1, len(degrees), calls,
                      len(rows), len(actual_seen), monotonic() - started, flush=True)
    result = {'status': 'Finite exact F4 flag synthesis with repaired endpoints; floating pressure is exploratory.',
              'n': n, 'ramification': P, 'modes': modes, 'geometry': g,
              'coefficient_field': f'F{args.field_size}', 'row_basis_field': 'F4',
              'flag_family': args.flag_family, 'flags': flags, 'degree_grid': grid,
              'common_offset': args.common_offset, 'top_degree': args.top_degree,
              'degree_vectors': len(degrees), 'row_priority_orders': priorities,
              'scale_grid': scales, 'charge_grid': charges,
              'exact_feasibility_queries': calls, 'returned_degree_flag_row_targets': len(rows),
              'distinct_actual_profile_pairs': len(actual_seen), 'rows': rows,
              'minimum_sampled_critical_pressure': min(row['actual_critical_grid'][0] for row in rows),
              'minimum_sampled_relative_log_cost': min(row['actual_relative_grid'][0] for row in rows),
              'gain_detection_tolerance': 1e-9,
              'negative_critical_pressure_found': any(row['actual_critical_grid'][0] < -1e-9 for row in rows),
              'same_gauge_gain_found': any(row['actual_relative_grid'][0] < -1e-9 for row in rows),
              'elapsed_seconds': monotonic() - started}
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print('DONE', len(rows), len(actual_seen), calls, result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
