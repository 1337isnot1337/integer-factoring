#!/usr/bin/env python3
"""Search joint beta/degree Pareto targets with residue-compressed equations."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from time import monotonic

from joint_binary_filter_synthesis import synthesize, local_endpoint_image
from prime_power_identity_filters import finite_difference_terms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, choices=(2, 3, 4), default=3)
    parser.add_argument('--degree-grid', default='0,1,32,33,64,65,96,97')
    parser.add_argument('--common-offset', type=int, default=0)
    parser.add_argument('--top-degree', type=int, default=0)
    parser.add_argument('--degree-family', choices=('occupation', 'individual'), default='occupation')
    parser.add_argument('--priority-limit', type=int)
    parser.add_argument('--degree-limit', type=int)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    modes, m = args.modes, args.n * args.ramification
    states = tuple(product(range(2), repeat=modes))
    size = len(states)
    lookup = {state: index for index, state in enumerate(states)}
    grid = tuple(map(int, args.degree_grid.split(',')))
    if args.degree_family == 'occupation':
        degrees = tuple(tuple(args.common_offset + (values + (args.top_degree,))[sum(state)]
                              for state in states) for values in product(grid, repeat=modes))
    else:
        degrees = tuple(tuple(args.common_offset + value for value in values + (args.top_degree,))
                        for values in product(grid, repeat=size - 1))
    if args.degree_limit:
        degrees = degrees[:args.degree_limit]
    natural = tuple(range(size))
    hamming = tuple(sorted(range(size), key=lambda index: (index.bit_count(), index)))
    base_orders = (natural, natural[::-1], hamming, hamming[::-1])
    priorities = tuple(dict.fromkeys(order[shift:] + order[:shift]
                                    for order in base_orders for shift in range(size)))
    if args.priority_limit:
        priorities = priorities[:args.priority_limit]
    transform = finite_difference_terms(2, modes)
    catalogue = {}
    feasible_targets = set()
    calls = 0
    started = monotonic()
    for index, caps in enumerate(degrees):
        maxima = tuple(max(caps[lookup[state]] + m * power for state, power, scalar in terms) + 1
                       for terms in transform.values())
        for priority in priorities:
            target = [0] * size
            for row in priority:
                low, high = target[row], maxima[row]
                while high - low > 1:
                    target[row] = middle = (low + high) // 2
                    calls += 1
                    if synthesize(modes, m, caps, target, construct=False):
                        low = middle
                    else:
                        high = middle
                target[row] = low
            if (caps, tuple(target)) in feasible_targets:
                continue
            feasible_targets.add((caps, tuple(target)))
            witness = synthesize(modes, m, caps, target)
            assert witness and witness['actual_orders'] == tuple(target)
            key = witness['actual_orders'], witness['actual_degrees']
            if key not in catalogue:
                alpha = [{'state': state, 'polynomial': sorted(polynomial.items())}
                         for state, polynomial in witness['alpha'].items()]
                catalogue[key] = {
                    'beta_orders': witness['actual_orders'], 'degrees': witness['actual_degrees'],
                    'source_degree_bounds': caps, 'source_beta_bounds': tuple(target),
                    'endpoint_values': witness['endpoint_values'],
                    'endpoint_image_dimension': witness['endpoint_image_dimension'],
                    'coefficient_field_order': witness['coefficient_field_order'],
                    'witness_polynomial_sha256': sha256(json.dumps(alpha, separators=(',', ':')).encode()).hexdigest(),
                    'multiplicity': 0,
                }
            catalogue[key]['multiplicity'] += 1
        if (index + 1) % 20 == 0:
            print('PROGRESS', index + 1, len(degrees), calls, len(feasible_targets),
                  len(catalogue), local_endpoint_image.cache_info(), monotonic() - started, flush=True)
    rows = tuple(catalogue.values())
    result = {'status': 'Exact residue-compressed polynomial witnesses; no pressure claim without a selected-index replay.',
              'n': args.n, 'ramification': args.ramification, 'modes': modes, 'm': m,
              'coefficient_field_order': 2 * size, 'degree_grid': grid,
              'degree_family': args.degree_family,
              'common_offset': args.common_offset, 'top_degree': args.top_degree,
              'degree_vectors': len(degrees), 'row_priority_orders': priorities,
              'exact_feasibility_queries': calls, 'distinct_feasible_bound_pairs': len(feasible_targets),
              'distinct_ordered_beta_degree_profiles': len(rows), 'catalogue': rows,
              'catalogue_sha256': sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest(),
              'elapsed_seconds': monotonic() - started}
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print('DONE', calls, len(feasible_targets), len(rows), result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
