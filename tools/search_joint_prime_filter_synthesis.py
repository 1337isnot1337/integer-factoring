#!/usr/bin/env python3
"""Finite Pareto search for residue-compressed joint p-state filters."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from time import monotonic

from joint_prime_filter_synthesis import difference_matrix, local_endpoint_image, synthesize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prime', type=int, default=3)
    parser.add_argument('--n', type=int, default=27)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, default=2)
    parser.add_argument('--degree-grid', default='0,1,27,28')
    parser.add_argument('--common-offset', type=int, default=0)
    parser.add_argument('--priority-limit', type=int)
    parser.add_argument('--degree-limit', type=int)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    p, r, m = args.prime, args.modes, args.n * args.ramification
    states = tuple(product(range(p), repeat=r)); size = len(states)
    grid = tuple(map(int, args.degree_grid.split(',')))
    degrees = tuple(tuple(args.common_offset + (values + (0,))[sum(state)] for state in states)
                    for values in product(grid, repeat=r * (p - 1)))
    if args.degree_limit:
        degrees = degrees[:args.degree_limit]
    natural = tuple(range(size)); hamming = tuple(sorted(range(size), key=lambda i: (sum(states[i]), i)))
    priorities = tuple(dict.fromkeys(order[shift:] + order[:shift]
                                    for order in (natural, natural[::-1], hamming, hamming[::-1])
                                    for shift in range(size)))
    if args.priority_limit:
        priorities = priorities[:args.priority_limit]
    matrix = difference_matrix(p, r)
    seen, catalogue = set(), {}
    calls = 0; started = monotonic()
    for index, caps in enumerate(degrees):
        maxima = tuple(max(caps[j] + m * max(value) for j, value in enumerate(row) if value) + 1
                       for row in matrix)
        for priority in priorities:
            target = [0] * size
            for row in priority:
                low, high = target[row], maxima[row]
                while high - low > 1:
                    target[row] = middle = (low + high) // 2
                    calls += 1
                    if synthesize(p, r, m, caps, target, construct=False):
                        low = middle
                    else:
                        high = middle
                target[row] = low
            if (caps, tuple(target)) in seen:
                continue
            seen.add((caps, tuple(target)))
            witness = synthesize(p, r, m, caps, target)
            assert witness and witness['actual_orders'] == tuple(target)
            key = witness['actual_orders'], witness['actual_degrees']
            if key not in catalogue:
                alpha = [{'state': state, 'polynomial': sorted(value.items())} for state, value in witness['alpha'].items()]
                catalogue[key] = {'beta_orders': witness['actual_orders'], 'degrees': witness['actual_degrees'],
                                  'source_degree_bounds': caps, 'source_beta_bounds': tuple(target),
                                  'endpoint_values': witness['endpoint_values'],
                                  'endpoint_image_dimension': witness['endpoint_image_dimension'],
                                  'coefficient_field_order': witness['coefficient_field_order'],
                                  'witness_polynomial_sha256': sha256(json.dumps(alpha, separators=(',', ':')).encode()).hexdigest(),
                                  'multiplicity': 0}
            catalogue[key]['multiplicity'] += 1
        if (index + 1) % 10 == 0:
            print('PROGRESS', index + 1, len(degrees), calls, len(seen), len(catalogue),
                  local_endpoint_image.cache_info(), monotonic() - started, flush=True)
    rows = tuple(catalogue.values())
    result = {'status': 'Exact joint-prime polynomial synthesis; pressures require the selected-index replay.',
              'prime': p, 'n': args.n, 'ramification': args.ramification, 'modes': r, 'm': m,
              'coefficient_field_order': p ** (r + 1), 'degree_grid': grid,
              'common_offset': args.common_offset, 'degree_vectors': len(degrees),
              'row_priority_orders': priorities, 'exact_feasibility_queries': calls,
              'distinct_feasible_bound_pairs': len(seen), 'distinct_ordered_beta_degree_profiles': len(rows),
              'catalogue': rows, 'catalogue_sha256': sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest(),
              'elapsed_seconds': monotonic() - started}
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print('DONE', calls, len(seen), len(rows), result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
