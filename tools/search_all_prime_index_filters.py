#!/usr/bin/env python3
"""Screen exact joint p-state polynomial profiles on every selected subset."""
import argparse
from itertools import combinations, permutations, product
import json
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.special import logsumexp

from joint_identity_filter_pressure import relative_log_cost
from prime_state_filter_pressure import gaussian_pressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--charges', default='-27,-9,-1,0,1,9,27')
    parser.add_argument('--scales')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.source.read_text())
    p, n, P, r = (source[key] for key in ('prime', 'n', 'ramification', 'modes'))
    states = tuple(product(range(p), repeat=r)); lookup = {state: i for i, state in enumerate(states)}
    catalogue = {}
    for row in source['catalogue']:
        for permutation in permutations(range(r)):
            indices = tuple(lookup[tuple(state[j] for j in permutation)] for state in states)
            beta = tuple(row['beta_orders'][i] for i in indices)
            degrees = tuple(row['degrees'][i] for i in indices)
            if (beta, degrees) not in catalogue:
                catalogue[beta, degrees] = {**row, 'beta_orders': beta, 'degrees': degrees,
                                             'source_mode_permutation': permutation}
    rows = tuple(catalogue.values())
    subsets = tuple(combinations(range(n), r))
    scales = tuple(map(float, args.scales.split(','))) if args.scales else tuple(v / n for v in (2, 4, 5, 6, 8))
    charges = tuple(map(float, args.charges.split(',')))
    digits = np.asarray(states, dtype=float)
    occupations = digits.sum(axis=1)
    beta = np.asarray([row['beta_orders'] for row in rows], dtype=float)
    degrees = np.asarray([row['degrees'][::-1] for row in rows], dtype=float)
    slack = 2 * n - 1 - 2 * np.asarray(subsets, dtype=float)
    best_relative = np.full(len(rows), np.inf); best_critical = np.full(len(rows), np.inf)
    rel_points = [None] * len(rows); crit_points = [None] * len(rows)
    summaries = []; started = monotonic()
    for scale, charge in product(scales, charges):
        lz = -scale * (slack + charge) @ digits.T; li = -scale * slack @ digits.T
        wz = np.exp(lz - logsumexp(lz, axis=1)[:, None]); wi = np.exp(li - logsumexp(li, axis=1)[:, None])
        az = 2 * scale * (n * occupations[None, :] - beta / P); ai = 2 * scale * degrees / P
        sz, si = az.max(axis=1), ai.max(axis=1)
        az, ai = np.exp(az - sz[:, None]), np.exp(ai - si[:, None])
        gaussian = gaussian_pressure(p, n, charge, scale)
        grid_best = np.inf
        for start in range(0, len(rows), 256):
            stop = min(start + 256, len(rows))
            values = np.log(wz @ az[start:stop].T) + np.log(wi @ ai[start:stop].T)
            values += (sz[start:stop] + si[start:stop])[None, :]
            indices = values.argmin(axis=0); minima = values[indices, np.arange(stop - start)]
            grid_best = min(grid_best, float(minima.min()))
            for at, (index, value) in enumerate(zip(indices, minima), start):
                point = (scale, charge, subsets[index])
                if value < best_relative[at]:best_relative[at], rel_points[at] = value, point
                if value + gaussian < best_critical[at]:best_critical[at], crit_points[at] = value + gaussian, point
            # Independent stabilized pressure evaluation at one genuine
            # selected subset per block protects the batched normalization.
            direct = relative_log_cost(p, n, P, charge, scale, subsets[0],
                                       rows[start]['beta_orders'], rows[start]['degrees'])
            assert abs(direct - values[0, 0]) < 1e-9
        summaries.append({'scale': scale, 'charge': charge, 'Gaussian_pressure': gaussian,
                          'minimum_relative_log_cost': grid_best,
                          'minimum_critical_pressure': grid_best + gaussian})
        print('GRID', scale, charge, grid_best, grid_best + gaussian, monotonic() - started, flush=True)
    selected = set(np.argsort(best_relative)[:30]) | set(np.argsort(best_critical)[:30])
    selected.update(index for index, value in enumerate(best_relative) if value < -1e-9)
    witnesses = [{**rows[i], 'minimum_relative_log_cost': float(best_relative[i]), 'relative_point': rel_points[i],
                  'minimum_critical_pressure': float(best_critical[i]), 'critical_point': crit_points[i]}
                 for i in sorted(selected)]
    result = {'status': 'Exact joint-prime coefficient profiles with a complete stated subset/grid screen; floating pressure is exploratory.',
              'source': str(args.source), 'prime': p, 'n': n, 'ramification': P, 'modes': r,
              'profiles_before_mode_permutations': len(source['catalogue']),
              'distinct_ordered_profiles': len(rows), 'all_mode_permutations_included': True,
              'selected_subsets': len(subsets), 'scale_grid': scales, 'charge_grid': charges,
              'total_profile_subset_grid_evaluations': len(rows) * len(subsets) * len(scales) * len(charges),
              'minimum_relative_log_cost': float(best_relative.min()),
              'minimum_critical_pressure': float(best_critical.min()),
              'gain_detection_tolerance': 1e-9, 'same_gauge_gain_found': bool((best_relative < -1e-9).any()),
              'negative_critical_pressure_found': bool((best_critical < -1e-9).any()),
              'grid_summaries': summaries, 'selected_profile_records': witnesses,
              'elapsed_seconds': monotonic() - started}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('DONE', result['minimum_relative_log_cost'], result['minimum_critical_pressure'], flush=True)


if __name__ == '__main__':
    main()
