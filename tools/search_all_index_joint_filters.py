#!/usr/bin/env python3
"""Reuse exact joint beta profiles over every selected subset at prime-power n."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import log
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.special import logsumexp

from joint_identity_filter_pressure import compiled_binary_beta_orders
from prime_state_filter_pressure import gaussian_pressure
from search_f4_repaired_identity_filters import candidate_exponents


def monomial_catalogue(n, P, modes, family):
    size = 2 ** modes
    if family == 'three-mode-f4':
        if modes != 3:
            raise ValueError('The complete colored family uses three modes.')
        exponents = candidate_exponents(n, P, modes)
        colors = tuple((1,) + tail for tail in product((1, 2, 3), repeat=size - 1))
    elif family == 'binary-indicators':
        exponents = ((0,) + tuple(sign * ((mask >> j) & 1) for j in range(size - 1))
                     for sign in (-1, 1) for mask in range(1 << (size - 1)))
        colors = ((1,) * size,)
    else:
        raise ValueError('Unknown complete coefficient family.')
    catalogue = {}
    cases = 0
    for index, exponents_row in enumerate(exponents):
        evaluate = compiled_binary_beta_orders(n, P, modes, exponents_row)
        degrees = tuple(-value for value in exponents_row)
        for coefficients in colors:
            orders = evaluate(coefficients)
            key = orders, degrees
            if key not in catalogue:
                catalogue[key] = {'beta_orders': orders, 'degrees': degrees,
                                  'exponents': exponents_row, 'colors': coefficients, 'multiplicity': 0}
            catalogue[key]['multiplicity'] += 1
            cases += 1
    return tuple(catalogue.values()), cases


def lattice_catalogue(path):
    """Transport exact row certificates; a strict gain is replayed literally."""
    from weyl_lattice_filter_synthesis import lattice_basis
    source = json.loads(Path(path).read_text())
    if source.get('catalogue_format') == 'joint-lattice-screen-rle-v1':
        from compact_joint_lattice_screen import expand_rows
        source_rows = expand_rows(source)
    else:
        source_rows = source['rows']
    n, P, modes, geometry = (source[key] for key in ('n', 'ramification', 'modes', 'geometry'))
    data = lattice_basis(n, P, geometry, modes)
    states = data['states']
    lookup = {state: index for index, state in enumerate(states)}
    rows = []
    for index in range(len(states)):
        support = [state for state, coefficients in zip(states, data['coefficients'])
                   if any(row == index for row, column in coefficients)]
        minimum = tuple(min(state[j] for state in support) for j in range(modes))
        k = tuple(1 - value for value in minimum)
        weight = P * sum(j * value for j, value in zip(data['selected'], minimum))
        rows.append((lookup[k], weight))
    assert len({index for index, weight in rows}) == len(states)
    for block in data['blocks']:
        for coefficients in data['coefficients']:
            for (row, column), value in coefficients.items():
                if block['row_offset'] <= row < block['row_offset'] + len(block['indices']):
                    assert row - block['row_offset'] == column
    catalogue = {}
    for source_index, row in enumerate(source_rows):
        beta = [0] * len(states)
        for target, (index, weight) in zip(row['row_targets'], rows):
            beta[index] = target - weight
        beta, degrees = tuple(beta), tuple(row['actual_degrees'])
        key = beta, degrees
        if key not in catalogue:
            catalogue[key] = {'beta_orders': beta, 'degrees': degrees,
                              'source_row': source_index, 'source_degree_bounds': row['degree_bounds'],
                              'source_row_targets': row['row_targets'],
                              'witness_polynomial_sha256': row['witness_polynomial_sha256'],
                              'multiplicity': 0}
        catalogue[key]['multiplicity'] += 1
    return tuple(catalogue.values()), len(source_rows), source


def permutation_closure(catalogue, modes):
    states = tuple(product(range(2), repeat=modes))
    lookup = {state: index for index, state in enumerate(states)}
    maps = tuple((permutation, tuple(lookup[tuple(state[j] for j in permutation)] for state in states))
                 for permutation in permutations(range(modes)))
    result = {}
    for row in catalogue:
        for permutation, indices in maps:
            beta = tuple(row['beta_orders'][i] for i in indices)
            degrees = tuple(row['degrees'][i] for i in indices)
            key = beta, degrees
            if key in result:
                continue
            value = dict(row)
            value.update({'beta_orders': beta, 'degrees': degrees, 'source_mode_permutation': permutation})
            for field in ('exponents', 'colors', 'endpoint_values'):
                if field in value:
                    value[field] = tuple(value[field][i] for i in indices)
            if 'multiplicity' in value:
                value['source_multiplicity'] = value.pop('multiplicity')
            result[key] = value
    return tuple(result.values())


def subset_weights(n, modes, subsets, scale, charge):
    states = np.asarray(tuple(product(range(2), repeat=modes)), dtype=float)
    slack = 2 * n - 1 - 2 * np.asarray(subsets, dtype=float)
    zero = -scale * (slack + charge) @ states.T
    infinity = -scale * slack @ states.T
    return np.exp(zero - logsumexp(zero, axis=1)[:, None]), np.exp(
        infinity - logsumexp(infinity, axis=1)[:, None])


def search(n, P, modes, catalogue, subsets, scales, charges):
    states = np.asarray(tuple(product(range(2), repeat=modes)), dtype=float)
    occupations = states.sum(axis=1)
    orders = np.asarray([row['beta_orders'] for row in catalogue], dtype=float)
    degrees = np.asarray([row['degrees'][::-1] for row in catalogue], dtype=float)
    best_relative = np.full(len(catalogue), np.inf)
    best_critical = np.full(len(catalogue), np.inf)
    relative_point = [None] * len(catalogue)
    critical_point = [None] * len(catalogue)
    summaries = []
    started = monotonic()
    for scale, charge in product(scales, charges):
        wz, wi = subset_weights(n, modes, subsets, scale, charge)
        az = 2 * scale * (n * occupations[None, :] - orders / P)
        ai = 2 * scale * degrees / P
        shift_z, shift_i = az.max(axis=1), ai.max(axis=1)
        az, ai = np.exp(az - shift_z[:, None]), np.exp(ai - shift_i[:, None])
        gaussian = gaussian_pressure(2, n, charge, scale)
        grid_best = np.inf
        for start in range(0, len(catalogue), 256):
            stop = min(start + 256, len(catalogue))
            values = np.log(wz @ az[start:stop].T) + np.log(wi @ ai[start:stop].T)
            values += (shift_z[start:stop] + shift_i[start:stop])[None, :]
            indices = np.argmin(values, axis=0)
            minima = values[indices, np.arange(stop - start)]
            grid_best = min(grid_best, float(minima.min()))
            for local, (index, value) in enumerate(zip(indices, minima)):
                at = start + local
                point = (float(scale), float(charge), tuple(subsets[index]))
                if value < best_relative[at]:
                    best_relative[at], relative_point[at] = value, point
                if value + gaussian < best_critical[at]:
                    best_critical[at], critical_point[at] = value + gaussian, point
        summaries.append({'scale': scale, 'charge': charge, 'Gaussian_pressure': gaussian,
                          'minimum_relative_log_cost': grid_best,
                          'minimum_critical_pressure': grid_best + gaussian})
        print('GRID', scale, charge, grid_best, grid_best + gaussian, monotonic() - started, flush=True)
    selected = set(np.argsort(best_relative)[:30]) | set(np.argsort(best_critical)[:30])
    selected.update(index for index, value in enumerate(best_relative) if value < -1e-9)
    records = []
    for index in sorted(selected):
        records.append({**catalogue[index], 'catalogue_index': int(index),
                        'minimum_relative_log_cost': float(best_relative[index]),
                        'relative_point': relative_point[index],
                        'minimum_critical_pressure': float(best_critical[index]),
                        'critical_point': critical_point[index]})
    return {'minimum_relative_log_cost': float(best_relative.min()),
            'minimum_critical_pressure': float(best_critical.min()),
            'same_gauge_gain_found': bool((best_relative < -1e-9).any()),
            'negative_critical_pressure_found': bool((best_critical < -1e-9).any()),
            'grid_summaries': summaries, 'selected_profile_records': records,
            'elapsed_seconds': monotonic() - started}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--modes', type=int, choices=(3, 4), default=3)
    parser.add_argument('--family', choices=('three-mode-f4', 'binary-indicators'), default='three-mode-f4')
    parser.add_argument('--lattice-source', type=Path,
                        help='Reuse exact polynomial row certificates from this complete lattice record.')
    parser.add_argument('--joint-source', type=Path,
                        help='Reuse an exact residue-compressed joint polynomial catalogue.')
    parser.add_argument('--scales', help='Comma-separated t=P*log(s) values.')
    parser.add_argument('--charges', default='-1,0,1')
    parser.add_argument('--subset-stride', type=int, default=1)
    parser.add_argument('--permute-modes', action='store_true',
                        help='Close a polynomial catalogue under all assignments of its modes to the selected subset.')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = monotonic()
    if args.joint_source:
        source = json.loads(args.joint_source.read_text())
        catalogue = tuple(source['catalogue'])
        cases = source['distinct_feasible_bound_pairs']
        args.n, args.ramification, args.modes = source['n'], source['ramification'], source['modes']
    elif args.lattice_source:
        catalogue, cases, source = lattice_catalogue(args.lattice_source)
        args.n, args.ramification, args.modes = source['n'], source['ramification'], source['modes']
    else:
        catalogue, cases = monomial_catalogue(args.n, args.ramification, args.modes, args.family)
    before_permutations = len(catalogue)
    if args.permute_modes:
        catalogue = permutation_closure(catalogue, args.modes)
    payload = json.dumps(catalogue, separators=(',', ':'))
    print('CATALOGUE', cases, len(catalogue), monotonic() - started, flush=True)
    subsets = tuple(combinations(range(args.n), args.modes))[::args.subset_stride]
    scales = tuple(map(float, args.scales.split(','))) if args.scales else tuple(
        value / args.n for value in (2, 4, 5, 6, 8))
    charges = tuple(map(float, args.charges.split(',')))
    result = {'status': 'Exact coefficient profiles and complete stated subset/grid screen; floating pressure is exploratory.',
              'n': args.n, 'ramification': args.ramification, 'modes': args.modes,
              'family': ('exact residue-compressed polynomial witnesses' if args.joint_source else
                         ('transported polynomial row certificates' if args.lattice_source else args.family)),
              'source_lattice_record': str(args.lattice_source) if args.lattice_source else None,
              'source_joint_record': str(args.joint_source) if args.joint_source else None,
              'coefficient_cases': cases,
              'all_mode_permutations_included': args.permute_modes,
              'ordered_profiles_before_mode_permutations': before_permutations,
              'distinct_ordered_beta_degree_profiles': len(catalogue),
              'catalogue_sha256': sha256(payload.encode()).hexdigest(),
              'selected_subsets': len(subsets), 'subset_stride': args.subset_stride,
              'scale_grid': scales, 'charge_grid': charges, 'gain_detection_tolerance': 1e-9,
              'total_profile_subset_grid_evaluations': len(catalogue) * len(subsets) * len(scales) * len(charges)}
    result.update(search(args.n, args.ramification, args.modes, catalogue, subsets, scales, charges))
    result['total_elapsed_seconds'] = monotonic() - started
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('DONE', result['minimum_relative_log_cost'], result['minimum_critical_pressure'], flush=True)


if __name__ == '__main__':
    main()
