#!/usr/bin/env python3
"""Exact affine Boolean orbits with complete full-state tiled pressure."""

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from math import lcm
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.special import logsumexp

from affine_symmetric_binary_filters import affine_relabelings, LiteralRelabeler
from search_symmetric_tiled_filters import TiledPressure, partitions, polynomial_hash
from symmetric_binary_filter_synthesis import synthesize


def build_orbit(source):
    assert source['block_size'] == 4 and source.get('prime', 2) == 2
    source_rows = source['catalogue']
    assert sha256(json.dumps(source_rows, separators=(',', ':')).encode()).hexdigest() == source['catalogue_sha256']
    mappings = affine_relabelings(4)
    common = lcm(*(row['source_m'] for row in source_rows))
    catalogue = {}
    dominated = cases = 0
    for source_index, row in enumerate(source_rows):
        witness = synthesize(4, row['source_m'], row['source_degree_bounds'], row['source_order_bounds'])
        assert witness and polynomial_hash(witness) == row['witness_polynomial_sha256']
        compiler = LiteralRelabeler(witness)
        multiplier = common // row['source_m']
        for map_index, mapping in enumerate(mappings):
            value = compiler.evaluate(mapping['weights'])
            cases += 1
            if all(order <= row['source_m'] * state.bit_count()
                   for state, order in enumerate(value['beta_orders'])):
                dominated += 1
                continue
            key = (tuple(order * multiplier for order in value['beta_orders']),
                   tuple(degree * multiplier for degree in value['degrees']))
            if key not in catalogue:
                catalogue[key] = {
                    'source_catalogue_index': source_index, 'affine_map_index': map_index,
                    'source_m': common, 'beta_orders': key[0], 'degrees': key[1],
                    'endpoint_values': value['endpoint_values'],
                    'common_shift_in_source': value['removed_common_monomial_order'],
                    'multiplicity': 0,
                }
            catalogue[key]['multiplicity'] += 1
    rows = tuple(catalogue.values())
    payload = json.dumps(rows, separators=(',', ':'))
    assert dominated + sum(row['multiplicity'] for row in rows) == cases
    manifest = {
        'source_catalogue_sha256': source['catalogue_sha256'],
        'source_profiles': len(source_rows), 'affine_maps': mappings,
        'complete_coefficient_cases': cases, 'baseline_dominated_cases': dominated,
        'distinct_retained_full_state_profiles': len(rows), 'common_m': common,
        'regenerated_catalogue_sha256': sha256(payload.encode()).hexdigest(),
    }
    return rows, manifest


def literal_profile(source, catalogue, manifest, index):
    row = catalogue[index]
    original = source['catalogue'][row['source_catalogue_index']]
    witness = synthesize(4, original['source_m'], original['source_degree_bounds'], original['source_order_bounds'])
    mapping = manifest['affine_maps'][row['affine_map_index']]
    value = LiteralRelabeler(witness).evaluate(mapping['weights'], literal=True)
    return {'catalogue_index': index, 'source_catalogue_index': row['source_catalogue_index'],
            'affine_map_index': row['affine_map_index'], 'original_m': original['source_m'],
            **value}


class AffineTiledPressure(TiledPressure):
    def __init__(self, n, catalogue, blocks):
        self.n, self.modes, self.catalogue, self.blocks = n, 4, catalogue, blocks
        self.energies = np.asarray(blocks, dtype=float) / n
        self.states = np.asarray([[(state >> bit) & 1 for bit in range(4)] for state in range(16)])
        self.weights = self.states.sum(axis=1)
        self.orders = np.asarray([[value / row['source_m'] for value in row['beta_orders']] for row in catalogue])
        self.degrees = np.asarray([[value / row['source_m'] for value in row['degrees']] for row in catalogue])
        self.calls = 0

    @lru_cache(maxsize=2)
    def exponential_arrays(self, tau):
        zero = 2 * tau * (self.weights[None, :] - self.orders)
        infinity = 2 * tau * self.degrees
        shift_zero, shift_infinity = zero.max(axis=1), infinity.max(axis=1)
        return (np.exp(zero - shift_zero[:, None]),
                np.exp(infinity - shift_infinity[:, None]), shift_zero + shift_infinity)

    def relative_cost(self, tau, charge):
        physical = -tau * self.energies @ self.states.T
        physical -= logsumexp(physical, axis=1)[:, None]
        shifted = physical - tau * charge * self.weights
        shifted -= logsumexp(shifted, axis=1)[:, None]
        physical, shifted = np.exp(physical), np.exp(shifted)
        zero, infinity, shifts = self.exponential_arrays(float(tau))
        best = np.zeros(len(self.blocks))
        selected = np.full(len(self.blocks), -1, dtype=int)
        for start in range(0, len(self.catalogue), 4096):
            stop = min(start + 4096, len(self.catalogue))
            values = (np.log(shifted @ zero[start:stop].T)
                      + np.log(physical @ infinity[start:stop].T) + shifts[start:stop][None, :])
            indices = np.argmin(values, axis=1)
            minima = values[np.arange(len(indices)), indices]
            improved = minima < best - 1e-12
            best[improved] = minima[improved]
            selected[improved] = start + indices[improved]
        return float(best.sum()), selected, best


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=root / 'evidence' / 'symmetric-tiled-filters-r4.json')
    parser.add_argument('--lengths', default='64,128')
    parser.add_argument('--iterations', type=int, default=30)
    parser.add_argument('--seed', type=int, default=720491)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.source.read_text())
    catalogue, manifest = build_orbit(source)
    print('ORBIT', manifest['complete_coefficient_cases'], len(catalogue),
          manifest['baseline_dominated_cases'], flush=True)
    samples = sorted({0, 1, 2, len(catalogue) // 4, len(catalogue) // 2,
                      3 * len(catalogue) // 4, len(catalogue) - 1})
    result = {
        'status': 'Exact orbit complete; listed full-rho searches are exploratory.',
        'source_file': args.source.name, 'orbit': manifest,
        'literal_representatives': [literal_profile(source, catalogue, manifest, index) for index in samples],
        'seed': args.seed, 'optimizer_iterations': args.iterations,
        'optimized_runs': [], 'selected_literal_profiles': [],
        'strict_crossing_candidate_found': False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    selected_profiles = set()
    for n in map(int, args.lengths.split(',')):
        if n < 4 or n & (n - 1):
            raise ValueError('Use a power-of-two physical length at least four.')
        for name, blocks in partitions(n, 4, args.seed).items():
            started = monotonic()
            run = AffineTiledPressure(n, catalogue, blocks)
            row = run.optimize(args.seed, args.iterations)
            row.update(n=n, partition=name, elapsed_seconds=monotonic() - started)
            selected_profiles.update(index for index in row['selected_catalogue_indices'] if index >= 0)
            result['optimized_runs'].append(row)
            result['selected_literal_profiles'] = [literal_profile(source, catalogue, manifest, index)
                                                   for index in sorted(selected_profiles)]
            result['strict_crossing_candidate_found'] |= row['full_critical_pressure'] < -1e-9
            args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
            print('TILED', n, name, row['rho'], row['nonbaseline_blocks'], row['elapsed_seconds'], flush=True)
            run.exponential_arrays.cache_clear()
    result['status'] = 'Exact orbit and all requested genuine-rho searches completed; no all-family exclusion is inferred.'
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')


if __name__ == '__main__':
    main()
