#!/usr/bin/env python3
"""Regenerate the affine family and independently replay original polynomials."""

import argparse
import json
from pathlib import Path
import random

from affine_symmetric_binary_filters import affine_relabelings, LiteralRelabeler
from check_symmetric_binary_filter_synthesis import full_transform, xor
from check_symmetric_filter_catalogues import physical_entropy
from prime_state_filter_pressure import gaussian_pressure
from search_affine_symmetric_tiles import AffineTiledPressure, build_orbit, literal_profile
from symmetric_binary_filter_synthesis import synthesize


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=root / 'evidence' / 'symmetric-tiled-filters-r4.json')
    parser.add_argument('--record', type=Path, default=root / 'evidence' / 'affine-symmetric-tiled-filters.json')
    parser.add_argument('--output', type=Path, default=root / 'evidence' / 'affine-symmetric-tiled-filter-check.json')
    args = parser.parse_args()
    source, record = json.loads(args.source.read_text()), json.loads(args.record.read_text())
    catalogue, manifest = build_orbit(source)
    assert json.loads(json.dumps(manifest)) == record['orbit']
    for stored in record['literal_representatives'] + record['selected_literal_profiles']:
        literal = literal_profile(source, catalogue, manifest, stored['catalogue_index'])
        assert json.loads(json.dumps(literal)) == stored
    generator = random.Random(720497)
    exact_states = 0
    for _ in range(64):
        source_row = generator.choice(source['catalogue'])
        mapping = generator.choice(manifest['affine_maps'])
        witness = synthesize(4, source_row['source_m'], source_row['source_degree_bounds'], source_row['source_order_bounds'])
        full = tuple(witness['f'][weight] for weight in mapping['weights'])
        beta = full_transform(full, 4, source_row['source_m'])
        common = min(min(value) for value in beta)
        normalized_f = tuple({degree - common: coefficient for degree, coefficient in value.items()} for value in full)
        normalized_beta = tuple({degree - common: coefficient for degree, coefficient in value.items()} for value in beta)
        compiled = LiteralRelabeler(witness).evaluate(mapping['weights'], literal=True)
        assert compiled['f'] == normalized_f and compiled['beta'] == normalized_beta
        assert tuple(map(min, beta)) == tuple(value + common for value in compiled['beta_orders'])
        assert tuple(xor(value.values()) for value in normalized_f) == compiled['endpoint_values']
        exact_states += 16

    # An old beta profile cannot be permuted with the occupation labels.
    witness = {'modes': 4, 'm': 4, 'coefficient_field_order': 2,
               'f': ({0: 1, 8: 1, 16: 1}, {0: 1}, {0: 1}, {0: 1}, {0: 1}),
               'actual_degrees': (16, 0, 0, 0, 0), 'endpoint_values': (1, 1, 1, 1, 1)}
    mapping = affine_relabelings(4)[1]
    actual = LiteralRelabeler(witness).evaluate(mapping['weights'], literal=True)
    old_orders = (0, 4, 24, 8, 8)
    incorrectly_permuted = tuple(old_orders[weight] for weight in mapping['weights'])
    assert actual['beta_orders'] != incorrectly_permuted
    assert actual['beta_orders'][0] == 0 and incorrectly_permuted[0] == 4

    roots = []
    for stored in record['optimized_runs']:
        blocks, n = stored['actual_mode_energy_blocks'], stored['n']
        assert sorted(value for block in blocks for value in block) == list(range(1, 2 * n, 2))
        evaluator = AffineTiledPressure(n, catalogue, blocks)
        relative, selected, pieces = evaluator.relative_cost(stored['tau'], stored['scaled_charge'])
        assert abs(relative - stored['relative_log_cost']) < 1e-8
        assert selected.tolist() == stored['selected_catalogue_indices']
        residual = gaussian_pressure(2, n, stored['charge'], stored['scale'], stored['rho']) + relative
        critical = gaussian_pressure(2, n, stored['charge'], stored['scale']) + relative
        assert abs(residual) < 1e-7 and abs(critical - stored['full_critical_pressure']) < 1e-8
        entropy = physical_entropy(2, n, stored['charge'], stored['scale'], stored['rho'])
        assert entropy > 1e-8
        roots.append({'n': n, 'partition': stored['partition'], 'rho': stored['rho'],
                      'root_residual': residual, 'physical_entropy': entropy,
                      'nonbaseline_blocks': int((selected >= 0).sum())})
        evaluator.exponential_arrays.cache_clear()
    result = {'status': 'Complete affine family regenerated, literal original-state controls and full roots replayed.',
              'scope': 'The stated finite affine family only; no all-affine or all-filter exclusion.',
              'coefficient_cases': manifest['complete_coefficient_cases'],
              'retained_profiles': len(catalogue), 'baseline_dominated_cases': manifest['baseline_dominated_cases'],
              'independent_original_polynomial_state_controls': exact_states,
              'old_beta_permutation_negative_control': {'actual_state_zero_order': 0, 'incorrect_permuted_order': 4},
              'complete_root_replays': roots}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
