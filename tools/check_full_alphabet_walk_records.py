#!/usr/bin/env python3
"""Replay complete saved tensors from finite target walks, including timed runs."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from finite_precision_smith_frames import build_lattice
from search_full_alphabet_two_cut_lattices import root_for_profiles
from check_full_alphabet_two_cut_lattices import independent_profiles


def replay(path, independent_reference=True):
    data = json.loads(path.read_text())
    args = data['arguments']
    assert 0 < data['completed_trials'] <= args['trials']
    assert data['completed_requested_trials'] == (data['completed_trials'] == args['trials'])
    for name, digest in data['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == digest, name
    reference_path = Path(args['reference'])
    if not reference_path.is_absolute():
        reference_path = Path(__file__).resolve().parents[1]/reference_path
    assert sha256(reference_path.read_bytes()).hexdigest() == data['reference_sha256']
    assert json.loads(reference_path.read_text()) == data['reference']
    assert len(data['classes']) == data['distinct_profile_classes']
    assert sum(row['samples'] for row in data['classes']) == data['complete_profile_queries']
    assert data['complete_profile_queries'] <= args['samples']*data['distinct_target_vectors']
    model = build_lattice(data['binary_reference_values'], data['physical_alphabet'],
        data['degree'], data['heights_infinity_then_zero'], raise_cap=args['upper_window'])
    assert model.variables == data['variables']
    assert model.orders == tuple(data['reference_row_targets'])
    assert model.caps == tuple(data['row_target_ceilings'])
    reference_profiles = model.profiles(data['binary_reference_values'])
    assert reference_profiles == tuple(map(tuple, data['reference_profiles']))
    toeplitz = 0
    if independent_reference:
        alternate, toeplitz = independent_profiles(model, data['binary_reference_values'])
        assert alternate == reference_profiles
    reference_cost = root_for_profiles(reference_profiles, *data['physical_paired_energies'])
    assert abs(reference_cost['rho']-data['reference_cost']['rho']) < 1e-9
    profiles_in_catalogue = {tuple(map(tuple, row['profiles'])) for row in data['classes']}
    costs = [reference_cost['rho']]
    negative = False
    for witness in data['retained_witnesses']:
        assert all(floor <= target <= cap for floor, target, cap in zip(
            data['row_target_floors'], witness['targets'], data['row_target_ceilings']))
        profiles = model.check_witness(witness['values'], witness['targets'])
        assert profiles == tuple(map(tuple, witness['profiles']))
        assert profiles in profiles_in_catalogue
        cost = root_for_profiles(profiles, *data['physical_paired_energies'])
        assert abs(cost['rho']-witness['rho']) < 1e-9
        assert abs(cost['critical_pressure']-witness['critical_pressure']) < 1e-8
        costs.append(cost['rho'])
        negative |= cost['critical_pressure'] < -1e-8
    assert abs(min(costs)-data['best_recorded_rho']) < 1e-9
    assert negative == data['strict_crossing_candidate_found']
    return {'record': path.name, 'physical_alphabet': model.q, 'degree': model.degree,
        'completed_trials': data['completed_trials'], 'requested_trials': args['trials'],
        'completed_requested_trials': data['completed_requested_trials'],
        'recorded_time_budget_exhausted': data['elapsed_seconds'] >= args['seconds'],
        'distinct_target_vectors': data['distinct_target_vectors'],
        'distinct_profile_classes': data['distinct_profile_classes'],
        'complete_profile_queries': data['complete_profile_queries'],
        'literal_retained_tensor_replays': len(data['retained_witnesses']),
        'independent_reference_toeplitz_ranks': toeplitz,
        'best_replayed_rho': data['best_recorded_rho'],
        'strict_crossing_candidate_found': negative}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--skip-independent-reference', action='store_true')
    args = parser.parse_args()
    result = {'status': 'Every retained tensor and full cost replayed; incomplete trial budgets remain explicit.',
        'records': [replay(path, not args.skip_independent_reference) for path in args.records]}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
