#!/usr/bin/env python3
"""Reconstruct retained coefficient-normal-form tensors and certify profiles."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_block_smith_crossings import binary_rank
from check_f4_recentered_walks import certify_profile, literal_matrix, same_cost
from search_f4_recentered_crossings import digest, five_letter_source, root_for_profiles


def replay(path):
    result = json.loads(path.read_text())
    source_path = Path(result['arguments']['directions'])
    assert sha256(source_path.read_bytes()).hexdigest() == result['source_record_sha256']
    source = json.loads(source_path.read_text())
    descriptors = result['slab_direction_descriptors']
    directions = [[value << shift for value in source['primitives'][i]['normalized_slab']]
                  for i, shift in descriptors]
    assert binary_rank([sum(value << (19 * j) for j, value in enumerate(row)) for row in directions]) == len(directions)
    assert len(directions) == result['independent_slab_directions']
    base = five_letter_source()[0]
    assert digest(base) == result['reference_digest']
    profiles_checked = endpoint_entries = witnesses = 0
    gauges = {}
    for name, gauge in result['gauges'].items():
        energy = gauge['energies']
        half, zero = tuple(value // 2 for value in energy), (0,) * 5
        heights = (zero, zero, half, half, half, half, zero, zero)
        for side in range(2):
            matrix, scalar = literal_matrix(base, 5, 18, heights, side)
            certify_profile(matrix, scalar, gauge['baseline_profiles'][side])
            profiles_checked += 1
        same_cost(root_for_profiles(gauge['baseline_profiles'], energy, energy), gauge['baseline_cost'])
        costs = [gauge['baseline_cost']['rho']]
        for witness in gauge['retained']:
            values = [list(pair) for pair in base]
            assert len(witness['slab_parameter_masks']) == 5
            for x, parameters in enumerate(witness['slab_parameter_masks']):
                assert parameters >= 0 and parameters.bit_length() <= len(directions)
                while parameters:
                    bit = parameters & -parameters
                    parameters ^= bit
                    direction = directions[bit.bit_length() - 1]
                    for j, value in enumerate(direction):
                        values[x + 5 * j][0] ^= value
            assert values == witness['values']
            assert digest(values) == witness['digest']
            for index, (low, high) in enumerate(values):
                x, y, u, v = ((index // 5 ** leg) % 5 for leg in range(4))
                assert low.bit_length() <= 19 and high == 0
                assert low.bit_count() % 2 == int(x == v and y == u)
                endpoint_entries += 1
            for side in range(2):
                matrix, scalar = literal_matrix(values, 5, 18, heights, side)
                certify_profile(matrix, scalar, witness['profiles'][side])
                profiles_checked += 1
            cost = root_for_profiles(witness['profiles'], energy, energy)
            same_cost(cost, witness)
            costs.append(cost['rho'])
            witnesses += 1
        assert abs(min(costs) - gauge['best_rho']) < 2e-8
        gauges[name] = {'retained_witnesses': len(gauge['retained']),
                        'best_rho': gauge['best_rho'],
                        'best_nonzero_witness_rho': min(w['rho'] for w in gauge['retained'])}
    assert result['completed_trials'] <= result['arguments']['trials']
    assert result['profile_queries'] + result['duplicate_tensors'] == result['completed_trials']
    return {'status': 'literal coefficient normal forms and complete polynomial profiles certified',
            'record_sha256': sha256(path.read_bytes()).hexdigest(),
            'source_record_sha256': result['source_record_sha256'],
            'independent_slab_directions': len(directions),
            'completed_trials': result['completed_trials'],
            'distinct_nonreference_profile_queries': result['profile_queries'],
            'retained_witnesses': witnesses, 'literal_endpoint_entries': endpoint_entries,
            'complete_profiles_certified_by_polynomial_frames': profiles_checked,
            'gauges': gauges, 'factoring_bound_improved': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path)
    args = parser.parse_args()
    print(json.dumps(replay(args.record), indent=2))
