#!/usr/bin/env python3
"""Replay exact module membership, row targets, and retained complete profiles."""

from hashlib import sha256
import json
from pathlib import Path

from check_f4_recentered_walks import certify_profile, literal_matrix, same_cost
from f4_recentered_lattices import multiply_matrix
from search_f4_recentered_crossings import digest, increment, root_for_profiles
from source_module_crossing_fiber import SourceModuleCrossingFiber


def replay(path):
    record = json.loads(path.read_text())
    source_path = Path(record['arguments']['directions'])
    assert sha256(source_path.read_bytes()).hexdigest() == record['directions_record_sha256']
    source = json.loads(source_path.read_text())
    degree = record['arguments']['direction_degree']
    source['direction_descriptors'] = [descriptor for descriptor in source['direction_descriptors']
        if source['primitives'][descriptor[0]]['normalized_degree'] + descriptor[2] <= degree]
    baseline, heights = record['reference_values'], record['heights_infinity_then_zero']
    assert digest(baseline) == record['reference_digest']
    fiber = SourceModuleCrossingFiber(baseline, source, heights, record['arguments']['raise_cap'])
    assert tuple(record['reference_targets']) == fiber.orders
    assert tuple(record['ceiling_targets']) == fiber.caps
    assert fiber.variables == record['direction_count'] == record['independent_direction_count']
    pivots = fiber.solve(fiber.orders)
    assert pivots is not None
    assert fiber.variables - len(pivots) == record['baseline_preserving_dimension'] == 0
    raises = [row for row in range(50) if increment(fiber, pivots, row, fiber.orders[row]) is not None]
    assert raises == record['complete_single_row_raises'] == []
    forced_controls = []
    for row, reported in enumerate(record['relaxed_forced_raise_feasibility'][:50]):
        assert reported['rows'] == [row]
        targets = [0] * 50
        targets[row] = fiber.orders[row] + 1
        feasible = fiber.solve(targets) is not None
        assert feasible == reported['feasible']
        forced_controls.append(feasible)
    ex, ey = record['paired_energies']
    profiles_checked = entries = targets_checked = 0
    for side in range(2):
        literal, scalar = literal_matrix(baseline, 5, 18, heights, side)
        certify_profile(literal, scalar, record['baseline_profiles'][side])
        profiles_checked += 1
    same_cost(root_for_profiles(record['baseline_profiles'], ex, ey), record['baseline_cost'])
    for witness in record['retained_witnesses']:
        values = witness['values']
        assert digest(values) == witness['digest'] and fiber.contains(values)
        for index, (low, high) in enumerate(values):
            x, y, u, v = ((index // 5 ** leg) % 5 for leg in range(4))
            assert low.bit_length() <= 19 and high.bit_length() <= 19
            assert low.bit_count() % 2 == int(x == v and y == u)
            assert high.bit_count() % 2 == 0
            entries += 1
        for side in range(2):
            literal, scalar = literal_matrix(values, 5, 18, heights, side)
            certify_profile(literal, scalar, witness['profiles'][side])
            profiles_checked += 1
            left, right = fiber.full_lattice.frames[side]
            precision = fiber.full_lattice.precisions[side]
            actual = multiply_matrix(multiply_matrix(left, literal, precision), right, precision)
            for row in range(25):
                target = witness['targets'][25 * side + row]
                mask = (1 << (2 * target)) - 1
                assert all(not value & mask for value in actual[row])
                targets_checked += 1
        same_cost(root_for_profiles(witness['profiles'], ex, ey), witness)
    assert record['completed_requested_schedule']
    assert record['completed_trials'] == record['arguments']['trials']
    assert record['distinct_target_vectors'] == record['complete_profile_queries']
    assert sum(row['samples'] for row in record['profile_classes']) == record['complete_profile_queries']
    assert abs(min([record['baseline_cost']['rho']] + [row['rho'] for row in record['retained_witnesses']])
               - record['best_recorded_rho']) < 1e-10
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(),
            'source_record_sha256': record['directions_record_sha256'],
            'status': 'passed', 'complete_module_dimension': fiber.variables,
            'baseline_preserving_dimension': 0, 'single_row_raises_feasible': raises,
            'fully_relaxed_single_raise_feasible_rows': [row for row, value in enumerate(forced_controls) if value],
            'completed_target_schedules': record['completed_trials'],
            'recorded_complete_profile_queries': record['complete_profile_queries'],
            'retained_literal_tensors': len(record['retained_witnesses']),
            'literal_endpoint_entries': entries, 'literal_frame_row_targets_checked': targets_checked,
            'complete_profiles_certified_by_polynomial_frames': profiles_checked,
            'best_recorded_rho': record['best_recorded_rho']}


if __name__ == '__main__':
    full_path = Path('evidence/aligned-gaussian-module-directions-q5.json')
    full = json.loads(full_path.read_text())
    low = json.loads(Path('evidence/aligned-gaussian-module-low-degree4-directions.json').read_text())
    assert low['original_module_record_sha256'] == sha256(full_path.read_bytes()).hexdigest()
    assert low['primitives'] == full['primitives']
    assert low['direction_descriptors'] == [item for item in full['direction_descriptors']
        if full['primitives'][item[0]]['normalized_degree'] + item[2] <= 4]
    assert len(low['direction_descriptors']) == low['direction_count_before_rank_reduction'] == 60
    paths = [Path('evidence/aligned-gaussian-module-lattice-' + name + '-d18.json')
             for name in ('natural', 'two-zero')]
    print(json.dumps({'status': 'exact complete-module finite-search controls passed',
                      'complete_degree_four_slice_directions': 60,
                      'records': [replay(path) for path in paths],
                      'factoring_bound_improved': False}, indent=2))
