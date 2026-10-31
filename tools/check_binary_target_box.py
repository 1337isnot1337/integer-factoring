#!/usr/bin/env python3
"""Independent catalogue and scalar-fiber replay for the binary target box.

This checks the recorded finite target schedule, its flag manifest, all
common-system rejections, and every point of each feasible affine line.
It deliberately does not repeat the full target-by-flag-pair search.
"""
import argparse
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import random

from check_binary_two_cut_lattices import endpoint, literal_profiles
from f4_lattice_linear import dot
from search_binary_f4_flags import FlagSearch, flags, insert
from search_binary_two_cut_lattices import TwoCutLattice


SCALE = .29277969091981376


def independent_targets(reference):
    targets = set()
    # Enumerate all eight coordinates; the search instead derived its
    # eighth coordinate from seven, so this is a separate catalogue path.
    for delta in product(range(-2, 3), repeat=8):
        if sum(delta):
            continue
        value = tuple(a + b for a, b in zip(reference, delta))
        if all(value[4 * side:4 * side + 4] == tuple(sorted(value[4 * side:4 * side + 4]))
               for side in range(2)):
            targets.add(value)
    for donor in range(8):
        for recipient in range(8):
            if donor == recipient:
                continue
            for amount in (1, 2, 3):
                value = list(reference)
                value[donor] -= amount
                value[recipient] += amount
                if all(value[4 * side:4 * side + 4] == sorted(value[4 * side:4 * side + 4])
                       for side in range(2)):
                    targets.add(tuple(value))
    return targets


def pressure(target, shifts, heights):
    def log_partition(energies):
        values = [-SCALE * value for value in energies]
        maximum = max(values)
        return maximum + math.log(math.fsum(math.exp(value - maximum) for value in values))
    physical = (heights[0] + heights[3] + heights[4] + heights[7],
                heights[1] + heights[2] + heights[5] + heights[6])
    cuts = math.fsum(log_partition([target[4 * side + i] - shifts[side] for i in range(4)])
                     for side in range(2))
    return cuts - math.sqrt(2) * math.fsum(
        log_partition((0, energy / math.sqrt(2))) for energy in physical)


def line_points(engine, pivots):
    free = [index for index in range(engine.lattice.variables) if index not in pivots]
    assert len(free) == 1
    points = []
    for color in range(4):
        planes = ((color & 1) << free[0], (color >> 1) << free[0])
        for pivot, (low, high, rhs) in sorted(pivots.items()):
            assert not ((planes[0] | planes[1]) >> pivot) & 1
            value = dot((low, high), planes) ^ rhs
            planes = planes[0] | ((value & 1) << pivot), planes[1] | ((value >> 1) << pivot)
        assert all(dot((low, high), planes) == rhs for low, high, rhs in pivots.values())
        points.append(engine.lattice.seed(*planes))
    return points


def shift_coefficients(values, exponent):
    if exponent >= 0:
        return tuple((low << exponent, high << exponent) for low, high in values)
    mask = (1 << -exponent) - 1
    assert all(not ((low | high) & mask) for low, high in values)
    return tuple((low >> -exponent, high >> -exponent) for low, high in values)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path,
                        default=root / 'evidence' / 'binary-two-cut-zero-sum-target-box.json')
    parser.add_argument('--output', type=Path,
                        default=root / 'evidence' / 'binary-target-box-replay.json')
    args = parser.parse_args()
    record = json.loads(args.record.read_text())
    assert record['degree'] == 12 and record['baseline_laurent_clearing'] == 4
    assert record['maximum_row_target_transfer'] == 3 and record['zero_sum_target_box_radius'] == 2
    assert record['random_seed'] == 310781 and record['complete_requested_targets'] is True
    engine = FlagSearch(12, 4)
    lattice = engine.lattice
    reference = lattice.orders
    assert reference == (5, 8, 10, 13, 5, 8, 10, 13)
    raw = independent_targets(reference)
    assert len(raw) == record['raw_sorted_target_count'] == 23099
    baseline_pressure = pressure(reference, lattice.scalar_shifts, lattice.heights)
    kept = {target for target in raw
            if pressure(target, lattice.scalar_shifts, lattice.heights) <= baseline_pressure + 1e-12}
    assert len(kept) == record['target_count'] == 7079
    variants = {tuple(row['targets']): row for row in record['variants']}
    assert len(variants) == len(record['variants']) == 7079 and set(variants) == kept
    assert tuple(record['variants'][0]['targets']) == reference

    frames = flags()
    identity = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
    identity_index = frames.index(identity)
    assert len(frames) == record['available_flags'] == 8925
    requested = random.Random(record['random_seed']).sample(range(len(frames)), 128)
    selected = sorted(set(requested + [identity_index]))
    assert selected == record['flag_indices'] and len(selected) == 129

    # Common constraints depend only on the two minima, not on any flag.
    common_cache = {}
    common_failures = 0
    for target, row in variants.items():
        key = min(target[:4]), min(target[4:])
        if key not in common_cache:
            common_cache[key] = engine.common_system(target)
        common = common_cache[key]
        assert row['common_constraint_rank'] == (None if common is None else len(common))
        assert row['flag_pairs_decided'] == 129 ** 2
        survivors = row['separately_feasible_cut_flags']
        assert len(survivors) == 2 and all(0 <= count <= 129 for count in survivors)
        assert 0 <= row['feasible_pairs'] <= survivors[0] * survivors[1]
        assert row['profile_queries'] == 8 * row['feasible_pairs']
        if common is None:
            assert survivors == [0, 0] and row['feasible_pairs'] == row['profile_queries'] == 0
            common_failures += 1
    assert common_failures == 1430 and len(common_cache) == 17

    feasible = {target for target, row in variants.items() if row['feasible_pairs']}
    expected = {tuple(value + (gain if i < 4 else -gain) for i, value in enumerate(reference)): gain
                for gain in range(-2, 3)}
    assert feasible == set(expected)
    families = {}
    fibers = []
    profile_checks = endpoint_coefficients = 0
    for target, gain in sorted(expected.items(), key=lambda item: item[1]):
        row = variants[target]
        assert row['feasible_pairs'] == 1 and row['separately_feasible_cut_flags'] == [1, 1]
        common = common_cache[min(target[:4]), min(target[4:])]
        pivots = dict(common)
        for side in range(2):
            system = engine.system(side, identity, target[4 * side:4 * side + 4], common)
            assert system is not None
            for low, high, rhs in system.values():
                assert insert(pivots, (low | ((rhs & 1) << lattice.variables),
                                       high | ((rhs >> 1) << lattice.variables)), lattice.variables)
        assert lattice.variables - len(pivots) == 1
        points = line_points(engine, pivots)
        expected_profiles = tuple(tuple(target[4 * side + i] - lattice.scalar_shifts[side]
                                         for i in range(4)) for side in range(2))
        family = set()
        for seed in points:
            values = seed['values']
            endpoint(values)
            endpoint_coefficients += 16
            assert literal_profiles(values) == expected_profiles
            profile_checks += 2
            # Increasing the infinity target by gain means multiplying
            # the source by w^-gain.  Multiplication by w^gain therefore
            # returns its complete coefficient table to the reference.
            family.add(shift_coefficients(values, gain))
        assert len(family) == 4
        families[gain] = family
        fibers.append({'infinity_target_shift': gain, 'target': target,
                       'affine_dimension_over_F4': 1, 'all_points_replayed': 4,
                       'complete_cut_profiles': expected_profiles})
    assert all(family == families[0] for family in families.values())
    original = TwoCutLattice(5)
    original_solution = original.solve(original.orders)
    assert original_solution is not None
    original_free = [index for index in range(original.variables) if index not in original_solution]
    assert len(original_free) == 1
    original_family = {
        shift_coefficients(original.seed(
            original.vector(original_solution, (color & 1) << original_free[0]),
            original.vector(original_solution, (color >> 1) << original_free[0], True))['values'], 4)
        for color in range(4)}
    assert original_family == families[0]
    assert len(record['profiles']) == 5
    represented = set()
    for row in record['profiles']:
        target = tuple(row['target'])
        assert target in expected and row['dimension'] == 1
        assert row['flag_indices'] == [identity_index, identity_index]
        values = tuple(map(tuple, row['values']))
        endpoint(values)
        endpoint_coefficients += 16
        assert literal_profiles(values) == tuple(map(tuple, row['profiles']))
        profile_checks += 2
        assert shift_coefficients(values, expected[target]) in families[0]
        represented.add(target)
    assert represented == feasible
    result = {
        'status': 'PASS: independent target reconstruction, exact common rejections, and complete scalar-shift fibers.',
        'record': str(args.record.relative_to(root)),
        'degree_after_clearing': 12, 'common_Laurent_clearing': 4,
        'raw_independently_enumerated_targets': len(raw), 'retained_fixed_pressure_targets': len(kept),
        'pressure_filter_scale': SCALE, 'pressure_filter_tolerance': 1e-12,
        'flags_per_cut': len(selected), 'declared_completion_checked': True,
        'distinct_common_minimum_systems_replayed': len(common_cache),
        'exact_common_inconsistent_targets': common_failures,
        'feasible_target_fibers': fibers,
        'normalized_four_point_families_exactly_equal': True,
        'also_equal_to_original_degree_five_line_after_common_clearing': True,
        'complete_weighted_cut_profiles_by_all_minors': profile_checks,
        'literal_endpoint_coefficients': endpoint_coefficients,
        'scope': 'The floating pressure filter is replayed as a finite schedule. Exact algebra verifies all common failures and every point of the five recorded feasible lines. The remaining flag-pair exclusions are manifest-checked, not independently rerun.',
        'source_sha256': {str(path.relative_to(root)): sha256(path.read_bytes()).hexdigest()
                          for path in (args.record, Path(__file__).resolve(),
                                       root / 'tools' / 'search_binary_f4_flags.py',
                                       root / 'tools' / 'search_binary_two_cut_lattices.py',
                                       root / 'tools' / 'check_binary_two_cut_lattices.py')},
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('source_sha256', 'feasible_target_fibers')}))


if __name__ == '__main__':
    main()
