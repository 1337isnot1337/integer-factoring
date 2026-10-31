#!/usr/bin/env python3
"""Exact polynomial identities, image ranks, and profile controls for saturation."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_block_smith_crossings import binary_rank
from check_f4_recentered_walks import certify_profile, literal_matrix
from saturated_dual_kernel_directions import literal_directions
from search_f4_recentered_crossings import five_letter_source
from search_saturated_dual_directions import GAUGES, check_endpoint, complete_profiles, specialize


ROOT = Path(__file__).resolve().parent.parent


def multiply(first, second):
    result = 0
    while second:
        bit = second & -second
        second ^= bit
        result ^= first << (bit.bit_length()-1)
    return result


def check(path):
    record = json.loads(path.read_text())
    assert record['dimensions'] == [5]*4 and record['degree'] == 18
    divisions = 0
    if record['constructor'] == 'saturated_bounded_dual_kernel_directions':
        source_path = ROOT/record['source_record']
        assert sha256(source_path.read_bytes()).hexdigest() == record['source_record_sha256']
        source = json.loads(source_path.read_text())
        for primitive in record['primitives']+record['excluded_by_degree']:
            old = source['kernel_slabs'][primitive['source_slab']]
            new = primitive['normalized_slab']
            common = primitive['entry_gcd']
            assert common > 0 and len(old) == len(new) == 125
            for a, b in zip(old, new):
                assert multiply(common, b) == a ^ (a << 1)
                divisions += 1
            assert max(value.bit_length()-1 for value in new) == primitive['normalized_degree']
    directions = literal_directions(record)
    expected = []
    for primitive, logical, shift in record['direction_descriptors']:
        slab = record['primitives'][primitive]['normalized_slab']
        assert 0 <= logical < 5 and shift >= 0 and len(slab) == 125
        direction = []
        for index in range(625):
            value = (slab[index//5] << shift) if index % 5 == logical else 0
            assert isinstance(value, int) and 0 <= value < 1 << 19
            assert value.bit_count() % 2 == 0
            direction.append((value, 0))
        expected.append(tuple(direction))
    assert tuple(expected) == directions
    encoded = [sum(low << (19*entry) for entry, (low, high) in enumerate(direction)) for direction in directions]
    rank = binary_rank(encoded)
    assert rank == len(directions)
    sparse = [tuple((entry, low) for entry, (low, high) in enumerate(direction) if low) for direction in directions]
    base = five_letter_source()[0]
    assignments = [((0, 1),), ((len(directions)-1, 2),),
                   tuple((i, c) for i, c in zip(sorted({0, len(directions)//2, len(directions)-1}), (1, 2, 3)))]
    profiles_checked = conjugate_profiles = 0
    for assignment in assignments:
        values = specialize(base, sparse, assignment)
        check_endpoint(values, 18)
        profiles = complete_profiles(values, GAUGES)
        conjugate = tuple((a ^ b, b) for a, b in values)
        assert complete_profiles(conjugate, GAUGES) == profiles
        conjugate_profiles += 4
        for name, energy in GAUGES.items():
            half, zero = tuple(value//2 for value in energy), (0,)*5
            heights = (zero, zero, half, half, half, half, zero, zero)
            for side in range(2):
                matrix, scalar = literal_matrix(values, 5, 18, heights, side)
                certify_profile(matrix, scalar, profiles[name][side])
                profiles_checked += 1
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(),
            'status': 'passed', 'exact_polynomial_normalization_identities': divisions,
            'endpoint_zero_directions': len(directions), 'direction_entries_checked': 625*len(directions),
            'independent_binary_rank': rank, 'dimension_after_F4_scalar_extension': rank,
            'literal_affine_endpoint_entries': 625*len(assignments),
            'complete_profiles_certified_by_polynomial_frames': profiles_checked,
            'Frobenius_profile_comparisons': conjugate_profiles}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('records', type=Path, nargs='+')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    results = []
    for path in args.records:
        item = check(path)
        results.append(item)
        print(json.dumps(item), flush=True)
    if args.output:
        args.output.write_text(json.dumps({'scope': 'Exact polynomial normalization, independent direction rank, and literal complete-profile controls.',
                                         'records': results, 'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2)+'\n')


if __name__ == '__main__':
    main()
