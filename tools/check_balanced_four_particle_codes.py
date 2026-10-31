#!/usr/bin/env python3
"""Independently certify balanced-code polynomial tensors and Smith profiles.

The tensor replay expands determinants in the closed 16 by 16 strict kernel,
without calling the native selected-minor constructor.  The profile replay
uses exact regular polynomial frames, without the native Smith oracle.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import time

from balanced_four_particle_codes import REFERENCE_CODE, construct
from check_aligned_gaussian_complements import expand_square, minor, multiply, strict_matrix
from check_f4_recentered_walks import certify_profile, literal_matrix, same_cost
from replay_source_provenance import resolve_source_hashes
from search_f4_recentered_crossings import digest, five_letter_source, root_for_profiles
from search_saturated_dual_directions import complete_profiles


def certify_tensor(record):
    code, values, common = record['code_masks'], record['values'], record['entry_gcd_in_z']
    assert len(code) == len(set(code)) == 5
    assert all(isinstance(mask, int) and 0 <= mask < 256 and mask.bit_count() == 4 for mask in code)
    grades = [sum(2*i+1 for i in range(8) if (mask >> i) & 1) for mask in code]
    assert list(record['paired_energies']) == [value-min(grades) for value in grades]
    assert common > 0 and common.bit_count() % 2 == 1
    assert common.bit_length()-1 == record['entry_gcd_degree']
    assert len(values) == 625 and max(low.bit_length()-1 for low, high in values) == record['degree']
    matrix = strict_matrix()
    for index, (polynomial, high) in enumerate(values):
        x, y, u, v = ((index//5**leg) % 5 for leg in range(4))
        assert high == 0 and polynomial.bit_count() % 2 == int(x == v and y == u)
        expected = expand_square(multiply(common, polynomial))
        actual = minor(matrix, code[v] | (code[y] << 8), code[u] | (code[x] << 8))
        assert actual == expected
    if 'digest' in record:
        assert digest(values) == record['digest']
    return 625


def certify_record_profiles(record):
    energy, degree = record['paired_energies'], record['degree']
    half = tuple(value//2 for value in energy)
    zero = (0,)*5
    heights = (zero, zero, half, half, half, half, zero, zero)
    for side in range(2):
        matrix, scalar = literal_matrix(record['values'], 5, degree, heights, side)
        certify_profile(matrix, scalar, record['profiles'][side])
    same_cost(root_for_profiles(record['profiles'], energy, energy), record)
    return 2


def controls():
    # The second alphabet has a repeated minimum grade. The final two have
    # no reflection-admissibility requirement and reach high physical modes.
    codes = (REFERENCE_CODE, (27, 39, 43, 51, 77),
             (27, 45, 54, 99, 170), (240, 232, 228, 226, 225))
    rows = []
    minors = profiles = 0
    for code in codes:
        tensor = construct(code)
        tensor['profiles'] = complete_profiles(tensor['values'], {'natural': tensor['paired_energies']})['natural']
        tensor.update(root_for_profiles(tensor['profiles'], tensor['paired_energies'], tensor['paired_energies']))
        minors += certify_tensor(tensor)
        profiles += certify_record_profiles(tensor)
        rows.append({key: value for key, value in tensor.items() if key != 'values'})
    baseline = construct(REFERENCE_CODE)
    reference = five_letter_source()[0]
    permutation = [REFERENCE_CODE.index(mask) for mask in baseline['code_masks']]
    for index, pair in enumerate(baseline['values']):
        old_index = sum(permutation[(index//5**leg) % 5]*5**leg for leg in range(4))
        assert pair == reference[old_index]
    rejected = 0
    for code in ((15, 15, 23, 43, 51), (1, 23, 43, 51, 77), (15, 23, 43, 51)):
        try:
            construct(code)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('An invalid physical alphabet was accepted.')
    return {'scope': 'Independent complete strict-kernel determinants and polynomial Smith frames.',
            'status': 'passed', 'control_alphabets': len(codes),
            'independent_complete_polynomial_minors': minors,
            'complete_profiles_certified_by_polynomial_frames': profiles,
            'baseline_permuted_coefficient_identities': 625,
            'invalid_alphabets_rejected': rejected, 'controls': rows}


def replay(path):
    search = json.loads(path.read_text())
    provenance, unresolved = resolve_source_hashes(search['source_sha256'])
    assert not unresolved, unresolved
    results, schedule = search['results'], search['planned_schedule']
    assert len(results) == search['complete_tensor_queries'] <= len(schedule)
    assert search['exact_endpoint_entries'] == 625*len(results)
    assert all(row['code_masks'] == item[0] and row['stage'] == item[1]
               for row, item in zip(results, schedule))
    assert len({tuple(row['code_masks']) for row in results}) == len(results)
    assert abs(min(row['rho'] for row in results)-search['best_recorded_rho']) < 1e-10
    if search['completed_requested_schedule']:
        assert len(results) == len(schedule)
    assert search['strict_pressure_candidate_found'] == any(row['critical_pressure'] < -1e-9 for row in results)
    minors = profiles = 0
    seen = set()
    lookup = {tuple(row['code_masks']): row for row in results}
    for row in search['retained']:
        key = tuple(row['code_masks'])
        assert key not in seen
        seen.add(key)
        assert {k: v for k, v in row.items() if k != 'values'} == lookup[key]
        minors += certify_tensor(row)
        profiles += certify_record_profiles(row)
    assert abs(min(row['rho'] for row in search['retained'])-search['best_recorded_rho']) < 1e-10
    return {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(),
            'status': 'passed', 'source_provenance': provenance,
            'completed_requested_schedule': search['completed_requested_schedule'],
            'recorded_complete_tensor_queries': len(results),
            'retained_literal_tensors': len(seen),
            'independent_complete_polynomial_minors': minors,
            'complete_profiles_certified_by_polynomial_frames': profiles,
            'best_recorded_rho': search['best_recorded_rho'],
            'strict_pressure_candidate_found': search['strict_pressure_candidate_found']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--search', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    report = replay(args.search) if args.search else controls()
    report.update(checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
