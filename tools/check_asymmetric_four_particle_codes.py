#!/usr/bin/env python3
"""Independent exact controls for two-alphabet Gaussian code searches.

Every retained coefficient is expanded as a determinant of the closed
strict kernel. Complete profiles are certified with regular polynomial
frames, independently of the native minor constructor and Smith oracle.
"""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from asymmetric_four_particle_codes import construct
from check_aligned_gaussian_complements import expand_square, minor, multiply, strict_matrix
from check_f4_recentered_walks import certify_profile, literal_matrix, same_cost
from replay_source_provenance import resolve_source_hashes
from search_f4_recentered_crossings import digest, five_letter_source, root_for_profiles


def grade(mask):
    return sum(2 * i + 1 for i in range(8) if mask >> i & 1)


def canonical(code):
    return tuple(sorted(code, key=lambda mask: (grade(mask), mask)))


def certify_tensor(row):
    code_x, code_y = row['code_x'], row['code_y']
    for code, energy in ((code_x, row['energies_x']), (code_y, row['energies_y'])):
        assert len(code) == len(set(code)) == 5
        assert all(isinstance(mask, int) and 0 <= mask < 256 and mask.bit_count() == 4 for mask in code)
        assert tuple(code) == canonical(code)
        grades = list(map(grade, code))
        assert tuple(energy) == tuple(value - min(grades) for value in grades)
    common, values = row['entry_gcd_in_z'], row['values']
    assert common > 0 and common.bit_count() % 2 == 1
    assert row['entry_gcd_degree'] == common.bit_length() - 1
    assert len(values) == 625 and row['degree'] == max(low.bit_length() - 1 for low, high in values)
    matrix = strict_matrix()
    for index, (polynomial, high) in enumerate(values):
        x, y, u, v = ((index // 5 ** leg) % 5 for leg in range(4))
        assert high == 0 and polynomial.bit_count() % 2 == int(x == v and y == u)
        actual = minor(matrix, code_x[v] | (code_y[y] << 8), code_y[u] | (code_x[x] << 8))
        expected = expand_square(multiply(common, polynomial))
        assert actual == expected
    assert digest(values) == row['digest']
    hx = tuple(value // 2 for value in row['energies_x'])
    hy = tuple(value // 2 for value in row['energies_y'])
    zero = (0,) * 5
    heights = (zero, zero, hy, hx, hx, hy, zero, zero)
    for side in range(2):
        literal, scalar = literal_matrix(values, 5, row['degree'], heights, side)
        certify_profile(literal, scalar, row['profiles'][side])
    same_cost(root_for_profiles(row['profiles'], row['energies_x'], row['energies_y']), row)
    return 625, 2


def replay(path, prior=None):
    search = json.loads(path.read_text())
    provenance, unresolved = resolve_source_hashes(search['source_sha256'])
    assert not unresolved, unresolved
    states = tuple(sum(1 << i for i in positions) for positions in combinations(range(8), 4))
    reference = (15, 23, 43, 51, 77)
    assert len(states) == 70
    if prior is None:
        planned = [(reference, reference, 'reference')]
        for removed in reference:
            kept = tuple(mask for mask in reference if mask != removed)
            planned.extend((reference, canonical(kept + (replacement,)), 'one_paired_symbol_replacement')
                           for replacement in states if replacement not in reference)
        assert len(planned) == 326
        assert search['planned_schedule'] == [[list(x), list(y), stage] for x, y, stage in planned]
        assert len(search['results']) == search['complete_tensor_queries'] == len(planned)
        for row, (code_x, code_y, stage) in zip(search['results'], planned):
            assert tuple(row['code_x']) == code_x and tuple(row['code_y']) == code_y and row['stage'] == stage
        complete_rows = search['results']
    else:
        prior_path, prior_search = prior
        assert search['reused_record_sha256'] == sha256(prior_path.read_bytes()).hexdigest()
        ordered = sorted(states, key=lambda mask: (grade(mask), mask))[:search['arguments']['states']]
        assert ordered == search['available_shell_states']
        planned = [canonical(code) for code in combinations(ordered, 5)]
        assert search['complete_shell_schedule'] == list(map(list, planned))
        previous = {tuple(row['code_y']): row for row in prior_search['results']}
        reused = [previous[code] for code in planned if code in previous]
        fresh = [code for code in planned if code not in previous]
        assert search['reused_results'] == reused
        assert search['fresh_schedule'] == list(map(list, fresh))
        assert [row['code_y'] for row in search['results']] == list(map(list, fresh))
        assert len(search['results']) == search['complete_new_tensor_queries'] == len(fresh)
        assert search['covered_shell_codes'] == len(planned) == len(reused) + len(fresh)
        complete_rows = reused + search['results']
    assert search['completed_requested_schedule']
    assert len({tuple(row['code_y']) for row in complete_rows}) == len(complete_rows)
    assert all(tuple(row['code_x']) == reference for row in complete_rows)
    assert abs(min(row['rho'] for row in complete_rows) - search['best_recorded_rho']) < 1e-10
    assert search['strict_pressure_candidate_found'] == any(row['critical_pressure'] < -1e-9 for row in complete_rows)
    lookup = {tuple(row['code_y']): row for row in search['results']}
    minors = profiles = 0
    seen = set()
    for row in search['retained']:
        key = tuple(row['code_y'])
        assert key not in seen
        seen.add(key)
        assert {name: value for name, value in row.items() if name != 'values'} == lookup[key]
        m, p = certify_tensor(row)
        minors += m
        profiles += p
    # Keep a nonreflection control even if the lowest costs are all admissible.
    additional = min((row for row in complete_rows if row['nonreflection_y_states']), key=lambda row: row['rho'])
    if tuple(additional['code_y']) not in seen:
        exact = construct(additional['code_x'], additional['code_y'])
        extra = additional | {'values': exact['values']}
        m, p = certify_tensor(extra)
        minors += m
        profiles += p
    return search, {'record': str(path), 'record_sha256': sha256(path.read_bytes()).hexdigest(),
                    'status': 'passed', 'source_provenance': provenance,
                    'complete_declared_alphabets': len(complete_rows),
                    'retained_literal_tensors': len(seen),
                    'independent_complete_polynomial_minors': minors,
                    'complete_profiles_certified_by_polynomial_frames': profiles,
                    'best_recorded_rho': search['best_recorded_rho'],
                    'best_nonreflection_control': {key: additional[key] for key in ('code_y', 'degree', 'digest', 'rho')},
                    'strict_pressure_candidate_found': search['strict_pressure_candidate_found']}


def main():
    first = Path('evidence/asymmetric-four-particle-codes-q5.json')
    shell = Path('evidence/asymmetric-four-particle-shell16-q5.json')
    previous, one_result = replay(first)
    _, shell_result = replay(shell, (first, previous))
    reference = five_letter_source()[0]
    base = previous['retained'][0]
    original_order = (15, 77, 43, 23, 51)
    permutation = [original_order.index(mask) for mask in base['code_x']]
    for index, pair in enumerate(base['values']):
        old_index = sum(permutation[(index // 5 ** leg) % 5] * 5 ** leg for leg in range(4))
        assert tuple(pair) == reference[old_index]
    print(json.dumps({'status': 'exact asymmetric Gaussian controls passed',
                      'baseline_permutation_identities': 625,
                      'records': [one_result, shell_result],
                      'factoring_bound_improved': False}, indent=2))


if __name__ == '__main__':
    main()
