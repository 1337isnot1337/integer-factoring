#!/usr/bin/env python3
"""Literal coherent hot-code construction and exact rho <18/25 certificate."""
from __future__ import annotations
import argparse
from hashlib import sha256
import json
from pathlib import Path
from certify_hot_charge_pressure import homogeneous_sum, integer_digest
from coherent_hot_boundary import compress
from hot_code_charge_profiles import hot_energy
from hot_compound_crossing import compound_crossing
from native_matrix_smith_f4 import MatrixOracle
from native_rectangular_smith_f4 import RectangularOracle
from physical_split_crossings import flatten, retract
from reflection_laurent_physical_frames import digest
from restrict_implicit_reflection import restricted_tensor


def logical(mask):
    return sum(1 << (7 - a) for a in range(8) if mask >> a & 1)


def full_maps(low, boundary, left, right):
    fixed = [[int(i == j) for j in range(low + boundary)] for i in range(low)]
    return (fixed + [[0] * low + list(row) for row in left],
            fixed + [[0] * low + list(row) for row in right])


def natural_profiles(seed, energies):
    q = seed['q_x']
    result = []
    for cut in (1, 0):
        matrix, rows, columns = flatten(seed, cut)
        if cut == 1:
            rw = [energies[u] / 2 for x, u in rows]
            cw = [energies[v] / 2 for y, v in columns]
        else:
            rw = [(energies[x] + energies[y]) / 2 for x, y in rows]
            cw = [0] * len(columns)
        orders = MatrixOracle(matrix, q * q).profile(rw, cw, infinity=bool(cut))['orders']
        assert all(value == round(value) for value in orders)
        result.append(sorted(map(int, orders)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    lower, boundary = [0, 1, 2, 3, 4, 5, 8], [6, 9]
    left, right = [[1, 1]], [[0, 1]]
    letters = lower + boundary
    raw = restricted_tensor(16, list(map(logical, letters)))
    seed = compress(raw, len(lower), left, right)
    full_left, full_right = full_maps(len(lower), len(boundary), left, right)
    independent = retract(raw, full_left, full_right)
    assert seed['values'] == independent['values']
    content = min((part & -part).bit_length() - 1
                  for value in seed['values'] for part in value if part)
    assert seed['values'][0] == (1 << 128, 0) and content == 104
    seed['values'] = [(a >> content, b >> content) for a, b in seed['values']]
    seed['degree'] = max(max(a.bit_length(), b.bit_length()) - 1 for a, b in seed['values'])
    assert seed['degree'] == 56
    assert digest(seed['values']) == '17254d71e2fece406c199fba1988424f3619911d5cc93f01c900d43a5c8bdbf9'
    energies = list(map(hot_energy, lower)) + [32]
    literal = natural_profiles(seed, energies)
    height = [value / 2 for value in energies[1:]]
    params = [0] * 14 + height + height + height + height + [0] * 14 + [1]
    native = RectangularOracle(seed['values'], 8, 8).prefixes(params)['orders']
    assert [sorted(native[:64]), sorted(native[64:])] == literal
    base = compound_crossing(letters)
    regular, coefficient_controls = [], 4096
    for swap in (True, False):
        compressed = compress(base, len(lower), left, right, swap_y=swap)
        reference = (retract(base, full_left, full_right, full_right, full_left)
                     if swap else retract(base, full_left, full_right))
        assert compressed['values'] == reference['values']
        coefficient_controls += len(reference['values'])
        matrix, _, _ = flatten(compressed, 0)
        regular.append(sorted(map(int, MatrixOracle(matrix, 64).profile(scale=2)['orders'])))
    assert literal == [[v - 24 for v in regular[0]], [v + 24 for v in regular[1]]]
    # A boundary at energy 36 mixes a singleton and a three-particle
    # state. Its full matrix, rather than old occupation blocks, is used.
    all_letters = sorted(range(256), key=lambda x: (hot_energy(x), x))
    mixed_lower = [x for x in all_letters if hot_energy(x) < 36]
    mixed_boundary = [x for x in all_letters if hot_energy(x) == 36]
    assert {x.bit_count() for x in mixed_boundary} == {1, 3}
    mixed_letters = mixed_lower + mixed_boundary
    mixed_raw = restricted_tensor(16, list(map(logical, mixed_letters)))
    mixed_seed = compress(mixed_raw, len(mixed_lower), [[1, 1]], [[1, 0]])
    mixed_literal = natural_profiles(mixed_seed, list(map(hot_energy, mixed_lower)) + [36])
    mixed_base = compound_crossing(mixed_letters)
    for index, swap in enumerate((True, False)):
        changed = compress(mixed_base, len(mixed_lower), [[1, 1]], [[1, 0]], swap_y=swap)
        matrix, _, _ = flatten(changed, 0)
        expected = sorted(map(int, MatrixOracle(matrix, 100).profile(scale=2)['orders']))
        assert mixed_literal[index] == [v + (-128 if swap else 128) for v in expected]
    # Equal-energy compatibility is necessary for this unchanged-gauge
    # shortcut. A split mixing energies 4 and 12 remains a strict tensor,
    # but declaring energy 12 for the output does not commute with it.
    bad_raw = restricted_tensor(16, list(map(logical, [0, 1, 2])))
    bad_seed = compress(bad_raw, 1, [[1, 1]], [[0, 1]])
    bad_literal = natural_profiles(bad_seed, [0, 12])[1]
    bad_base = compound_crossing([0, 1, 2])
    bad_compressed = compress(bad_base, 1, [[1, 1]], [[0, 1]])
    bad_matrix, _, _ = flatten(bad_compressed, 0)
    bad_shortcut = [int(v) + 128 for v in sorted(MatrixOracle(bad_matrix, 4).profile(scale=2)['orders'])]
    assert bad_literal != bad_shortcut
    A, B, numerator, denominator = 36, 25, 141, 142
    p0, d0 = homogeneous_sum([A * x for x in regular[0]], numerator, denominator)
    p1, d1 = homogeneous_sum([A * x for x in regular[1]], numerator, denominator)
    Q, dq = homogeneous_sum([B * x for x in energies], numerator, denominator)
    lhs, rhs = (p0 * p1) ** B, Q ** (2 * A)
    difference = B * (d0 + d1) - 2 * A * dq
    if difference >= 0:
        rhs *= denominator ** difference
    else:
        lhs *= denominator ** (-difference)
    assert 0 < lhs < rhs and lhs.bit_length() == rhs.bit_length() == 720752
    names = (Path(__file__).name, 'coherent_hot_boundary.py', 'hot_compound_crossing.py',
             'hot_code_charge_profiles.py', 'hot_code_charge_profiles.cpp',
             'physical_split_crossings.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.py',
             'native_rectangular_smith_f4.cpp', 'restrict_implicit_reflection.py',
             'certify_hot_charge_pressure.py', 'reflection_dual_rail_gaussian.py',
             'check_schur_gaussian_reduction.py', 'selected_gaussian_coefficients.cpp',
             'native_weighted_smith_f4.py')
    report = {'scope': 'Explicit coherent eight-letter crossing and exact uniform rho <18/25; no factoring score gain.',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'recipe': {'interval': 16, 'lower_letters': lower, 'boundary_letters': boundary,
                         'boundary_left_rows': left, 'boundary_right_columns': right},
              'seed': seed, 'monomial_content': content, 'tensor_sha256': digest(seed['values']),
              'energies': energies, 'literal_profiles_infinity_zero': literal,
              'regular_profiles_infinity_zero': regular,
              'certificate': {'base': [numerator, denominator], 'support_parameter': [A, B],
                              'balanced_strict_bound': [18, 25], 'left_integer_bits': lhs.bit_length(),
                              'right_integer_bits': rhs.bit_length(), 'left_sha256': integer_digest(lhs),
                              'right_sha256': integer_digest(rhs), 'positive_gap_sha256': integer_digest(rhs - lhs),
                              'type_count_exponents': [[B * value, dq - B * value] for value in energies],
                              'type_count_formula': '141^first * 142^second', 'type_count_sum_hex': hex(Q)},
              'counts': {'independent_complete_coefficients': coefficient_controls,
                         'independent_native_profile_channels': 128,
                         'literal_to_compound_profile_channels': 328,
                         'mixed_occupation_controls': 1, 'unequal_energy_countercontrols': 1},
              'unequal_energy_countercontrol': {'literal_orders': bad_literal, 'invalid_shortcut': bad_shortcut},
              'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'bound': [18, 25], 'integer_bits': lhs.bit_length(), 'counts': report['counts']}))


if __name__ == '__main__':
    main()
