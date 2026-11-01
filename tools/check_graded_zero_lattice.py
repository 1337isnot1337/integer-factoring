#!/usr/bin/env python3
"""Independent complete fraction-free profiles for the graded zero construction."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_joint_energy_conjugations import cut_matrix
from graded_zero_lattice import construct, positive_dress
from joint_energy_conjugations import model
from search_f4_recentered_crossings import five_letter_source
from weighted_smith_oracle_f4 import WeightedSmithOracle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--numerical', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    seed = construct()
    source, _, ex, ey = five_letter_source()
    denominator, clearing = seed['ramification'], seed['common_monomial_clearing']
    data = model()
    infinity_weights = [denominator*(ex[i//5]//2) for i in range(25)]
    infinity, _ = WeightedSmithOracle(cut_matrix(seed['values'], True), infinity=True).profile(
        infinity_weights, infinity_weights)
    zero, _ = WeightedSmithOracle(cut_matrix(seed['values'], False)).profile(
        [denominator*energy for energy in data['energies']], [0]*25)
    normalized_zero = [(value-clearing)/denominator for value in zero]
    normalized_infinity = [(value+clearing)/denominator for value in infinity]
    expected = [8,10,10,12,13,14,14,16,16,17,18,18,18,19,19,20,20,20,21,22,22,22,23,24,28]
    assert normalized_zero == expected
    assert list(infinity) == [-20,-20,-20,-19,-19,-19,-19,-19,-19,-18,-18,-18,-17,-17,-17,-17,-16,-16,-15,-15,-11,-9,-2,4,6]
    assert sum(zero) == 444*denominator + 25*clearing
    dressings = []
    for energy, power in ((1, 0), (9, 2)):
        block = next(block for block in data['blocks'] if data['energies'][block[0]] == energy)
        generator = [[0]*25 for _ in range(25)]
        local = ((2, 1), (3, 2))
        for i, r in enumerate(block):
            for j, c in enumerate(block):
                generator[r][c] = local[i][j]
        dressed = positive_dress(seed['values'], generator, power)
        actual, _ = WeightedSmithOracle(cut_matrix(dressed, False)).profile(
            [denominator*value for value in data['energies']], [0]*25)
        assert actual == zero
        dressings.append({'energy_block': energy, 'generator': local, 'power': power,
                          'independent_complete_zero_orders': actual})
    report = {'status': 'pass', 'ramification': denominator, 'common_monomial_clearing': clearing,
              'polynomial_degree': seed['polynomial_degree'], 'basis_columns': seed['basis_columns'],
              'joint_potentials': seed['ramified_joint_potentials'],
              'changed_energy_blocks': seed['changed_energy_blocks'],
              'all_strict_endpoint_entries_checked': 625,
              'independent_fraction_free_integer_orders': 100,
              'complete_integer_zero_orders': zero, 'complete_integer_infinity_orders': infinity,
              'normalized_zero_orders': normalized_zero,
              'normalized_infinity_orders': normalized_infinity,
              'polynomial_zero_preserving_dressings': dressings,
              'scope': 'An explicit finite polynomial crossing with an improved zero profile and '
              'a worse opposite profile. The complete paired cost does not improve factoring.'}
    if args.numerical:
        from search_f4_recentered_crossings import root_for_profiles
        original_weights = [ex[i//5]//2 for i in range(25)]
        original_infinity, _ = WeightedSmithOracle(cut_matrix(source, True), infinity=True).profile(
            original_weights, original_weights)
        report['actual_numerical_pressure'] = root_for_profiles(
            (normalized_infinity, normalized_zero), ex, ey)
        report['unrealized_reference_infinity_comparison'] = root_for_profiles(
            (original_infinity, normalized_zero), ex, ey)
        report['comparison_limit'] = 'The reference infinity profile belongs to a different kernel and cannot be substituted in the actual cost.'
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/check_graded_zero_lattice.py', 'tools/graded_zero_lattice.py',
               'tools/joint_energy_conjugations.py', 'tools/check_joint_energy_conjugations.py',
               'tools/weighted_smith_oracle_f4.py', 'tools/weighted_pascal_flags.py',
               'tools/tensor_pascal_box_crossings.py', 'tools/search_f4_recentered_crossings.py']
    sources += ['tools/f4_pascal_flags.py', 'tools/cube_f4_bivariate_smith.py']
    report['source_sha256'] = {name: sha256((root / name).read_bytes()).hexdigest() for name in sources}
    report['polynomial_seed_sha256'] = sha256(json.dumps(seed['values']).encode()).hexdigest()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass', 'independent_complete_orders': 100,
                      'zero_profile_start': normalized_zero[:3], 'polynomial_degree': seed['polynomial_degree']}))


if __name__ == '__main__':
    main()
