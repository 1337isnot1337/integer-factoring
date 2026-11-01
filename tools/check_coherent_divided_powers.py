#!/usr/bin/env python3
"""Literal derivative ranks, symplectic identities and unequal-cut controls."""
import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path

from check_binary_coefficient_entropy import binary_rank
from check_cancellation_module_automata import pmul
from check_coherent_hot_compounds import determinant
from check_joint_energy_conjugations import polynomial
from coherent_divided_powers import disjoint_rank, graded_dimensions, paired_top
from coherent_hot_compounds import height, one_particle, profiles, support
from coherent_inverse_systems import derivative, inverse_system
from coherent_unequal_compounds import profiles as unequal_profiles
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]


def certificate():
    disjoint_cases = 0
    for n in range(10):
        for a in range(n+1):
            rows = [sum(1 << i for i in c) for c in combinations(range(n), a)]
            for b in range(n-a+1):
                columns = [sum(1 << i for i in c) for c in combinations(range(n), b)]
                matrix = [sum(1 << j for j, column in enumerate(columns) if not row & column)
                          for row in rows]
                assert binary_rank(matrix) == disjoint_rank(n, a, b)
                disjoint_cases += 1
    derivative_cases, derivative_vectors = [], 0
    for modes in range(2, 11, 2):
        for power in range(modes//2+1):
            top = sum(1 << word for word in paired_top(modes, power))
            groups = defaultdict(list)
            for deleted in range(1 << modes):
                value = derivative(top, deleted)
                if value:
                    words = support(value)
                    groups[(words[0].bit_count(), height(words[0]))].append(value)
                    derivative_vectors += 1
            actual = {key: binary_rank(rows) for key, rows in groups.items()}
            expected = graded_dimensions(modes, power)
            assert actual == expected
            derivative_cases.append(dict(modes=modes, power=power, q=sum(expected.values())))
    large = []
    for modes, power in ((16, 8), (32, 4), (64, 3), (128, 2), (512, 2), (2048, 1)):
        counts = graded_dimensions(modes, power)
        total = sum(counts.values())
        assert all(counts.get((2*power-d, 2*modes*power-h)) == value
                   for (d, h), value in counts.items())
        if power == 2:
            assert total == comb(modes, 2)+2*modes+2-(modes//2 % 2)
        if power == 1:
            assert total == modes+2
        if power == modes//2:
            assert total == 1 << modes
        large.append(dict(modes=modes, power=power, q=total,
                          occupied_height_cells=len(counts)))
    symplectic_entries, smith_orders = 0, 0
    for modes in (2, 4, 8, 16, 32):
        matrix = one_particle(modes)
        mate = lambda i: (i//modes)*modes+modes-1-i%modes
        for i, j in product(range(2*modes), repeat=2):
            value = 0
            for k in range(2*modes):
                value ^= pmul(matrix[i][k], matrix[j][mate(k)])
            assert value == ((1 << (4*modes)) if j == mate(i) else 0)
            symplectic_entries += 1
        if modes <= 8:
            orders, _ = WeightedSmithOracle([[polynomial((v, 0)) for v in row] for row in matrix]).profile(
                [0]*(2*modes), [0]*(2*modes))
            assert list(orders) == list(range(1, 4*modes, 2))
            smith_orders += len(orders)

    # Actually unequal alphabets, including a nonorthogonal Y map so that
    # exchanging its embedding and retraction changes the second cut.
    y_code = inverse_system([sum(1 << word for word in paired_top(4, 1))], 4)
    coordinate = [1 << word for word in (0, 1, 2, 4, 8)]
    literal_entries, literal_orders = 0, 0
    for modes, px, lx, py, ly in (
            (2, [1, 2], [1, 2], [1, 2, 4, 8], [1, 2, 4, 8]),
            (4, coordinate, coordinate, y_code['embedding'], y_code['retraction'])):
        result = unequal_profiles(px, lx, py, ly, modes, keep_blocks=True)
        ambient = one_particle(modes)
        for cut in result['cuts']:
            for block in cut['blocks']:
                pairs = block['pairs']
                literal = []
                for i, (x, y) in enumerate(pairs):
                    row = []
                    for j, (u, v) in enumerate(pairs):
                        value = 0
                        left = (support(lx[x]), support(ly[y] if cut['cut'] == 'zero' else py[y]))
                        right = (support(py[v] if cut['cut'] == 'zero' else ly[v]), support(px[u]))
                        for a, b, c, d in product(*left, *right):
                            value ^= determinant(ambient, a | (b << modes), c | (d << modes))
                        assert value == block['matrix'][i][j]
                        row.append(polynomial((value, 0)))
                        literal_entries += 1
                    literal.append(row)
                orders, _ = WeightedSmithOracle(literal).profile([0]*len(pairs), [0]*len(pairs))
                assert list(orders) == block['orders_t']
                literal_orders += len(orders)
    code = inverse_system([sum(1 << word for word in paired_top(4, 1))], 4)
    equal = unequal_profiles(code['embedding'], code['retraction'], code['embedding'], code['retraction'], 4)
    original = profiles(code['embedding'], code['retraction'], 4)
    for key in ('zero_orders_t', 'infinity_orders_t'):
        assert equal[key] == original[key]
    low_jet_cases = []
    family_path = ROOT/'evidence/coherent-inverse-f4-paired-family.json'
    retained = json.loads(family_path.read_text())['shards'][0]['report']['records'][0]
    for modes, power, computed in ((4, 1, original), (8, 2, retained)):
        cutoff = min(2*modes, (power+1)**2)
        expected = [0]
        for weight in range(1, cutoff, 2):
            expected += [value+weight for value in list(expected) if value+weight < cutoff]
        expected.sort()
        for key in ('zero_orders_t', 'infinity_orders_t'):
            assert [value for value in computed[key] if value < cutoff] == expected
        low_jet_cases.append(dict(modes=modes, power=power, cutoff=cutoff,
                                  orders_below_cutoff=expected))
    names = (Path(__file__).name, 'coherent_divided_powers.py', 'coherent_unequal_compounds.py',
             'check_binary_coefficient_entropy.py', 'coherent_inverse_systems.py',
             'coherent_hot_compounds.py', 'coherent_hot_compounds.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'check_coherent_hot_compounds.py', 'check_cancellation_module_automata.py',
             'weighted_smith_oracle_f4.py', 'check_joint_energy_conjugations.py')
    return dict(disjoint_rank_cases=disjoint_cases, literal_derivative_cases=derivative_cases,
                nonzero_literal_derivatives=derivative_vectors, large_closed_formula_cases=large,
                exact_symplectic_polynomial_entries=symplectic_entries,
                independent_one_particle_Smith_orders=smith_orders,
                unequal_literal_polynomial_entries=literal_entries,
                unequal_independent_Smith_orders=literal_orders, equal_alphabet_controls=2,
                low_jet_cases=low_jet_cases,
                retained_profile_sha256=sha256(family_path.read_bytes()).hexdigest(),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names},
                scope='Exact finite controls for the stated all-field dimension and symplectic formulas; no new crossing rate or factoring score.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/coherent-divided-powers-controls.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'literal_derivative_cases', 'large_closed_formula_cases')}))


if __name__ == '__main__':
    main()
