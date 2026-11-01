#!/usr/bin/env python3
"""Binary specialization and independent literal controls for F4 inverse systems."""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path

from check_coherent_hot_compounds import determinant
from check_joint_energy_conjugations import polynomial
from coherent_hot_compounds import one_particle, profiles as binary_profiles, support
from coherent_hot_compounds_f4 import axpy, basis, mul, profiles
from coherent_inverse_systems import derivative as binary_derivative, inverse_system as binary_inverse
from coherent_inverse_systems_f4 import inverse_system
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]


def quartic_catalecticant():
    """Symbolic height blocks over F2[a,b], plus all sixteen F4 fibers."""
    words = [sum(1 << i for i in pair) for pair in combinations(range(8), 2)]
    height = lambda word: sum(2*i+1 for i in range(8) if word >> i & 1)
    top = {word: (0, 0) for word in (60, 90, 102, 153, 165, 195)}
    top.update({105: (1, 0), 150: (0, 1)})
    blocks = []
    for h in sorted({height(word) for word in words}):
        rows = [word for word in words if height(word) == h]
        columns = [word for word in words if height(word) == 32-h]
        matrix = [[top.get(row | column) if not row & column else None
                   for column in columns] for row in rows]
        determinant_terms = set()
        for permutation in permutations(range(len(rows))):
            selected = [matrix[i][j] for i, j in enumerate(permutation)]
            if all(term is not None for term in selected):
                term = tuple(sum(value[k] for value in selected) for k in range(2))
                determinant_terms.symmetric_difference_update((term,))
        expected = {(0, 0), (1, 1)} if h in (8, 12, 14, 18, 20, 24) else {(0, 0)}
        assert determinant_terms == expected
        # Every variable block is I plus one explicit 2-by-2 block
        # [[1,a],[b,1]], up to simultaneous row/column permutations.
        if len(expected) == 2:
            variable_rows = [i for i, row in enumerate(matrix)
                             if any(value in ((1, 0), (0, 1)) for value in row)]
            assert len(variable_rows) == 2
            sub = [[matrix[i][j] for j in variable_rows] for i in variable_rows]
            assert sub[0][0] == sub[1][1] == (0, 0)
            assert {sub[0][1], sub[1][0]} == {(1, 0), (0, 1)}
            assert all(matrix[i][j] == ((0, 0) if i == j else None)
                       for i in range(len(rows)) for j in range(len(rows))
                       if i not in variable_rows or j not in variable_rows)
        blocks.append(dict(derivative_height=h, row_words=rows, column_words=columns,
                           entries=[[list(term) if term is not None else None for term in row]
                                    for row in matrix],
                           determinant_monomials=[list(term) for term in sorted(determinant_terms)]))
    fibers = []
    for a, b in product(range(4), repeat=2):
        form = {word: 1 for word in (60, 90, 102, 153, 165, 195)}
        form.update({word: value for word, value in ((105, a), (150, b)) if value})
        derivatives = [{word ^ deleted: value for word, value in form.items()
                        if word & deleted == deleted} for deleted in words]
        middle = len(basis(derivatives))
        code = inverse_system([form], 8)
        hilbert = [code['charges'].count(k) for k in range(5)]
        assert middle == (22 if mul(a, b) == 1 else 28)
        assert hilbert == [1, 8, middle, 8, 1]
        fibers.append(dict(a=a, b=b, hilbert=hilbert))
    return dict(symbolic_height_blocks=blocks, F4_fibers=fibers,
                determinant='(1+a*b)^6',
                all_characteristic_two_hilbert='(1,8,22,8,1) if a*b=1; (1,8,28,8,1) otherwise')


def certificate():
    specializations = []
    for modes, forms in ((4, [{6: 1, 9: 1}]), (5, [{21: 1, 14: 1}]),
                         (6, [{21: 1, 14: 1}, {9: 1}])):
        old = binary_inverse([sum(1 << word for word in form) for form in forms], modes)
        new = inverse_system(forms, modes)
        assert new['embedding'] == [{word: 1 for word in support(v)} for v in old['embedding']]
        assert new['retraction'] == [{word: 1 for word in support(v)} for v in old['retraction']]
        specializations.append(dict(modes=modes, q=len(new['embedding'])))
    forms = [{6: 1, 9: 2}, {21: 1, 14: 2}]
    reconstructions = 0
    for modes, top in zip((4, 5), forms):
        code = inverse_system([top], modes)
        for scalar in (1, 2, 3):
            changed = inverse_system([{word: mul(scalar, value) for word, value in top.items()}], modes)
            assert changed['embedding'] == code['embedding']
            assert changed['retraction'] == code['retraction']
        conjugate = inverse_system([{word: mul(value, value) for word, value in top.items()}], modes)
        assert conjugate['embedding'] == [{word: mul(value, value) for word, value in v.items()}
                                         for v in code['embedding']]
        assert conjugate['retraction'] == code['retraction']
        planes = [sum(((value >> p) & 1) << word for word, value in top.items()) for p in range(2)]
        for deleted in range(1 << modes):
            derivative_planes = [binary_derivative(value, deleted) for value in planes]
            expected = {word: sum(((value >> word) & 1) << p for p, value in enumerate(derivative_planes))
                        for word in range(1 << modes)}
            expected = {word: value for word, value in expected.items() if value}
            actual = {}
            for pivot, vector in zip(code['pivot_words'], code['embedding']):
                actual = axpy(actual, expected.get(pivot, 0), vector)
            assert actual == expected
            reconstructions += 1
    first = {6: 1, 9: 2}
    second = {6: 2, 9: 1}
    combined = axpy(first, 1, second)
    original = inverse_system([first, second], 4)
    presented = inverse_system([second, combined], 4)
    assert original['embedding'] == presented['embedding']
    assert original['retraction'] == presented['retraction']

    code = inverse_system([{6: 1, 9: 2}], 4)
    result = profiles(code['embedding'], code['retraction'], 4, keep_blocks=True)
    p, l, charges = code['embedding'], code['retraction'], code['charges']
    ambient = one_particle(4)

    @lru_cache(None)
    def literal_minor(row, column):
        return determinant(ambient, row, column)

    entries, nonbinary, channels = 0, 0, 0
    for cut in result['cuts']:
        levels = sorted({a+b for a in charges for b in charges})
        for charge, block in zip(levels, cut['blocks']):
            pairs = [(x, y) for x in range(len(p)) for y in range(len(p)) if charges[x]+charges[y] == charge]
            literal = []
            for x, y in pairs:
                row = []
                for u, v in pairs:
                    value = [0, 0]
                    left = (l[x], l[y] if cut['cut'] == 'zero' else p[y])
                    right = (p[v] if cut['cut'] == 'zero' else l[v], p[u])
                    for (a, ca), (b, cb), (c, cc), (d, cd) in product(
                            left[0].items(), left[1].items(), right[0].items(), right[1].items()):
                        scalar = mul(mul(ca, cb), mul(cc, cd))
                        coefficient = literal_minor(a | (b << 4), c | (d << 4))
                        if scalar & 1:
                            value[0] ^= coefficient
                        if scalar & 2:
                            value[1] ^= coefficient
                    value = tuple(value)
                    assert value == block['matrix'][len(literal)][len(row)]
                    row.append(polynomial(value))
                    entries += 1
                    nonbinary += bool(value[1])
                literal.append(row)
            orders, _ = WeightedSmithOracle(literal).profile([0]*len(pairs), [0]*len(pairs))
            assert list(orders) == block['orders_t']
            channels += len(orders)
    assert entries == 908 and channels == 72 and nonbinary
    binary_code = binary_inverse([(1 << 6) ^ (1 << 9)], 4)
    old = binary_profiles(binary_code['embedding'], binary_code['retraction'], 4)
    new = inverse_system([{6: 1, 9: 1}], 4)
    binary_result = profiles(new['embedding'], new['retraction'], 4)
    for name in ('zero_orders_t', 'infinity_orders_t'):
        assert old[name] == binary_result[name]
    names = ('check_coherent_inverse_systems_f4.py', 'coherent_inverse_systems_f4.py',
             'coherent_inverse_systems.py', 'coherent_hot_compounds_f4.py',
             'coherent_hot_compounds_f4.cpp', 'coherent_hot_compounds.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'check_coherent_hot_compounds.py', 'weighted_smith_oracle_f4.py',
             'check_joint_energy_conjugations.py')
    return dict(binary_specializations=specializations, nonbinary_derivative_reconstructions=reconstructions,
                literal_complete_F4_entries=entries, nonbinary_entries=nonbinary,
                independently_replayed_fraction_free_orders=channels,
                binary_cut_profile_controls=2, presentation_and_Frobenius_controls=True,
                quartic_catalecticant=quartic_catalecticant(),
                six_letter_zero_orders=result['zero_orders_t'],
                six_letter_infinity_orders=result['infinity_orders_t'],
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names},
                scope='Exact constructor, closure, endpoint-entry and complete small-profile controls; no numerical rate certificate or field-wide optimum.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/coherent-inverse-f4-controls.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items() if k != 'source_sha256'}))


if __name__ == '__main__':
    main()
