#!/usr/bin/env python3
"""Exact pairing, basis-independence and two-cut Gram-equivalence controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from coherent_hot_compounds_f4 import axpy, mul, profiles
from coherent_inverse_systems_f4 import inverse_system
from coherent_orthogonal_f4 import orthogonal_profiles, paired_retraction

ROOT = Path(__file__).resolve().parents[1]


def certificate():
    cases = []
    for modes, forms, dual_forms in (
            (4, [{6: 1, 9: 2}], [{6: 1, 9: 1}]),
            (5, [{21: 1, 14: 2}], [{21: 1, 14: 1}]),
            (6, [{12: 1, 18: 1}, {18: 1, 33: 2}],
             [{12: 1, 18: 2}, {18: 1, 33: 2}])):
        embedding = inverse_system(forms, modes)['embedding']
        dual = inverse_system(dual_forms, modes)['embedding']
        maps = paired_retraction(embedding, dual, modes)
        g, inverse = maps['gram'], maps['inverse_gram']
        size = len(g)
        for left, right in ((g, inverse), (inverse, g)):
            for i in range(size):
                for j in range(size):
                    value = 0
                    for k in range(size):
                        value ^= mul(left[i][k], right[k][j])
                    assert value == int(i == j)
        changed = [{word: mul(2, value) for word, value in row.items()} for row in dual]
        # A non-diagonal basis change inside a homogeneous cell.
        same_grade = next(((i, j) for i in range(size) for j in range(i)
                           if (maps['heights'][i], maps['charges'][i]) ==
                           (maps['heights'][j], maps['charges'][j])), None)
        if same_grade:
            i, j = same_grade
            changed[i] = axpy(changed[i], 3, changed[j])
        assert paired_retraction(embedding, changed, modes)['retraction'] == maps['retraction']
        self_maps = paired_retraction(embedding, embedding, modes)
        ordinary = profiles(embedding, self_maps['retraction'], modes)
        derived = orthogonal_profiles(embedding, modes)
        assert ordinary['zero_orders_t'] == derived['zero_orders_t']
        assert ordinary['infinity_orders_t'] == derived['infinity_orders_t']
        assert ordinary['zero_orders_t'] == ordinary['infinity_orders_t']
        cases.append(dict(modes=modes, q=size, non_diagonal_basis_change=bool(same_grade),
                          non_diagonal_self_Gram=any(self_maps['gram'][i][j]
                                                   for i in range(size) for j in range(size) if i != j),
                          zero_orders_t=ordinary['zero_orders_t'],
                          infinity_orders_t=ordinary['infinity_orders_t']))
    assert any(case['non_diagonal_self_Gram'] for case in cases)
    p = inverse_system([{6: 1, 9: 2}], 4)['embedding']
    singular_dual = inverse_system([{6: 1, 9: 3}], 4)['embedding']
    try:
        paired_retraction(p, singular_dual, 4)
    except ValueError as error:
        assert 'singular' in str(error)
    else:
        raise AssertionError('An orthogonal pair was accepted as nondegenerate.')
    names = ('check_coherent_gram_f4.py', 'coherent_orthogonal_f4.py',
             'coherent_inverse_systems_f4.py', 'coherent_hot_compounds_f4.py',
             'coherent_hot_compounds_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp', 'coherent_hot_compounds.py')
    return dict(cases=cases, exact_Gram_inverse_products=6,
                independent_dual_basis_changes=3, singular_pair_rejected=True,
                computed_both_cut_orders=sum(2*case['q']**2 for case in cases),
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in names},
                scope='Finite exact controls for ordinary bilinear pairings and the constant-unit two-cut theorem; no field-wide pressure optimum.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/coherent-gram-f4-controls.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({key: value for key, value in result.items() if key not in ('source_sha256', 'cases')}))


if __name__ == '__main__':
    main()
