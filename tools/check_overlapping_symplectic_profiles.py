#!/usr/bin/env python3
"""Check actual overlapping-form code jets against both Ward rank bounds."""

import argparse
import hashlib
import json
from pathlib import Path

from check_double_endpoint_retracts import bit_indices, optimize
from check_modular_crossing import rank_two
from explore_overlapping_symplectic_codes import overlapping_code
from explore_symplectic_boundary_codes import coefficient_basis
from ward_rectangular_bounds import ward_score


def submasks(value):
    answer = []
    current = value
    while True:
        answer.append(current)
        if not current:
            return answer
        current = (current - 1) & value


def audit(data, row, column):
    jets = {}
    for i, polynomial in enumerate(data['basis_scalar_polynomials_p_q_r']):
        for p, q, _ in polynomial:
            for a in submasks(p):
                for b in submasks(q):
                    jets[a, b] = jets.get((a, b), 0) ^ (1 << i)
    code = data['code']
    bound = ward_score(code['pairs'], code['length'], code['form_count'], code['labels'])
    ranks = {(a, b): (r, c) for a, b, r, c in bound['profile']}
    exact = []
    for (a, b), mask in sorted(jets.items()):
        if not mask:
            continue
        actual = []
        for cut in (row, column):
            blocks = [[0] * size for size in data['block_sizes']]
            for i in bit_indices(mask):
                for target, source in zip(blocks, cut[i]):
                    for j, value in enumerate(source):
                        target[j] ^= value
            actual.append(sum(rank_two(block) for block in blocks))
        assert b >= code['form_count']
        assert all(x <= y for x, y in zip(actual, ranks.get((a, b), (0, 0)))), (a, b, actual)
        exact.append([a, b, *actual])
    finite = data.get('F2_optimization') or optimize(data, row, column)
    return {
        'form_column_types': code['form_column_types'],
        'pairs': code['pairs'], 'length': code['length'], 'forms': code['form_count'],
        'labels': code['labels'], 'coefficient_dimension': data['basis_dimension'],
        'full_product_image_dimension': code['full_product_image_homogeneous_dimension'],
        'zero_syndrome_dimension': code['zero_syndrome_dimension'],
        'maximum_input_terms': code['maximum_input_terms'],
        'common_q_plus_r_power': data['common_q_plus_r_power'],
        'common_p_plus_r_power': data['common_p_plus_r_power'],
        'p_endpoint_zero': data['p_endpoint_zero'], 'q_endpoint_ranks': data['q_endpoint_ranks'],
        'F2_optimization': finite,
        'literal_centered_jet_ranks': exact,
        'ward_rank_bound': bound,
        'all_literal_jets_satisfy_both_Ward_bounds': True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', action='append', type=Path, default=[])
    parser.add_argument('--large', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = []
    if args.input:
        for path in args.input:
            data = json.loads(path.read_text())
            row, column = ([[[int(v, 16) for v in block] for block in matrix]
                            for matrix in data[key]]
                           for key in ('basis_row_blocks', 'basis_column_blocks'))
            cases.append(audit(data, row, column))
    else:
        for columns, length in [([1, 1, 2, 2, 3, 3], 6)] + (
                [([1, 1, 1, 2, 2, 3, 3], 7)] if args.large else []):
            code = overlapping_code(columns, length)
            data, row, column = coefficient_basis(code)
            data['code'] = code
            cases.append(audit(data, row, column))
    sources = [Path(__file__), Path(__file__).with_name('explore_overlapping_symplectic_codes.py'),
               Path(__file__).with_name('explore_symplectic_boundary_codes.py'),
               Path(__file__).with_name('ward_rectangular_bounds.py')]
    result = {'scope': 'The two recorded overlapping-form codes; complete binary coefficient spaces.',
              'cases': cases,
              'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(json.dumps({**result, 'cases': [{k: v for k, v in case.items()
                                       if k not in ('F2_optimization', 'ward_rank_bound', 'literal_centered_jet_ranks')}
                                      for case in cases]}))


if __name__ == '__main__':
    main()
