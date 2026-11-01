#!/usr/bin/env python3
"""Exact module, Laurent-space and independent profile controls."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_joint_energy_conjugations import polynomial
from check_joint_hot_retraction_module import product as matrix_product, at_two, rank_f4
from joint_hot_retraction_module import degree, dependency, polynomial_right_inverse
from joint_six_to_five_retractions import H, HEIGHTS, MU_ZERO, MU_INFINITY, crossing, module, profiles, row_space, source
from joint_small_rational_retractions import saturated_basis as small_basis, source_rows as small_source, row_space as small_space
from polynomial_retraction_spaces import bounded_laurent_space, row_reduced_basis
from search_f4_recentered_crossings import five_letter_source
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]


def span(rows):
    pivots = {}
    for v in rows:
        while v:
            k = v.bit_length()-1
            if k in pivots:
                v ^= pivots[k]
            else:
                pivots[k] = v
                break
    return pivots


def packed_rows(rows, width):
    return [sum(p << (width*j) for j, p in enumerate(row)) for row in rows]


def compare_spaces(a, b):
    aa, bb = packed_rows(a, 17), packed_rows(b, 17)
    assert len(span(aa)) == len(aa) == len(bb) == len(span(bb))
    assert len(span(aa+bb)) == len(aa)


def python_profiles(code):
    clearing = max(code['row_clearing'])
    matrix = [[polynomial((a << (clearing-c), b << (clearing-c))) for a, b in row]
              for row, c in zip(code['matrix'], code['row_clearing'])]
    opposite = [[0]*25 for _ in range(25)]
    for x, y, u, v in product(range(5), repeat=4):
        opposite[5*x+u][5*v+y] = matrix[5*x+y][5*v+u]
    zero, _ = WeightedSmithOracle(matrix).profile(H, [0]*25)
    infinity, _ = WeightedSmithOracle(opposite, infinity=True).profile(
        [y for _ in HEIGHTS for y in HEIGHTS], [x for x in HEIGHTS for _ in HEIGHTS])
    return [v-clearing for v in zero], [v+clearing for v in infinity]


def certificate():
    # Regression against the separately implemented earlier polynomial spaces.
    small, _ = row_reduced_basis(small_source())
    assert small == small_basis()[0]
    comparisons = 0
    for k, h in product(range(7), range(9)):
        result = bounded_laurent_space(small, 8+k-h,
                                       [8+a+b-k for a in (0, 1, 3, 4) for b in (0, 1, 3, 4)])
        assert result['clearing'] == 0
        compare_spaces(result['rows'], small_space(h, k))
        comparisons += 1
    data = source()
    assert len(data['selected']) == 25 and len(data['excluded']) == 11
    assert all(p.bit_count() % 2 == int(i == j) for i, row in enumerate(data['matrix']) for j, p in enumerate(row))
    old, _, _, _ = five_letter_source()
    permutation = (0, 3, 2, 4, 1)
    for x, y, u, v in product(range(5), repeat=4):
        index = permutation[x]+5*permutation[y]+25*permutation[u]+125*permutation[v]
        assert old[index][1] == 0
        assert data['retained'][5*x+y][5*v+u] == old[index][0] << 4
    basis, steps = module()
    assert [degree(r) for r in basis] == [0, 7, 4, 7, 4, 9, 6, 9, 16]
    assert not dependency([sum((p & 1) << j for j, p in enumerate(row)) for row in basis])
    assert not dependency([sum(((p >> degree(row)) & 1) << j for j, p in enumerate(row)) for row in basis])
    inverse, bound = polynomial_right_inverse(basis)
    assert bound == 5 and matrix_product(basis, inverse) == [[int(i == j) for j in range(9)] for i in range(9)]
    conversion = matrix_product(data['source_rows'], inverse)
    assert matrix_product(conversion, basis) == data['source_rows']
    evaluated = at_two(conversion)
    assert rank_f4(evaluated) == 9
    pivot_rows = []
    for i, row in enumerate(evaluated):
        if rank_f4([evaluated[j] for j in pivot_rows]+[row]) > len(pivot_rows):
            pivot_rows.append(i)
    assert len(pivot_rows) == 9
    for i, j in product(range(25), repeat=2):
        p = data['retained'][i][j]
        if p:
            assert (p & -p).bit_length()-1 >= MU_ZERO-H[i]
            assert p.bit_length()-1 <= H[j]-MU_INFINITY
    dims, space_checks = {}, 0
    for k in range(7):
        dims[str(k)] = []
        for h in sorted(set(H)):
            result = row_space(h, k)
            assert not dependency(packed_rows(result['rows'], 64))
            assert all(p.bit_count() % 2 == 0 for row in result['rows'] for p in row)
            dims[str(k)].append([h, result['dimension'], result['clearing']])
            space_checks += 1
    cases = [[], [(0, 0, 1)], [(0, row_space(0)['dimension']-1, 2)],
             [(24, 0, 2)], [(24, row_space(H[24])['dimension']-1, 3)],
             [(24, 0, 1), (3, row_space(H[3])['dimension']//2, 3)],
             [(7, row_space(H[7])['dimension']//2, 2), (13, 0, 1)],
             [(0, 0, 3), (24, 0, 3), (12, 0, 2)]]
    witnesses = []
    for assignments in cases:
        code = crossing(assignments)
        native = profiles(code)
        zero, infinity = python_profiles(code)
        assert zero == native['raw_zero_orders'] and infinity == native['raw_infinity_orders']
        assert native['normalized_zero_orders'] == [v-MU_ZERO for v in zero]
        assert native['normalized_infinity_orders'] == [v-MU_INFINITY for v in infinity]
        assert native['declared_deficit'] == 6
        witnesses.append(dict(assignments=[list(a) for a in assignments], **native,
                              actual_deficit=-min(zero)-min(infinity)))
    names = (Path(__file__).name, 'joint_six_to_five_retractions.py', 'polynomial_retraction_spaces.py',
             'joint_hot_retraction_module.py', 'check_joint_hot_retraction_module.py',
             'joint_small_rational_retractions.py', 'restrict_implicit_reflection.py',
             'reflection_dual_rail_gaussian.py', 'search_f4_recentered_crossings.py',
             'weighted_smith_oracle_f4.py', 'native_matrix_smith_f4.py',
             'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    return dict(scope='Constructive rational six-to-five Laurent output space at the stated new paired energies and deficit; no favorable profile or improved factoring exponent claimed.',
                source_letters=list(data['letters']), source_common_w_shift=data['common_w_shift'],
                source_endpoint_entries=1296, retained_specialization_entries=625,
                source_row_count=11, module_rank=9, basis=basis,
                polynomial_right_inverse=inverse, right_inverse_degree=bound,
                source_to_basis=conversion, independent_source_rows_at_omega=pivot_rows,
                saturation=steps, polynomial_space_regressions=comparisons,
                paired_energies=[2*h for h in HEIGHTS], declared_deficit=6,
                mu_zero=MU_ZERO, mu_infinity=MU_INFINITY, bounded_spaces=space_checks,
                bounded_dimensions=dims, literal_profile_cases=len(cases),
                independent_complete_orders=50*len(cases), endpoint_entries=625*len(cases),
                witnesses=witnesses,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in names})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = certificate()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads((ROOT/'evidence/joint-six-to-five-retractions.json').read_text())
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ('basis', 'polynomial_right_inverse', 'source_to_basis', 'saturation',
                       'bounded_dimensions', 'witnesses', 'source_sha256')}))


if __name__ == '__main__':
    main()
