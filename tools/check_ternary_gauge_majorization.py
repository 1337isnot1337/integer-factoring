#!/usr/bin/env python3
"""Exact balanced-minor certificates for uniform-type ternary gauge profiles.

The default checker uses only exact rational and finite-field arithmetic.
The optional --propose mode uses SciPy LPs to find small rational frames;
every proposed frame is then checked independently by this same verifier.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_all_field_pascal_crossings import determinant, trim
from check_multipoint_hermite_crossings import local_polynomial
from check_ternary_basis_mixed import transform_seed, ramified_profile


BASES = {'inverse_evaluation': [[1, 0, 0], [0, 2, 1], [2, 2, 2]],
         'inverse_pascal': [[1, 0, 0], [2, 1, 0], [1, 1, 1]]}
ROW_PROFILE = tuple(Fraction(value, 2) for value in (0, 1, 1, 2, 2, 2, 3, 3, 4))
SLOPES = ((0, 0, Fraction(1, 2), Fraction(1, 2)),
          (-Fraction(1, 2), -Fraction(1, 2), 0, 0), (-1, -1, 1, 1))
LOCATIONS = ((2, 0), (None, 0), (0, 1))


def encode(value):
    value = Fraction(value)
    return [value.numerator, value.denominator]


def labels(row, column, cut):
    x, second = divmod(row, 3)
    third, v = divmod(column, 3)
    return (x, second, third, v) if cut == 0 else (x, third, second, v)


def minor_counts(rows, columns, cut):
    counts = [0] * 8
    row_legs, column_legs = ((0, 1), (2, 3)) if cut == 0 else ((0, 2), (1, 3))
    for indices, legs in ((rows, row_legs), (columns, column_legs)):
        for index in indices:
            for leg, symbol in zip(legs, divmod(index, 3)):
                if symbol:
                    counts[2 * leg + symbol - 1] += 1
    return counts


def minor_polynomial(seed, rows, columns, cut):
    return determinant([[trim(seed.get(labels(row, column, cut), ()), 3) for column in columns]
                        for row in rows], 3)


def propose():
    from itertools import combinations
    import numpy as np
    from scipy.optimize import linprog
    from search_ternary_physical_bases import change_basis, oracles_for
    cases = []
    for name, basis in BASES.items():
        seed = transform_seed(basis, basis)
        oracles, _ = oracles_for(change_basis(np.array(basis), np.array(basis)))
        places = []
        for index, ((point, cut), oracle) in enumerate(zip(LOCATIONS, oracles)):
            reference = ROW_PROFILE if index < 2 else tuple(2 * value for value in ROW_PROFILE)
            frames = []
            for size in range(1, 10):
                subsets = list(combinations(range(9), size))
                masks = np.array([sum(1 << i for i in subset) for subset in subsets])
                row_counts = np.array([oracle.row_counts[list(subset)].sum(axis=0) for subset in subsets])
                col_counts = np.array([oracle.column_counts[list(subset)].sum(axis=0) for subset in subsets])
                table = oracle.valuations[masks[:, None], masks[None, :]]
                ri, ci = np.nonzero(table < 10000)
                features = 3 * (row_counts[ri] + col_counts[ci]) - size
                values = table[ri, ci]
                answer = linprog(np.r_[np.zeros(8), -1.0],
                                 A_ub=np.c_[-features, np.ones(len(values))], b_ub=3 * values,
                                 bounds=[(None, None)] * 9, method='highs')
                if not answer.success:
                    raise ArithmeticError(answer.message)
                selected = np.flatnonzero(-answer.ineqlin.marginals > 1e-8)
                minors = []
                for position in selected:
                    weight = Fraction(float(-answer.ineqlin.marginals[position])).limit_denominator(100000)
                    rows, columns = subsets[int(ri[position])], subsets[int(ci[position])]
                    polynomial = minor_polynomial(seed, rows, columns, cut)
                    minors.append({'rows': list(rows), 'columns': list(columns), 'weight': encode(weight),
                                   'valuation': int(values[position]), 'polynomial': list(polynomial)})
                frames.append({'size': size, 'prefix_bound': encode(sum(reference[:size])), 'minors': minors})
            places.append({'point': point, 'cut': cut, 'effective_profile': list(map(encode, reference)),
                           'frames': frames})
        cases.append({'name': name, 'basis_x': basis, 'basis_y': basis, 'places': places})
    return {'schema_version': 1,
            'scope': 'universal rational physical gauges on two specified bases and uniform paired types only',
            'cases': cases}


def check(certificate):
    totals = {'bases': 0, 'balanced_frames': 0, 'independent_minor_determinants': 0,
              'exact_attaining_profiles': 0, 'endpoint_entries': 0}
    details = []
    assert {case['name'] for case in certificate['cases']} == set(BASES)
    for case in certificate['cases']:
        assert case['basis_x'] == case['basis_y'] == BASES[case['name']]
        seed = transform_seed(case['basis_x'], case['basis_y'])
        for x, y, u, v in product(range(3), repeat=4):
            assert sum(seed.get((x, y, u, v), ())) % 3 == int(x == v and y == u)
            totals['endpoint_entries'] += 1
        assert len(case['places']) == 3
        for index, location in enumerate(case['places']):
            point, cut = LOCATIONS[index]
            assert location['point'] == point and location['cut'] == cut
            reference = ROW_PROFILE if index < 2 else tuple(2 * value for value in ROW_PROFILE)
            assert location['effective_profile'] == list(map(encode, reference))
            actual, _ = ramified_profile(seed, point, cut, SLOPES[index])
            mean = sum(Fraction(slope) for slope in SLOPES[index])
            restored = [Fraction(order, 2) - mean for order in actual]
            assert restored == list(reference)
            totals['exact_attaining_profiles'] += 1
            assert len(location['frames']) == 9
            for size, frame in enumerate(location['frames'], 1):
                assert frame['size'] == size
                target = Fraction(*frame['prefix_bound'])
                assert target == sum(reference[:size])
                total_weight = Fraction()
                average_counts = [Fraction()] * 8
                average_valuation = Fraction()
                for minor in frame['minors']:
                    rows, columns = minor['rows'], minor['columns']
                    assert len(rows) == len(columns) == size
                    assert rows == sorted(set(rows)) and columns == sorted(set(columns))
                    assert all(0 <= value < 9 for value in rows + columns)
                    weight = Fraction(*minor['weight'])
                    assert weight > 0
                    polynomial = minor_polynomial(seed, rows, columns, cut)
                    assert polynomial and list(polynomial) == minor['polynomial']
                    if size == 9:
                        expected = [1] + [0] * 8 + [1] if cut == 0 else [0] * 18 + [2]
                        assert list(polynomial) == expected
                    local = local_polynomial(polynomial, point, 2 * size, 3)
                    valuation = next(i for i, coefficient in enumerate(local) if coefficient)
                    assert valuation == minor['valuation']
                    counts = minor_counts(rows, columns, cut)
                    total_weight += weight
                    average_valuation += weight * valuation
                    average_counts = [a + weight * count for a, count in zip(average_counts, counts)]
                    totals['independent_minor_determinants'] += 1
                assert total_weight == 1
                assert average_counts == [Fraction(size, 3)] * 8
                assert average_valuation == target
                totals['balanced_frames'] += 1
            details.append({'basis': case['name'], 'point': point, 'cut': cut,
                            'sharp_effective_profile': list(map(encode, reference))})
        totals['bases'] += 1
    return {'status': 'pass', 'counts': totals, 'sharp_profiles': details,
            'scope': certificate['scope'], 'no_global_basis_or_type_optimality_claim': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('evidence/ternary-uniform-gauge-majorization.json'))
    parser.add_argument('--propose', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.propose:
        certificate = propose()
        result = check(certificate)
        args.input.parent.mkdir(parents=True, exist_ok=True)
        args.input.write_text(json.dumps(certificate, indent=2, sort_keys=True) + '\n')
    else:
        certificate = json.loads(args.input.read_text())
        result = check(certificate)
    directory = Path(__file__).parent
    sources = (Path(__file__).name, 'check_all_field_pascal_crossings.py',
               'check_ternary_basis_mixed.py', 'check_multipoint_hermite_crossings.py')
    result['source_sha256'] = {name: sha256((directory / name).read_bytes()).hexdigest() for name in sources}
    result['certificate_sha256'] = sha256(args.input.read_bytes()).hexdigest()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result['counts']))


if __name__ == '__main__':
    main()
