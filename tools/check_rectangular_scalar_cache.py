#!/usr/bin/env python3
"""Exact two-base-matrix synthesis of Laurent quartic perturbations."""
import argparse
import json
from pathlib import Path

from check_rectangular_defect_core_basis import data_for
from check_rectangular_coordinate_transforms import compare_exact, truncate_fixture
from rectangular_defect_coordinates import coordinate_matrix
from rectangular_coordinate_transforms import linear_combination_coordinate_matrices


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records = []
    for place in ('infinity', 'zero'):
        data = data_for(place, 4)
        one = coordinate_matrix(data, quartic_terms=[(1, 0)])
        w = coordinate_matrix(data, quartic_terms=[(1, 1)])
        sign = -1 if place == 'infinity' else 1
        for P, Q in (([(1, 0), (2, -1), (2, 0)], [(2, -1)]),
                     ([(1, 3)], [(1, 0), (1, 1), (1, 2)])):
            local_Q = [(coefficient, sign*exponent) for coefficient, exponent in Q]
            synthesized = linear_combination_coordinate_matrices([
                (one, [(1, 0)]+local_Q), (w, local_Q)])
            native = coordinate_matrix(data, quartic_terms=P)
            assert synthesized[2]['exact_polynomial_matrix'] and native[2]['exact_polynomial_matrix']
            count = compare_exact(synthesized, native)
            shortened = linear_combination_coordinate_matrices([
                (truncate_fixture(one), [(1, 0)]+local_Q), (w, local_Q)])
            cutoff = shortened[2]['known_below']
            compared = compare_exact(shortened, synthesized, below=cutoff)
            records.append({'place': place, 'M': 4, 'quartic_P': P, 'quartic_Q': Q,
                            'exact_rational_entries_checked': count,
                            'truncated_entries_checked': compared, 'truncated_known_below': cutoff,
                            'combination_seconds': synthesized[2]['scalar_combination_seconds']})
            print(json.dumps({'place': place, 'P': P, 'entries': count}), flush=True)
    result = {'status': 'exact affine quartic cache and conservative cutoff controls', 'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
