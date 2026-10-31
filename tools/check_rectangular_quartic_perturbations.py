#!/usr/bin/env python3
"""Exact polynomial crossing and reduced-profile controls for Laurent quartics."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_rectangular_f4_repair import build
from check_rectangular_defect_core_basis import data_for
from check_rectangular_defect_coordinates import original_cut, coordinate_weights
from rectangular_defect_coordinates import coordinate_matrix
from native_matrix_smith_f4 import MatrixOracle


def profile_control(seed, terms, place, gauges):
    data = data_for(place, 4)
    matrix, shifts, metadata = coordinate_matrix(data, quartic_terms=terms)
    rw, cw = coordinate_weights(metadata, gauges)
    expected = MatrixOracle(matrix, len(matrix), shifts).profile(rw, cw)['orders']
    if place == 'infinity':
        spectator = -3+gauges[0]+gauges[1]
        expected = sorted(expected+[value+spectator for value in expected])
    matrix, rw, cw = original_cut(seed, place, gauges)
    actual = MatrixOracle(matrix, len(matrix)).profile(rw, cw, infinity=place == 'infinity')['orders']
    # The exact seed clears negative Laurent exponents by a global power of w.
    restore = (seed['denominator_degree']+seed['global_w_power'] if place == 'infinity'
               else -seed['global_w_power'])
    actual = sorted(value+restore for value in actual)
    error = max(abs(a-b) for a, b in zip(actual, sorted(expected)))
    assert len(actual) == len(expected) and error < 1e-9, (terms, place, error)
    return {'place': place, 'quartic_terms': terms, 'global_w_power_cleared_in_seed': seed['global_w_power'],
            'lambda': gauges[0], 'mu': gauges[1], 'v': gauges[2], 'y': gauges[3], 'u': gauges[4],
            'orders_checked': len(actual), 'maximum_error': error,
            'orders_sha256': sha256(json.dumps(actual, separators=(',', ':')).encode()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    terms_list = [[(1, exponent)] for exponent in range(-4, 5)]
    terms_list += [[(1, 0), (2, exponent), (2, exponent+1)] for exponent in (-2, 0, 3)]
    exact = []
    for terms in terms_list:
        seed = build((0, 1, 2), 4, 1, 2, independent_exterior=True, quartic_terms=terms)
        exact.append({key: seed[key] for key in ('quartic_terms', 'global_w_power', 'degree',
                                               'denominator_degree', 'ramified_endpoint_checks',
                                               'two_copy_type_checks', 'independent_exterior_checks',
                                               'coefficient_sha256')})
    gauges = (0.75, -1.25, [0, 0.5, -1, 1.5], [0, -0.75, 1.25], [0, 0.25, -0.5])
    profiles = []
    for terms in ([(1, -3)], [(1, 2)], [(1, 0), (2, -1), (2, 0)]):
        seed = build(tuple(range(4)), 4, 1, 2, check_types=False, quartic_terms=terms)
        for place in ('infinity', 'zero'):
            profiles.append(profile_control(seed, terms, place, gauges))
    result = {'status': 'exact Laurent-quartic crossing controls and independent dense profile replay',
              'polynomial_crossing_cases': exact, 'dense_profile_cases': profiles}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
