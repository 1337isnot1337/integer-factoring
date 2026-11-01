#!/usr/bin/env python3
"""Exact controls for the full physical-letter fused place oracle."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from fused_weighted_place_profiles import FusedPlaceOracle
from check_fused_closed_profiles import from_canonical, quotient_rank
from check_fused_multipoint_profiles import local_polynomial, toeplitz_rank
from f16_identity_operator_profiles import multiply


def weighted_rank_controls(oracle, cut, point, heights):
    orders = oracle.profile(cut, point, heights)
    q = oracle.alphabet
    all_heights = [sum(heights[leg * q + key[leg]] for leg in range(4))
                   for block in oracle.blocks[cut] for row in block['labels'] for key in row]
    clearing = max(0, -min(all_heights))
    uniformizer = {1: 1} if point is None else from_canonical(point)
    maximum = max(all_heights) + clearing
    powers = [{0: 1}]
    for _ in range(maximum):
        powers.append(multiply(powers[-1], uniformizer))
    weighted = []
    for block in oracle.blocks[cut]:
        matrix = []
        for i, row in enumerate(block['matrix']):
            out = []
            for j, value in enumerate(row):
                key = block['labels'][i][j]
                exponent = sum(heights[leg * q + key[leg]] for leg in range(4)) + clearing
                local = local_polynomial(value, None, oracle.degree) if point is None else value
                out.append(multiply(local, powers[exponent]))
            matrix.append(out)
        weighted.append(matrix)
    count = 0
    for length in range(1, max(orders) + clearing + 2):
        actual = sum(toeplitz_rank(matrix, length) if point is None else
                     quotient_rank(matrix, uniformizer, length) for matrix in weighted)
        residue_degree = 1 if point is None else len(point) - 1
        expected = residue_degree * sum(max(0, length - order - clearing) for order in orders)
        assert actual == expected
        count += 1
    return count


def check_instance(m, theta):
    oracle = FusedPlaceOracle(m, theta)
    q = oracle.alphabet
    zero = (0,) * (4 * q)
    expected = {tuple(row['factor']): row['cuts'] for row in oracle.record['closed_points']}
    expected[None] = next(row['cuts'] for row in oracle.record['profiles'] if row['point'] is None)
    gauges = [zero,
              oracle.linear_heights((1, -1, 0, 0)),
              oracle.linear_heights((0, 1, -1, 0)),
              oracle.linear_heights((-1, 0, 1, 0)),
              tuple(((leg + a) ** 2 % 3) - 1 for leg in range(4) for a in range(q)),
              tuple(Fraction(((leg + 2 * a) % 5) - 2, 2) for leg in range(4) for a in range(q))]
    x_type = tuple(Fraction(a + 1, q * (q + 1) // 2) for a in range(q))
    y_type = tuple(reversed(x_type))
    exact_profiles = float_controls = shift_controls = ramification_controls = charge_controls = 0
    subgradient_controls = 0
    for point in oracle.points:
        for cut in (0, 1):
            assert oracle.profile(cut, point, zero) == tuple(expected[point][cut])
            charge = oracle.linear_heights((1, 1, -1, -1))
            assert oracle.profile(cut, point, charge) == oracle.profile(cut, point, zero)
            assert oracle.type_shift(charge, x_type, y_type) == 0
            charge_controls += 1
            for heights in gauges:
                exact = oracle.profile(cut, point, heights)
                floating = oracle.profile_float(cut, point, [float(value) for value in heights])
                assert all(abs(float(a) - b) < 1e-10 for a, b in zip(exact, floating))
                exact_profiles += 1
                float_controls += 1
                if heights == zero or heights == gauges[4]:
                    current = np.asarray([float(value) for value in heights])
                    value, gradient, derivative_scale, returned = oracle.log_profile_with_gradient(
                        cut, point, current, scale=1.0)
                    assert np.max(np.abs(returned - floating)) < 1e-10
                    for step in (-0.37, 0.29):
                        direction = np.asarray([((3 * i + 1) % 7 - 3) / 5
                                                for i in range(4 * q)])
                        changed = current + step * direction
                        changed_scale = 1.0 + step / 3
                        changed_value = oracle.log_profile_with_gradient(
                            cut, point, changed, changed_scale)[0]
                        supporting = value + step * float(gradient @ direction) + step / 3 * derivative_scale
                        assert changed_value >= supporting - 1e-9
                        subgradient_controls += 1
                doubled = oracle.profile(cut, point, tuple(2 * value for value in heights), scale=2)
                assert doubled == tuple(2 * value for value in exact)
                ramification_controls += 1
                translated = tuple(value + (2 if index < q else 0) for index, value in enumerate(heights))
                shifted = oracle.profile(cut, point, translated)
                assert shifted == tuple(value + 2 for value in exact)
                assert oracle.type_shift(translated, x_type, y_type) - oracle.type_shift(heights, x_type, y_type) == 2
                shift_controls += 1
    rank_controls = 0
    if q <= 4:
        selected = [oracle.points[0], None]
        nonlinear = next((point for point in oracle.points if point is not None and len(point) > 2), None)
        if nonlinear is not None:
            selected.append(nonlinear)
        for point in selected:
            for cut in (0, 1):
                for leg, sign in ((0, 1), (3, -1)):
                    heights = tuple(sign * int(index == leg * q + 1) for index in range(4 * q))
                    rank_controls += weighted_rank_controls(oracle, cut, point, heights)
    return {'m': m, 'theta': theta, 'alphabet': q,
            'closed_places_including_infinity': len(oracle.points),
            'exact_weighted_profiles': exact_profiles,
            'floating_api_exact_controls': float_controls,
            'constant_shift_controls': shift_controls,
            'rational_ramification_controls': ramification_controls,
            'conserved_charge_controls': charge_controls,
            'active_minor_subgradient_controls': subgradient_controls,
            'independent_weighted_quotient_ranks': rank_controls}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [check_instance(m, theta) for m, theta in ((2, 2), (3, 4), (4, 4), (5, 4))]
    keys = ('exact_weighted_profiles', 'floating_api_exact_controls', 'constant_shift_controls',
            'rational_ramification_controls', 'conserved_charge_controls',
            'active_minor_subgradient_controls',
            'independent_weighted_quotient_ranks')
    totals = {key: sum(case[key] for case in cases) for key in keys}
    directory = Path(__file__).parent
    sources = ('fused_weighted_place_profiles.py', 'check_fused_weighted_places.py',
               'check_fused_closed_profiles.py', 'check_fused_multipoint_profiles.py',
               'fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py', 'f4_lattice_linear.py')
    result = {'schema_version': 1,
              'claim': 'Exact all-letter weighted fused profiles and bookkeeping controls; no improved factoring score.',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in sources},
              'cases': cases, 'totals': totals}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(totals, sort_keys=True))


if __name__ == '__main__':
    main()
