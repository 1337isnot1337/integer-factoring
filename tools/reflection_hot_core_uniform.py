#!/usr/bin/env python3
"""Bivariate complete matrices for uniform one-hot reflection profiles.

T represents z^(2M). Multiplication by the local unit 1+z² removes
the geometric-series denominator. A global degree bound on ALL minors
then lets a sufficiently separated specialization decode their exact
lexicographic (T,z) orders.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import itertools
import json
from pathlib import Path

from check_rectangular_f4_repair import fmul
from reflection_hot_core_profiles import compile_interval
from search_reflection_hot_core_filters import polynomial_specs


ONE = {(0, 0): 1}


def add(first, second):
    answer = dict(first)
    for degree, value in second.items():
        value ^= answer.get(degree, 0)
        if value:
            answer[degree] = value
        else:
            answer.pop(degree, None)
    return answer


def multiply(first, second):
    answer = {}
    for (a, b), x in first.items():
        for (c, d), y in second.items():
            degree = a + c, b + d
            value = answer.get(degree, 0) ^ fmul(x, y)
            if value:
                answer[degree] = value
            else:
                answer.pop(degree, None)
    return answer


def monomial(z=0, t=0):
    return {(z, t): 1}


def act(matrix, creation, coefficients):
    answer = [[{} for _ in range(8)] for _ in range(8)]
    for mode, coefficient in enumerate(coefficients):
        if not coefficient:
            continue
        bit = 1 << mode
        for source in range(8):
            if bool(source & bit) == creation:
                continue
            target = source ^ bit
            for column, value in enumerate(matrix[source]):
                if value:
                    answer[target][column] = add(answer[target][column], multiply(coefficient, value))
    return answer


def components(place):
    c = {(0, 0): 1, (2, 0): 1}
    # The right-first operator is scaled by c. Other monomials get
    # the same scalar after all actions, so every returned component
    # is precisely c times the original component.
    images = {
        'left_first': (True, ({}, ONE, {})),
        'left_second': (True, (ONE, c, multiply(c, monomial(4)))),
        'right_first': (False, (multiply(c, monomial(-2)),
                               multiply(multiply(c, c), monomial(-2, -1)),
                               {(-2, -1): 1, (-6, 0): 1})),
        'right_second': (False, ({}, monomial(-2, -1), {})),
    }
    initial = [[{} for _ in range(8)] for _ in range(8)]
    for state in range(8):
        a, b, d = ((state >> i) & 1 for i in range(3))
        initial[state][state] = monomial(2 * a + 2 * b + 6 * d, b + d)
    left = ('x', 'y') if place == 'zero' else ('u', 'x')
    right = ('u', 'v') if place == 'zero' else ('v', 'y')
    names = {left[0]: 'left_first', left[1]: 'left_second',
             right[0]: 'right_first', right[1]: 'right_second'}
    gauged = (0, 1) if place == 'zero' else (2, 3)
    answer = {}
    for mask in (0, 3, 5, 6, 9, 10, 12, 15):
        matrix = initial
        has_scaled = False
        for letter in right + left:
            if mask >> 'xyuv'.index(letter) & 1:
                name = names[letter]
                has_scaled |= name == 'right_first'
                matrix = act(matrix, *images[name])
        scalar = monomial(2 * sum((mask >> leg) & 1 for leg in gauged))
        if not has_scaled:
            scalar = multiply(scalar, c)
        answer[mask] = [[multiply(scalar, value) for value in row] for row in matrix]
    return answer


def sources(place):
    base = components(place)
    sign = 1 if place == 'zero' else -1
    answer = [base[0]]
    for mask, polynomial in polynomial_specs(1):
        scalar = {(sign * power, 0): value for power, value in polynomial.items()}
        answer.append([[multiply(scalar, value) for value in row] for row in base[mask]])
    return answer


def specialize(polynomial, separation):
    answer = {}
    for (z, t), value in polynomial.items():
        degree = z + separation * t
        value ^= answer.get(degree, 0)
        if value:
            answer[degree] = value
        else:
            answer.pop(degree, None)
    return answer


def degree_bounds(place):
    powers = [degree for source in sources(place) for row in source for value in row for degree in value]
    zmin, zmax = min(a for a, b in powers), max(a for a, b in powers)
    tmin, tmax = min(b for a, b in powers), max(b for a, b in powers)
    zclear, tclear = max(0, -zmin), max(0, -tmin)
    width = zmax + zclear
    return {'place': place, 'entry_z_range': (zmin, zmax), 'entry_T_range': (tmin, tmax),
            'uniform_z_clearing': zclear, 'uniform_T_clearing': tclear,
            'all_k_minor_z_interval': [(0, k * width) for k in range(1, 9)],
            'sufficient_strict_separation_bound': 8 * width}


def decode_profile(orders, separation, bounds):
    assert separation > bounds['sufficient_strict_separation_bound']
    zclear, tclear = bounds['uniform_z_clearing'], bounds['uniform_T_clearing']
    prefix = 0
    old_z = old_t = 0
    result = []
    for k, order in enumerate(orders, 1):
        prefix += order
        cleared = prefix + k * (zclear + separation * tclear)
        t, z = divmod(cleared, separation)
        assert 0 <= z <= bounds['all_k_minor_z_interval'][k - 1][1]
        assert 0 <= t <= k * (bounds['entry_T_range'][1] + tclear)
        result.append((z - old_z - zclear, t - old_t - tclear))
        old_z, old_t = z, t
    return tuple(result)


def pair_affines(profiles):
    return [(a + c, b + d) for a, b in profiles[0] for c, d in profiles[1]]


def verify(record_path):
    record = json.loads(record_path.read_text())
    assert record['complete'] and record['ramification'] == 1
    assert record['completed_candidates'] == 4 ** 8 and not record['unresolved']
    M = record['interval']
    assert isinstance(M, int) and M >= 4 and M & (M - 1) == 0
    assert sum(row['rows'] for row in record['profile_classes']) == 4 ** 8
    assert len({tuple(tuple(profile) for profile in row['profiles'])
                for row in record['profile_classes']}) == len(record['profile_classes'])
    for row in record['profile_classes']:
        assert len(row['profiles']) == 2
        for profile in row['profiles']:
            assert len(profile) == 8 and profile == sorted(profile)
            assert all(isinstance(order, int) for order in profile)
    bounds = [degree_bounds(place) for place in ('zero', 'infinity')]
    controls = 0
    # The bivariate constructor is separate from the geometric-sum engine.
    for interval in (4, 8, 16, 64):
        for place in ('zero', 'infinity'):
            old = compile_interval(interval, place)['components']
            for mask, matrix in components(place).items():
                for i, row in enumerate(matrix):
                    for j, value in enumerate(row):
                        expected = {}
                        for degree, coefficient in old[mask][i][j].items():
                            for shift in (0, 2):
                                power = degree + shift
                                expected[power] = expected.get(power, 0) ^ coefficient
                        expected = {degree: value for degree, value in expected.items() if value}
                        assert specialize(value, 2 * interval) == expected
                        controls += 1
    baseline = next(row for row in record['profile_classes'] if not any(row['representative']))
    gaussian = [decode_profile(profile, 2 * M, bound)
                for profile, bound in zip(baseline['profiles'], bounds)]
    base_pairs = sorted(pair_affines(gaussian), key=lambda x: (x[1], x[0]))
    certificates = []
    for row in record['profile_classes']:
        affine = [decode_profile(profile, 2 * M, bound)
                  for profile, bound in zip(row['profiles'], bounds)]
        pairs = sorted(pair_affines(affine), key=lambda x: (x[1], x[0]))
        # This fixed pairing proves the inequality for every larger M.
        differences = [(b[0] - a[0], b[1] - a[1]) for a, b in zip(pairs, base_pairs)]
        assert all(slope >= 0 and intercept + 2 * M * slope >= 0 for intercept, slope in differences)
        certificates.append({'representative': row['representative'], 'rows': row['rows'],
                             'affine_profiles_z_intercept_T_slope': affine,
                             'minimum_slope_margin': min(b for a, b in differences),
                             'minimum_value_margin_at_base': min(a + 2 * M * b for a, b in differences)})
    return {'status': 'Exact bivariate lifting and all-larger-interval order domination.',
            'base_interval': M, 'coefficient_rows': sum(row['rows'] for row in certificates),
            'profile_classes': len(certificates), 'coefficientwise_component_controls': controls,
            'degree_bounds': bounds, 'certificates': certificates,
            'record_sha256': sha256(record_path.read_bytes()).hexdigest(),
            'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                              for name in (Path(__file__).name, 'reflection_hot_core_profiles.py',
                                           'search_reflection_hot_core_filters.py')},
            'scope': 'Every power-of-two M at least the base interval, every s>1, the stated F4^8 polynomial family at natural gauges.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path, nargs='?', default=Path('evidence/reflection-hot-grid-M128.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify(args.record)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'certificates'}))


if __name__ == '__main__':
    main()
