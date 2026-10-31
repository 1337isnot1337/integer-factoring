#!/usr/bin/env python3
"""Exact charge-gauge identities and single-sign profile majorization controls."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
import random

from check_rectangular_core_monomials import (
    core_monomial_seed, positive_monomials, occupation_charge,
)
from check_rectangular_defect_coordinates import original_cut, physical_masks
from check_rectangular_f4_repair import times_binary, binary_power, scale, fmul
from native_matrix_smith_f4 import MatrixOracle


def tensor_charges(seed):
    M, A, B = seed['ambient'], seed['q_x'], seed['q_u']
    assert A == B
    counts = [mask.bit_count() for mask in physical_masks(M, A)]
    return [counts[v]+counts[y]-counts[x]-counts[u]
            for v in range(A) for u in range(B) for y in range(B) for x in range(A)]


def check_gauge(original, transformed, coefficient=1, power=0):
    """Check K_transformed=k^charge*w^(power*charge)*K_original exactly."""
    def fpow(value, exponent):
        result = 1
        for _ in range(exponent % 3):
            result = fmul(result, value)
        return result
    checked = 0
    for left, right, charge in zip(transformed['values'], original['values'], tensor_charges(original)):
        left = times_binary(left, binary_power(3, original['denominator_degree']))
        right = times_binary(right, binary_power(3, transformed['denominator_degree']))
        left_power = original.get('global_w_power', 0)
        right_power = transformed.get('global_w_power', 0)+power*charge
        offset = -min(0, left_power, right_power)
        left = tuple(plane << (left_power+offset) for plane in left)
        right = tuple(plane << (right_power+offset) for plane in right)
        right = scale(right, fpow(coefficient, charge))
        assert left == right, (checked, charge, coefficient, power)
        checked += 1
    return checked


def profile(seed, place, row_weights, column_weights):
    matrix, _, _ = original_cut(seed, place, (0, 0, [0]*4, [0]*4, [0]*4))
    oracle = MatrixOracle(matrix, len(matrix))
    values = oracle.profile(row_weights, column_weights, infinity=place == 'infinity')['orders']
    global_w = seed.get('global_w_power', 0)
    correction = seed['denominator_degree']+global_w if place == 'infinity' else -global_w
    return [value+correction for value in values]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cubics = positive_monomials(degree=3)
    minus = [p for p in cubics if occupation_charge(p) == -1]
    plus = [p for p in cubics if occupation_charge(p) == 1]
    base = core_monomial_seed([], (0, 1), 2, full_alphabet_gamma=2, check_types=False)
    records = []
    cases = [[(*pattern, [(2, -1)], 0)] for pattern in cubics]
    cases += [[(*pattern, [(1+(i % 3), i-2)], 0) for i, pattern in enumerate(patterns)]
              for patterns in (minus, plus)]
    rng = random.Random(763419)
    for index, extras in enumerate(cases):
        changed = core_monomial_seed(extras, (0, 1), 2, full_alphabet_gamma=2, check_types=False)
        for place in ('infinity', 'zero'):
            for trial in range(2):
                rw = [rng.randrange(-12, 13)/4 for _ in range(16)] if trial else [0]*16
                cw = [rng.randrange(-12, 13)/4 for _ in range(16)] if trial else [0]*16
                baseline, perturbed = [profile(seed, place, rw, cw) for seed in (base, changed)]
                difference, prefixes = 0, []
                for old, new in zip(baseline, perturbed):
                    difference += new-old
                    prefixes.append(difference)
                assert max(prefixes) <= 1e-10 and abs(prefixes[-1]) < 1e-10
                records.append({'case': index, 'place': place, 'trial': trial,
                                'orders_checked': 16, 'minimum_prefix_difference': min(prefixes),
                                'maximum_prefix_difference': max(prefixes),
                                'total_difference': prefixes[-1]})
    gauges = []
    for negative, positive in product(minus, plus):
        for k in (2, 3):
            original = core_monomial_seed([(*negative, [(k, -1)], 0), (*positive, [(2, 2)], 0)],
                                           (0, 1), 2, full_alphabet_gamma=2, check_types=False)
            canonical = core_monomial_seed([(*negative, [(1, -1)], 0),
                                            (*positive, [(fmul(k, 2), 2)], 0)],
                                           (0, 1), 2, full_alphabet_gamma=2, check_types=False)
            gauges.append({'minus': negative, 'plus': positive, 'kind': 'field', 'coefficient': k,
                           'entries_checked': check_gauge(original, canonical, coefficient=k)})
        for a, b in ((-2, 3), (3, -1)):
            original = core_monomial_seed([(*negative, [(1, a)], 0), (*positive, [(2, b)], 0)],
                                           (0, 1), 2, full_alphabet_gamma=2, check_types=False)
            canonical = core_monomial_seed([(*negative, [(1, 0)], 0), (*positive, [(2, a+b)], 0)],
                                           (0, 1), 2, full_alphabet_gamma=2, check_types=False)
            gauges.append({'minus': negative, 'plus': positive, 'kind': 'Laurent', 'powers': [a, b],
                           'entries_checked': check_gauge(original, canonical, power=a)})
    result = {'status': 'exact finite controls of uniform full-alphabet charge reductions',
              'scope': 'Full physical alphabet and occupation-preserving base correction only.',
              'single_sign_profile_controls': records, 'charge_gauge_controls': gauges}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
