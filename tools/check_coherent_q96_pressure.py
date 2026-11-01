#!/usr/bin/env python3
"""Exact q96 coherent pressure witness and a scale-interval improvement.

The default checks saved complete profiles, exact maps, and block coverage.
Optional full-profile replay is expensive and explicit. Every pressure
decision uses bounded-mantissa outward integer intervals, without floats.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from check_coherent_q32_pressure import counting_polynomial, encode_interval, interval_controls, multiply, ratio_interval
from check_large_interval_pressure import BinaryIntervals
from coherent_hot_compounds import height, profiles as generic_profiles, validate_maps
from coherent_hot_compounds_fast import profiles as fast_profiles
from search_coherent_hot_retractions import _lift, q96_start


ROOT = Path(__file__).resolve().parents[1]
P, Q = 7100058652501, 10000000000000
ROOT_DENOMINATOR = 27855864849641


def full_profile_checks(profile):
    embedding, retraction, heights, charges, checks = validate_maps(
        profile['embedding'], profile['retraction'], 12)
    assert len(embedding) == len(retraction) == profile['q'] == 96
    assert list(heights) == profile['heights'] and list(charges) == profile['charges']
    expected = {(x, y) for x in range(96) for y in range(96)}
    for cut, key in zip(profile['cuts'], ('zero_orders_t', 'infinity_orders_t')):
        seen, orders = set(), []
        for block in cut['blocks']:
            pairs = set(map(tuple, block['pairs']))
            assert len(pairs) == block['size'] == len(block['orders_t'])
            assert not (seen & pairs)
            assert all(charges[x] + charges[y] == block['charge'] for x, y in pairs)
            assert block['endpoint_entries'] == block['size'] ** 2
            assert block['precision_attempts'][-1]['status'] == 0
            assert block['precision_attempts'][-1]['precision'] == block['precision_t']
            seen |= pairs
            orders.extend(block['orders_t'])
        assert seen == expected
        assert sorted(orders) == profile[key] and len(orders) == 9216
    return checks


def comparison(coordinate, coherent):
    assert sorted(coordinate['heights']) == sorted(coherent['heights'])
    first = multiply(counting_polynomial(coordinate['zero_orders_t']),
                     counting_polynomial(coordinate['infinity_orders_t']))
    second = multiply(counting_polynomial(coherent['zero_orders_t']),
                      counting_polynomial(coherent['infinity_orders_t']))
    difference = [(first[i] if i < len(first) else 0) - (second[i] if i < len(second) else 0)
                  for i in range(max(len(first), len(second)))]
    while difference and not difference[-1]:
        difference.pop()
    shift = next(i for i, value in enumerate(difference) if value)
    quotient, factors = difference[shift:], 0
    while sum(quotient) == 0:
        running, divided = 0, []
        for value in quotient[:-1]:
            running += value
            divided.append(running)
        assert running + quotient[-1] == 0
        quotient = divided
        factors += 1
    rebuilt = quotient
    for _ in range(factors):
        rebuilt = multiply([1, -1], rebuilt)
    assert [0] * shift + rebuilt == difference
    degree = len(quotient) - 1
    assert (shift, factors, degree, sum(quotient)) == (63, 1, 176, -8441856)
    # This is 5^degree Q(4x/5), in the ordinary power basis.
    scaled = [value * 4 ** j * 5 ** (degree - j) for j, value in enumerate(quotient)]
    bernstein = [sum(scaled[j] * comb(degree - j, i - j) for j in range(i + 1))
                 for i in range(degree + 1)]
    assert all(value >= 0 for value in bernstein)
    assert bernstein[0] > 0 and bernstein[-1] > 0
    expanded = [sum(bernstein[i] * (-1) ** (j - i) * comb(degree - i, j - i)
                    for i in range(j + 1)) for j in range(degree + 1)]
    assert expanded == scaled
    return {'difference_factorization': 'r^63*(1-r)*Q(r)',
            'quotient_power_coefficients': quotient,
            'scaled_quotient': '5^176*Q(4*x/5)',
            'nonnegative_unnormalized_bernstein_coefficients': bernstein,
            'strict_gain_interval': '0 < r <= 4/5',
            'Q_at_one': sum(quotient),
            'coordinate_better_in_some_neighborhood_below_one': True,
            'all_scale_domination_claimed': False}


def strict_greater_than_one(interval):
    lower, _, exponent = interval
    if exponent >= 0:
        return lower > 1 or (lower == 1 and exponent > 0)
    return lower > 1 << -exponent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coordinate', type=Path,
                        default=ROOT / 'evidence/coherent-q96-coordinate-reference.json')
    parser.add_argument('--candidate', type=Path,
                        default=ROOT / 'evidence/coherent-q96-four-boundary-candidate.json')
    parser.add_argument('--recompute-profiles', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    coordinate_receipt = json.loads(args.coordinate.read_text())
    candidate_receipt = json.loads(args.candidate.read_text())
    assert coordinate_receipt['completed'] and candidate_receipt['completed']
    coordinate, coherent = coordinate_receipt['profile'], candidate_receipt['profile']
    for receipt in (coordinate_receipt, candidate_receipt):
        for name, digest in receipt['source_sha256'].items():
            assert sha256((ROOT / 'tools' / name).read_bytes()).hexdigest() == digest, name
    words = sorted(range(1 << 12), key=lambda word: (height(word), word))[:96]
    expected = [1 << word for word in words]
    assert coordinate['embedding'] == coordinate['retraction'] == expected
    base, boundary, _, _ = q96_start()
    p = base + [_lift(value, boundary) for value in [3, 5, 9]]
    l = base + [_lift(value, boundary) for value in [13, 11, 8]]
    assert coherent['embedding'] == candidate_receipt['embedding'] == p
    assert coherent['retraction'] == candidate_receipt['retraction'] == l
    counts = Counter(deletion_checks=full_profile_checks(coordinate) + full_profile_checks(coherent),
                     complete_saved_profiles=4, complete_saved_smith_orders=36864,
                     independent_interval_controls=interval_controls())
    if args.recompute_profiles:
        for receipt, constructor in ((coordinate, generic_profiles), (coherent, fast_profiles)):
            current = constructor(receipt['embedding'], receipt['retraction'], 12)
            for key in ('zero_orders_t', 'infinity_orders_t'):
                assert current[key] == receipt[key]
                counts['complete_profiles_recomputed'] += 1
    coherent_ratio = ratio_interval(coherent, P, Q, ROOT_DENOMINATOR - 1, ROOT_DENOMINATOR)
    coordinate_ratio = ratio_interval(coordinate, P, Q, ROOT_DENOMINATOR - 1, ROOT_DENOMINATOR)
    lower, upper = BinaryIntervals.fraction_bounds(coherent_ratio)
    assert 0 < lower <= upper < Fraction(7, 20000) < 1
    assert strict_greater_than_one(coordinate_ratio)
    arithmetic = BinaryIntervals(256)
    root = arithmetic.rational(ROOT_DENOMINATOR - 1, ROOT_DENOMINATOR)
    r_lower, r_upper = arithmetic.fraction_bounds(arithmetic.power(root, P))
    assert 0 < r_lower <= r_upper < Fraction(4, 5)
    assert Fraction(P, Q) < Fraction(710005874881191, 10 ** 15)
    names = (Path(__file__).name, 'check_coherent_q32_pressure.py', 'check_large_interval_pressure.py',
             'coherent_hot_compounds.py', 'coherent_hot_compounds.cpp',
             'coherent_hot_compounds_fast.py', 'coherent_hot_compounds_fast.cpp',
             'search_coherent_hot_retractions.py', 'selected_gaussian_coefficients.cpp',
             'native_rectangular_smith_f4.cpp')
    report = {'status': 'pass', 'scope': __doc__, 'counts': dict(counts),
              'all_saved_profiles_recomputed': args.recompute_profiles,
              'rho_strict_upper_bound': str(Fraction(P, Q)),
              'root': [ROOT_DENOMINATOR - 1, ROOT_DENOMINATOR],
              'counting_variable_formula': f'r=root^{P}',
              'rational_physical_type_formula': f'p_j=root^({Q}*h_j)/sum_k root^({Q}*h_k)',
              'strict_pressure_inequality': f'(Z0(r)*Zinfinity(r))^{Q}<A(root^{Q})^{4 * P}',
              'coherent_ratio_interval': encode_interval(coherent_ratio),
              'coherent_ratio_rational_upper_bound': '7/20000',
              'coordinate_ratio_interval': encode_interval(coordinate_ratio),
              'coordinate_fails_same_rational_rate_and_scale': True,
              'fixed_mantissa_bits': 256,
              'huge_cleared_integer_powers_materialized': False,
              'scale_interval_comparison': comparison(coordinate, coherent),
              'input_sha256': {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                               for path in (args.coordinate, args.candidate)},
              'source_sha256': {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in names},
              'floating_arithmetic_used_for_certification': False,
              'global_optimality_claimed': False,
              'post_half_target_reached': 2 * P * P < Q * Q,
              'factoring_exponent_improved': False}
    args.output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps({'status': 'pass', 'rho_strict_upper_bound': report['rho_strict_upper_bound'],
                      'counts': dict(counts), 'strict_gain_interval': '0 < r <= 4/5',
                      'coherent_ratio_rational_upper_bound': '7/20000'}, indent=2))


if __name__ == '__main__':
    main()
