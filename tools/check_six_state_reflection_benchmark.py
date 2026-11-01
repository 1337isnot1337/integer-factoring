#!/usr/bin/env python3
"""Exact q6 benchmark replay, historical identification, and Gibbs root bracket.

The bracket is for the displayed rational natural scale, not a global
minimum over scale or physical heights. The result reproduces an earlier
reflection-code benchmark and does not improve the factoring exponent.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import time

from balanced_particle_extensions import construct
from charged_gaussian_floor import Dyadic
from check_balanced_particle_extensions import certify_tensor, local_integer_matrix
from check_f4_recentered_walks import certify_profile
from search_f4_recentered_crossings import digest


def log_positive(arithmetic, value):
    assert value[0] >= arithmetic.scale
    power = 0
    while value[0] >= (arithmetic.scale << (power+1)):
        power += 1
    reduced = arithmetic.divide_integer(value, 1 << power)
    logarithm = arithmetic.log_one_plus(arithmetic.subtract(reduced, arithmetic.one))
    return arithmetic.add(logarithm, arithmetic.multiply(arithmetic.rational(power),
                                                        arithmetic.log_one_plus(arithmetic.one)))


def partition(arithmetic, energy, scale):
    minimum = min(energy)
    total = arithmetic.zero
    for value in energy:
        total = arithmetic.add(total, arithmetic.exp_negative(
            arithmetic.multiply(arithmetic.rational(value-minimum), scale)))
    log_value = arithmetic.subtract(log_positive(arithmetic, total),
                                    arithmetic.multiply(arithmetic.rational(minimum), scale))
    return total, log_value


def pressure(arithmetic, profiles, energies_x, energies_y, scale, rho):
    rate = arithmetic.rational(2*rho.numerator, rho.denominator)
    result = arithmetic.zero
    for profile in profiles:
        result = arithmetic.add(result, partition(arithmetic, profile, scale)[1])
    physical = arithmetic.zero
    for energies in (energies_x, energies_y):
        physical = arithmetic.add(physical, partition(arithmetic, energies,
                                  arithmetic.divide(scale, rate))[1])
    return arithmetic.subtract(result, arithmetic.multiply(rate, physical))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--natural', type=Path, default=Path('evidence/balanced-particle-extension-natural-search.json'))
    parser.add_argument('--historical', type=Path, default=Path('evidence/reflection-code-alphabet-search.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    natural = json.loads(args.natural.read_text())
    row = next(item for item in natural['results'] if item['family'] == 'symmetric_6x6' and item['added_mask'] == 85)
    seed = construct(row['code_x'], row['code_y'])
    assert digest(seed['values']) == row['digest']
    minors = certify_tensor(seed)
    for side in (0, 1):
        matrix, scalar = local_integer_matrix(seed, side)
        certify_profile(matrix, scalar, row['profiles'][side])

    historical = json.loads(args.historical.read_text())
    old = next(item for item in historical['results'] if set(item['alphabet']) == {0, 2, 4, 8, 10, 12})
    masks = [sum(1 << (7-i if (letter >> i) & 1 else i) for i in range(4))
             for letter in old['alphabet']]
    assert set(masks) == set(seed['code_x'])
    old_profiles = [sorted(int(value) for value, count in profile.items() for _ in range(count))
                    for profile in old['profiles']]
    assert old_profiles[0] == [2*value+8 for value in row['profiles'][1]]
    assert old_profiles[1] == [2*value-8 for value in row['profiles'][0]]
    assert sorted(old['energies']) == sorted(2*value for value in seed['energies_x'])

    arithmetic = Dyadic(100)
    scale_value = Fraction(551453777, 10**9)
    scale = arithmetic.rational(scale_value.numerator, scale_value.denominator)
    low, high = Fraction(7226932545,10**10), Fraction(7226932547,10**10)
    low_pressure = pressure(arithmetic, row['profiles'], seed['energies_x'], seed['energies_y'], scale, low)
    high_pressure = pressure(arithmetic, row['profiles'], seed['energies_x'], seed['energies_y'], scale, high)
    assert low_pressure[0] > 0 > high_pressure[1]
    rate_interval = (arithmetic.rational(2*low.numerator, low.denominator)[0],
                     arithmetic.rational(2*high.numerator, high.denominator)[1])
    physical_scale = arithmetic.divide(scale, rate_interval)
    total, _ = partition(arithmetic, seed['energies_x'], physical_scale)
    added_index = seed['code_x'].index(85)
    numerator = arithmetic.exp_negative(arithmetic.multiply(
        arithmetic.rational(seed['energies_x'][added_index]), physical_scale))
    probability = arithmetic.divide(numerator, total)
    assert probability[0] > arithmetic.rational(623,100000)[1]
    assert probability[1] < arithmetic.rational(625,100000)[0]
    report = {'status': 'passed',
              'scope': 'Exact original Gaussian minors, complete natural Smith profiles, and fixed-rational-scale Gibbs bracket for a previously retained q6 reflection benchmark; no new exponent or global minimum.',
              'natural_record': str(args.natural), 'natural_record_sha256': sha256(args.natural.read_bytes()).hexdigest(),
              'historical_record': str(args.historical), 'historical_record_sha256': sha256(args.historical.read_bytes()).hexdigest(),
              'historical_name': old['name'], 'historical_reflection_alphabet': old['alphabet'],
              'physical_masks_in_historical_order': masks, 'historical_profile_identity':
                  'old_zero=2*new_zero+8; old_infinity=2*new_infinity-8; old_energy=2*new_energy after permutation',
              'historical_numerical_rho': old['rho'], 'literal_gaussian_minors': minors,
              'complete_profiles_certified_by_regular_frames': 2,
              'literal_seed': seed, 'coefficient_digest': row['digest'], 'profiles': row['profiles'],
              'natural_log_scale': str(scale_value), 'gibbs_rho_bracket': [str(low), str(high)],
              'dyadic_denominator': str(arithmetic.scale),
              'lower_endpoint_pressure_interval_numerators': list(map(str, low_pressure)),
              'upper_endpoint_pressure_interval_numerators': list(map(str, high_pressure)),
              'added_mask': 85, 'added_mass_strict_bracket': ['623/100000','625/100000'],
              'added_mass_interval_numerators': list(map(str, probability)),
              'factoring_exponent_changed': False,
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'elapsed_seconds': time.monotonic()-started}
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: value for key, value in report.items() if key != 'literal_seed'}), flush=True)


if __name__ == '__main__':
    main()
