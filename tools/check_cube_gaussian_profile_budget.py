#!/usr/bin/env python3
"""Exact endpoint-minor budgets for arbitrary homogeneous eight-mode kernels."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import random
import time

from cube_homogeneous_gaussian import (
    ZERO, axis_parts, determinantal_orders, minor, monomial, valuation,
)


def integer_product(first, second):
    result = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        if a:
            for j, b in enumerate(second):
                if b:
                    result[i + j] += a * b
    return result


def floor_certificate():
    partition = [1]
    for degree in (6, 12, 27, 27):
        factor = [0] * (degree + 1)
        factor[0] = factor[degree] = 1
        partition = integer_product(partition, factor)
    difference = [1]
    for _ in range(5):
        difference = integer_product(difference, partition)
    for i in range(13):
        difference[10 * i] -= comb(12, i)
    order = 0
    while difference[0] == 0:
        difference = difference[1:]
        order += 1
    assert order == 6
    degree = len(difference) - 1
    transformed = [sum(coefficient * comb(degree - i, j - i)
                       for i, coefficient in enumerate(difference[:j + 1]))
                   for j in range(degree + 1)]
    assert degree == 354 and min(transformed) == 5
    assert all(value > 0 for value in transformed)
    return {'rho_floor': '3/5', 'quotient_degree': degree, 'vanishing_order': order,
            'difference': '[(1+s^6)(1+s^12)(1+s^27)^2]^5 - (1+s^10)^12',
            'transform': '(1+u)^354 * (difference/s^6)(u/(1+u))',
            'positive_integer_coefficients': transformed}


def sampled_kernel(degree, rng):
    monomials = [(a, b, degree - a - b)
                 for a in range(degree + 1) for b in range(degree - a + 1)]
    source = [[ZERO] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(i + 1, 8):
            value = frozenset(term for term in monomials if rng.randrange(2))
            if len(value) % 2 != int(i ^ j == 7):
                value ^= monomial(degree, 0, 0)
            source[i][j] = source[j][i] = value
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='evidence/cube-gaussian-profile-budget.json')
    args = parser.parse_args()
    started = time.monotonic()
    rng, records, minor_checks = random.Random(2026110103), [], 0
    for degree in range(1, 5):
        for case in range(8):
            source = sampled_kernel(degree, rng)
            assert all(len(source[i][j]) % 2 == int(i ^ j == 7)
                       for i in range(8) for j in range(8))
            profiles = []
            for axis in range(3):
                first, second = axis_parts(axis)
                block = [[source[i][j] for j in second] for i in first]
                profile, checked = determinantal_orders(block, axis)
                profiles.append(profile)
                minor_checks += checked
            first_witness = source[0][7]
            assert first_witness and len(first_witness) % 2 == 1
            assert sum(profile[0] for profile in profiles) <= sum(valuation(first_witness, axis) for axis in range(3)) <= degree
            pair_witnesses = []
            for omitted in range(3):
                axes = [axis for axis in range(3) if axis != omitted]
                other = 1 << omitted
                witness = minor(source, [0, other], [7, other ^ 7])
                assert witness and len(witness) % 2 == 1
                left = sum(sum(profiles[axis][:2]) for axis in axes)
                right = sum(valuation(witness, axis) for axis in axes)
                assert left <= right <= 2 * degree
                pair_witnesses.append({'axes': axes, 'profile_prefix_total': left,
                                       'minor_valuation_total': right})
            assert sum(sum(profile[:2]) for profile in profiles) <= 3 * degree
            assert sum(sum(profile) for profile in profiles) <= 12 * degree
            records.append({'degree': degree, 'sample': case, 'profiles': profiles,
                            'pair_witnesses': pair_witnesses})
    pressure = floor_certificate()
    paths = [Path(__file__), Path(__file__).with_name('cube_homogeneous_gaussian.py')]
    report = {'schema': 'cube-gaussian-profile-budget-v1',
              'scope': 'A bound on common unweighted Smith-profile pressure; not a rank lower bound for arbitrary pieces.',
              'sampled_common_kernels': len(records), 'literal_square_minors': minor_checks,
              'records': records, 'pressure_floor_certificate': pressure,
              'new_factoring_score': False,
              'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in paths},
              'elapsed_seconds': time.monotonic() - started}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'kernels': len(records), 'minors': minor_checks,
                      'positive_coefficients': len(pressure['positive_integer_coefficients']),
                      'seconds': report['elapsed_seconds']}))


if __name__ == '__main__':
    main()
