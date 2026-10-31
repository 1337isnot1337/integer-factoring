#!/usr/bin/env python3
"""Replay the compiled F4 Weyl filter engine against complete active cores."""
from __future__ import annotations

import argparse
from itertools import product
import json
from math import exp, log, sqrt
from pathlib import Path
from random import Random

from coupled_identity_operator_filters import _geometry
from f4_repaired_identity_filters import (
    compiled_monomial_profiles, monomial_filters, zero_profiles,
)
from prime_state_filter_pressure import gaussian_pressure
from search_f4_repaired_identity_filters import pressure_function, pressure_grid


def profile_key(blocks):
    return tuple((block['charge'], block['orders']) for block in blocks)


def compiled_controls():
    rng = Random(8361127)
    records = []
    for n, ramification, modes, geometry in ((3, 1, 2, 'initial'), (8, 3, 2, 'mixed'),
                                             (8, 1, 2, 'terminal'), (4, 1, 3, 'initial'),
                                             (8, 1, 3, 'terminal')):
        for variant in range(6):
            exponents = (0,) + tuple(rng.randrange(-3, 4) for _ in range(2 ** modes - 1))
            coefficients = (1,) + tuple(rng.randrange(1, 4) for _ in range(2 ** modes - 1))
            alpha = monomial_filters(modes, exponents, coefficients)
            actual = compiled_monomial_profiles(n, ramification, modes, exponents, geometry)(coefficients)
            expected = zero_profiles(n, ramification, modes, alpha, geometry)
            assert profile_key(actual) == profile_key(expected)
            records.append({'n': n, 'ramification': ramification, 'modes': modes,
                            'geometry': geometry, 'exponents': exponents, 'coefficients': coefficients,
                            'full_active_matrix_order': 2 ** (2 * modes),
                            'reduced_matrix_order': 2 ** modes,
                            'maximum_exact_jet_precision': max(
                                block['exact_jet_precision'] for block in actual)})
    exponents, coefficients = (0, 1, 600, 601), (1, 2, 3, 1)
    actual = compiled_monomial_profiles(2, 1, 2, exponents)(coefficients)
    expected = zero_profiles(2, 1, 2, monomial_filters(2, exponents, coefficients))
    assert profile_key(actual) == profile_key(expected)
    assert max(block['exact_jet_precision'] for block in actual) >= 1024
    return {'complete_active_profile_replays': records,
            'late_exponent_complete_source_replay': {'exponents': exponents, 'coefficients': coefficients,
                'maximum_exact_jet_precision': max(block['exact_jet_precision'] for block in actual)}}


def pressure_controls():
    records = 0
    maximum_error = 0.
    for n, modes, geometry in ((3, 2, 'initial'), (4, 2, 'mixed'), (4, 3, 'initial'),
                                (4, 3, 'terminal')):
        states = tuple(product(range(2), repeat=modes))
        exponents = tuple(sum(state) ** 2 - 2 * sum(state) for state in states)
        coefficients = tuple(1 + i % 3 for i in range(len(states)))
        blocks = compiled_monomial_profiles(n, 1, modes, exponents, geometry)(coefficients)
        evaluate = pressure_function(n, 1, modes, geometry, exponents, blocks)
        data = _geometry(2, n, modes, geometry)
        selected, retained = data['selected'], data['retained']
        for charge, scale in product((-2, 0, 2, 3), (.03, .15, .6)):
            active_zero = [value + (1 - charge) * block['charge']
                           for block in blocks for value in block['orders']]
            active_infinity = [sum((2 * j + 1 - 2 * n) * digit for j, digit in zip(selected, state))
                               + 2 * exponent for state, exponent in zip(states, exponents)]
            zero = active_zero
            for i in range(n):
                if i not in selected:
                    zero = [value + digit * (2 * i + 1 - charge) for value in zero for digit in range(2)]
                if i not in retained:
                    zero = [value + digit * (2 * i + 1 + charge) for value in zero for digit in range(2)]
            infinity = active_infinity
            for i in range(2 * n):
                if i not in selected:
                    infinity = [value + digit * (2 * i + 1 - 2 * n)
                                for value in infinity for digit in range(2)]
            physical = [0]
            for i in range(n):
                for sign in (-1, 1):
                    energy = 4 * i - 2 * (n - 1) + sign * charge
                    physical = [value + digit * energy for value in physical for digit in range(2)]
            direct = (log(sum(exp(-scale * value) for value in zero))
                      + log(sum(exp(-scale * value) for value in infinity))
                      - sqrt(2) * log(sum(exp(-scale * value / sqrt(2)) for value in physical)))
            actual = evaluate(log(scale), charge)
            grid, relative = pressure_grid(n, 1, geometry, (scale,), (charge,), modes)(exponents, blocks)
            maximum_error = max(maximum_error, abs(actual - direct), abs(grid[0] - direct))
            assert abs(actual - direct) < 1e-10 and abs(grid[0] - direct) < 1e-10
            assert abs(relative[0] - actual + gaussian_pressure(2, n, charge, scale)) < 1e-10
            records += 1
    return {'complete_profile_and_physical_partition_comparisons': records,
            'maximum_absolute_floating_error': maximum_error,
            'comparison_tolerance': 1e-10}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'status': 'Exact compiled F4 Weyl profiles and independent pressure controls; no witness claim.',
              'compiled_profiles': compiled_controls(), 'pressure_assembly': pressure_controls()}
    encoded = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
