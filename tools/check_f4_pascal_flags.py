#!/usr/bin/env python3
"""Exact field, kernel inverse, binary specialization, and original-Fock controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_weighted_pascal_flags import subset_sums
from cube_f4_bivariate_smith import F4_MUL
from f4_pascal_flags import change_flag, clifford_data, inverse_flag, inverse_kernel, matrix_product
from search_weighted_pascal_flags import WeightedPascalPressure
from weighted_pascal_flags import (
    change_flag as binary_flag, clifford_data as binary_data, weighted_kernel,
)
from weighted_smith_oracle_f4 import WeightedSmithOracle, polynomial_multiply


def packed(pair):
    out = 0
    for plane, value in enumerate(pair):
        while value:
            bit = value & -value
            out ^= 1 << (2 * (bit.bit_length() - 1) + plane)
            value ^= bit
    return out


def scalar_product(first, second):
    n = len(first)
    result = [[0] * n for _ in range(n)]
    for i, j, k in product(range(n), repeat=3):
        result[i][j] ^= F4_MUL[first[i][k]][second[k][j]]
    return result


def full_tensor(seed):
    n = seed['n']
    matrix = [[packed(value) for value in row] for row in seed['matrix']]

    def minor(rows, columns):
        rows = [i for i in range(2 * n) if rows >> i & 1]
        columns = [i for i in range(2 * n) if columns >> i & 1]
        if len(rows) != len(columns):
            return 0
        dp = {0: 1}
        for row in rows:
            following = {}
            for mask, value in dp.items():
                for bit, column in enumerate(columns):
                    if not mask >> bit & 1 and matrix[row][column]:
                        following[mask | (1 << bit)] = following.get(mask | (1 << bit), 0) ^ polynomial_multiply(value, matrix[row][column])
            dp = following
        return dp.get((1 << len(columns)) - 1, 0)

    q = 1 << n
    return {(b, x, v, y): minor(v | (y << n), b | (x << n))
            for b, x, v, y in product(range(q), repeat=4)}


def original_profile_controls(seed, rng):
    n, q = seed['n'], 1 << seed['n']
    pressure = WeightedPascalPressure(seed)
    parameters = pressure.natural(1.)
    parameters[:-1] += rng.integers(-2, 3, size=len(parameters) - 1)
    pressure.evaluate(parameters)
    tensor = full_tensor(seed)
    controls = []
    for place, legs in enumerate(((0, 2, 1, 3), (0, 1, 2, 3))):
        matrix = [[0] * (q * q) for _ in range(q * q)]
        for labels, value in tensor.items():
            matrix[labels[legs[0]] + q * labels[legs[1]]][labels[legs[2]] + q * labels[legs[3]]] = value
        heights = parameters[place * 4 * n:(place + 1) * 4 * n].reshape(4, n).astype(int)
        symbol = [[sum(int(heights[leg, bit]) for bit in range(n) if word >> bit & 1)
                   for word in range(q)] for leg in range(4)]
        rows = [symbol[legs[0]][a] + symbol[legs[1]][b] for b in range(q) for a in range(q)]
        columns = [symbol[legs[2]][a] + symbol[legs[3]][b] for b in range(q) for a in range(q)]
        exact, _ = WeightedSmithOracle(matrix, infinity=bool(place)).profile(rows, columns)
        clifford = pressure.last['clifford_profiles'][place]['orders']
        assert all(value == int(value) for value in clifford)
        positive = list(map(int, clifford[2 * n:]))
        twice_mean = (2 * seed['outer'] * n if place == 0 else 0) + int(heights.sum())
        assert (twice_mean - sum(positive)) % 2 == 0
        minimum = (twice_mean - sum(positive)) // 2
        expected = [minimum + value for value in subset_sums(positive)]
        assert list(exact) == expected
        controls.append({'place': place, 'physical_mode_heights': heights.tolist(),
                         'complete_original_fock_orders': exact, 'clifford_orders': clifford})
    return {'constant_flag': seed['constant_flag'], 'controls': controls}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    original = weighted_kernel([1], 2)
    invertible, singular, binary_comparisons = [], 0, 0
    for values in product(range(4), repeat=4):
        constant = [list(values[:2]), list(values[2:])]
        try:
            inverse = inverse_flag(constant)
        except ValueError:
            singular += 1
            continue
        assert scalar_product(constant, inverse) == [[1, 0], [0, 1]]
        seed = change_flag(original, constant)
        identity = [[(int(i == j), 0) for j in range(4)] for i in range(4)]
        assert matrix_product(seed['matrix'], inverse_kernel(seed)) == identity
        assert matrix_product(inverse_kernel(seed), seed['matrix']) == identity
        if all(value < 2 for value in values):
            binary = binary_flag(original, constant)
            assert seed['matrix'] == [[(value, 0) for value in row] for row in binary['matrix']]
            for place in (0, 1):
                first, second = clifford_data(seed, place), binary_data(binary, place)
                assert first[1:] == second[1:]
                assert first[0] == [[(value, 0) for value in row] for row in second[0]]
                binary_comparisons += 1
        invertible.append(constant)
    assert len(invertible) == 180 and singular == 76
    rng = np.random.default_rng(2026110103)
    samples = [invertible[i] for i in range(0, 180, 13)]
    profiles = [original_profile_controls(change_flag(original, constant), rng)
                for constant in samples]
    numerical = []
    for constant in samples[:3]:
        pressure = WeightedPascalPressure(change_flag(original, constant))
        fitted = pressure.fitted_natural_start()
        pressure.evaluate(fitted)
        fitted_rho = pressure.last['rho']
        control = fitted.copy()
        control[:-1] += rng.normal(size=len(control) - 1) * .001
        _, derivative = pressure.evaluate(control, True)
        error = 0.
        for coordinate in (0, 3, 7, 8 * original['n']):
            delta = np.zeros_like(control)
            delta[coordinate] = 1e-7
            finite = (pressure.evaluate(control + delta) - pressure.evaluate(control - delta)) / 2e-7
            error = max(error, abs(float(finite - derivative[coordinate])))
        assert error < 1e-5
        numerical.append({'constant_flag': constant, 'fitted_rho': fitted_rho,
                          'gradient_max_error': error})
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/check_f4_pascal_flags.py', 'tools/f4_pascal_flags.py',
               'tools/search_weighted_pascal_flags.py', 'tools/weighted_pascal_flags.py',
               'tools/cube_f4_bivariate_smith.py', 'tools/weighted_smith_oracle_f4.py',
               'tools/check_weighted_pascal_flags.py', 'tools/native_matrix_smith_f4.py',
               'tools/native_matrix_smith_f4.cpp', 'tools/native_rectangular_smith_f4.cpp',
               'tools/native_weighted_smith_f4.py']
    report = {'status': 'pass', 'all_invertible_two_by_two_f4_flags': len(invertible),
              'all_singular_two_by_two_flags_rejected': singular,
              'binary_complete_clifford_specializations': binary_comparisons,
              'independent_original_profile_controls': profiles,
              'original_orders_compared': sum(len(cut['complete_original_fock_orders'])
                  for case in profiles for cut in case['controls']),
              'numerical_gradient_controls': numerical,
              'scope': 'Exact field/kernel identities and independent full-polynomial F4 '
              'Fock profiles; numerical scale and gradient controls are separate. '
              'No improved factoring exponent is claimed.',
              'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest()
                                for name in sources}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'all_invertible_two_by_two_f4_flags', 'original_orders_compared')}))


if __name__ == '__main__':
    main()
