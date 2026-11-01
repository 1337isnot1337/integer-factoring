#!/usr/bin/env python3
"""Literal endpoint, integer Smith, and numerical derivative controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from check_weighted_pascal_flags import exact_weighted_profile
from pascal_letter_gauges import PascalLetterPressure
from permuted_pascal_kernels import permuted_pascal_kernel
from search_weighted_pascal_flags import WeightedPascalPressure
from tensor_pascal_box_crossings import exterior_tensor


def endpoint_controls(seed):
    tensor = exterior_tensor(seed)
    q = len(tensor['alphabet'])
    endpoint = {key: value.bit_count() & 1 for key, value in tensor['values'].items()}
    energy = [sum(i for i in range(seed['n']) if symbol >> i & 1) for symbol in range(q)]
    nonzero = 0
    for b in range(q):
        for x in range(q):
            for v in range(q):
                for y in range(q):
                    value = endpoint.get((b, x, v, y), 0)
                    base = endpoint.get((0, x, v, 0), 0)
                    assert value == int(b == y) * base
                    if value:
                        assert x == v or energy[v] > energy[x]
                        nonzero += 1
    assert all(endpoint.get((0, x, x, 0), 0) == 1 for x in range(q))
    return {'alphabet': q, 'all_endpoint_entries_checked': q ** 4,
            'nonzero_entries': nonzero, 'strict_off_diagonal_energy': energy}


def integer_profile_controls(seed, alphabet):
    pressure = PascalLetterPressure(seed, alphabet=alphabet)
    q = pressure.q
    heights = [[[((11 * i + 5 * place + 3 * leg) % 13) - 6 for i in range(q)]
                for leg in range(4)] for place in range(2)]
    parameters = np.r_[np.asarray(heights).ravel(), 1.]
    pressure.evaluate(parameters)
    tensor = exterior_tensor(seed)
    index = {symbol: i for i, symbol in enumerate(alphabet)}
    controls = []
    for place, legs in enumerate(((0, 2, 1, 3), (0, 1, 2, 3))):
        size = q * q
        matrix = [[0] * size for _ in range(size)]
        for labels, value in tensor['values'].items():
            if all(symbol in index for symbol in labels):
                row = index[labels[legs[0]]] + q * index[labels[legs[1]]]
                column = index[labels[legs[2]]] + q * index[labels[legs[3]]]
                matrix[row][column] = value
        rows = [heights[place][legs[0]][a] + heights[place][legs[1]][b]
                for b in range(q) for a in range(q)]
        columns = [heights[place][legs[2]][a] + heights[place][legs[3]][b]
                   for b in range(q) for a in range(q)]
        exact = exact_weighted_profile(matrix, [[0] * size for _ in matrix],
                                       rows, columns, bool(place))
        assert pressure.last['complete_original_profiles'][place]['orders'] == exact
        controls.append({'place': place, 'integer_orders': exact})
    return {'images': seed['images'], 'alphabet': alphabet, 'symbol_heights': heights,
            'profiles': controls}


def numerical_controls(seed, rng):
    mode, letter = WeightedPascalPressure(seed), PascalLetterPressure(seed)
    cases = []
    for scale in (.07, .31):
        parameters = mode.natural(scale)
        parameters[:-1] += rng.normal(size=len(parameters) - 1) * .003
        lifted = letter.lift_mode_parameters(parameters)
        first, second = mode.evaluate(parameters), letter.evaluate(lifted)
        assert abs(first - second) < 1e-10
        changed = lifted.copy()
        changed[:-1] += rng.normal(size=len(changed) - 1) * .001
        value, derivative = letter.evaluate(changed, True)
        error = 0.
        for coordinate in (0, 1, letter.q, 2 * letter.q + 1, 4 * letter.q, 8 * letter.q):
            delta = np.zeros_like(changed)
            delta[coordinate] = 1e-7
            finite = (letter.evaluate(changed + delta) - letter.evaluate(changed - delta)) / 2e-7
            error = max(error, abs(float(finite - derivative[coordinate])))
        assert error < 1e-5
        cases.append({'scale': scale, 'mode_rho': first, 'lifted_rho': second,
                      'nonadditive_symbol_rho': value, 'gradient_max_error': error})
    return {'images': seed['images'], 'cases': cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    seeds = [permuted_pascal_kernel(images) for images in
             ([0, 1], [1, 0], [0, 1, 2, 3], [0, 1, 3, 2])]
    endpoint = [endpoint_controls(seed) for seed in seeds]
    integer = [integer_profile_controls(seeds[0], [0, 1, 2, 3]),
               integer_profile_controls(seeds[1], [0, 1, 3]),
               integer_profile_controls(seeds[2], list(range(16))),
               integer_profile_controls(seeds[2], [0, 3, 5, 6, 9, 10, 12, 15]),
               integer_profile_controls(seeds[3], [0, 1, 2, 4, 8])]
    rng = np.random.default_rng(2026110101)
    numerical = [numerical_controls(seed, rng) for seed in seeds]
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/check_pascal_letter_gauges.py', 'tools/pascal_letter_gauges.py',
               'tools/search_pascal_letter_gauges.py', 'tools/permuted_pascal_kernels.py',
               'tools/weighted_pascal_flags.py', 'tools/search_weighted_pascal_flags.py',
               'tools/tensor_pascal_box_crossings.py', 'tools/check_weighted_pascal_flags.py',
               'tools/cube_symmetric_exterior_blocks.py', 'tools/native_matrix_smith_f4.py',
               'tools/native_matrix_smith_f4.cpp', 'tools/native_rectangular_smith_f4.cpp',
               'tools/native_weighted_smith_f4.py']
    report = {'status': 'pass', 'endpoint_controls': endpoint,
              'exact_integer_profiles': integer, 'numerical_controls': numerical,
              'exact_original_orders_compared': sum(len(c['integer_orders'])
                  for case in integer for c in case['profiles']),
              'scope': 'Exact finite endpoints and independent original-polynomial integer '
              'Smith profiles; separate floating mode/letter and derivative controls. '
              'No stronger factoring exponent is claimed.',
              'source_sha256': {path: sha256((root / path).read_bytes()).hexdigest()
                                for path in sources}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass', 'exact_original_orders_compared':
                      report['exact_original_orders_compared'], 'numerical_cases': 8}))


if __name__ == '__main__':
    main()
