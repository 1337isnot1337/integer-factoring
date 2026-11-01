#!/usr/bin/env python3
"""Independent determinant, full-profile, scale, and derivative controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from pascal_coordinate_sectors import coordinate_exterior_tensor, SectorPascalLetterPressure
from pascal_letter_gauges import PascalLetterPressure
from permuted_pascal_kernels import permuted_pascal_kernel
from tensor_pascal_box_crossings import exterior_tensor, polynomial_product


def independent_minor(matrix, rows, columns):
    rows = [i for i in range(len(matrix)) if rows >> i & 1]
    columns = [j for j in range(len(matrix)) if columns >> j & 1]
    if len(rows) != len(columns):
        return 0
    dp = {0: 1}
    for row in rows:
        following = {}
        for used, value in dp.items():
            for bit, column in enumerate(columns):
                if not used >> bit & 1 and matrix[row][column]:
                    mask = used | (1 << bit)
                    following[mask] = following.get(mask, 0) ^ polynomial_product(value, matrix[row][column])
        dp = following
    return dp.get((1 << len(columns)) - 1, 0)


def entry_controls():
    cases = []
    for images in ([0, 1], [0, 1, 2, 3], [0, 1, 3, 2]):
        seed = permuted_pascal_kernel(images)
        original = exterior_tensor(seed)
        selected = coordinate_exterior_tensor(seed, original['alphabet'])
        assert all(selected['values'].get(labels, 0) == value
                   for labels, value in original['values'].items())
        cases.append({'images': images, 'all_original_polynomial_entries': len(original['values'])})
    seed = permuted_pascal_kernel(list(range(8)))
    alphabet = [0] + [1 << i for i in range(8)]
    selected = coordinate_exterior_tensor(seed, alphabet)
    count = 0
    for b, x, v, y in product(alphabet, repeat=4):
        expected = independent_minor(seed['matrix'], v | (y << 8), b | (x << 8))
        assert selected['values'].get((b, x, v, y), 0) == expected
        count += 1
    cases.append({'images': list(range(8)), 'alphabet': alphabet,
                  'independent_subset_dynamic_program_entries': count})
    return cases


def profile_controls(rng):
    cases = []
    for images, alphabet in [([0, 1], list(range(4))),
                             ([0, 1, 2, 3], list(range(16))),
                             ([0, 1, 3, 2], [0, 1, 2, 4, 8])]:
        seed = permuted_pascal_kernel(images)
        dense = PascalLetterPressure(seed, alphabet=alphabet)
        sectors = SectorPascalLetterPressure(seed, alphabet)
        comparisons = []
        for trial in range(3):
            heights = rng.integers(-5, 6, size=8 * len(alphabet)).astype(float)
            parameters = np.r_[heights, 1.]
            first, second = dense.evaluate(parameters), sectors.evaluate(parameters)
            assert abs(first - second) < 1e-10
            count = 0
            for actual, expected in zip(sectors.last['complete_original_profiles'],
                                        dense.last['complete_original_profiles']):
                assert actual['orders'] == expected['orders']
                count += len(actual['orders'])
            comparisons.append({'integer_orders_compared': count, 'dense_rho': first,
                                'sector_rho': second, 'parameters': parameters.tolist()})
        start = sectors.natural(.19)
        start[:-1] += rng.normal(size=len(start) - 1) * .001
        value, derivative = sectors.evaluate(start, True)
        error = 0.
        for coordinate in (0, sectors.q, 4 * sectors.q, 8 * sectors.q):
            delta = np.zeros_like(start)
            delta[coordinate] = 1e-7
            finite = (sectors.evaluate(start + delta) - sectors.evaluate(start - delta)) / 2e-7
            error = max(error, abs(float(finite - derivative[coordinate])))
        assert error < 1e-5
        cases.append({'images': images, 'alphabet': alphabet, 'integer_comparisons': comparisons,
                      'nonadditive_rho': value, 'gradient_max_error': error,
                      'sector_sizes': [[len(rows) for _, rows, _, _ in oracle.blocks]
                                       for oracle, _ in sectors.cuts]})
    return cases


def scale_controls():
    seed = permuted_pascal_kernel([0, 1, 2, 3])
    dense = PascalLetterPressure(seed)
    sectors = SectorPascalLetterPressure(seed, list(range(16)))
    first, second = dense.fitted_natural_start(), sectors.fitted_natural_start()
    assert max(abs(first - second)) < 1e-8
    dense.evaluate(first)
    sectors.evaluate(second)
    assert abs(dense.last['rho'] - sectors.last['rho']) < 1e-10
    assert abs(sectors.last['rho'] - .7310809359736) < 1e-10
    base = sectors.natural(1.)
    sectors.evaluate(base)
    profiles = [cut['orders'] for cut in sectors.last['complete_original_profiles']]
    for scale in (.125, .5, 2.):
        sectors.evaluate(base * scale)
        assert all(actual['orders'] == [scale * value for value in expected]
                   for actual, expected in zip(sectors.last['complete_original_profiles'], profiles))
    return {'matched_natural_rho': dense.last['rho'], 'positive_scale_profile_controls': 3}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    entries = entry_controls()
    profiles = profile_controls(np.random.default_rng(2026110102))
    scale = scale_controls()
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/check_pascal_coordinate_sectors.py', 'tools/pascal_coordinate_sectors.py',
               'tools/pascal_letter_gauges.py', 'tools/search_pascal_letter_gauges.py',
               'tools/permuted_pascal_kernels.py', 'tools/tensor_pascal_box_crossings.py',
               'tools/weighted_pascal_flags.py', 'tools/native_matrix_smith_f4.py',
               'tools/native_matrix_smith_f4.cpp', 'tools/native_rectangular_smith_f4.cpp',
               'tools/native_weighted_smith_f4.py']
    report = {'status': 'pass', 'entry_controls': entries, 'profile_controls': profiles,
              'scale_controls': scale,
              'integer_orders_compared': sum(c['integer_orders_compared']
                  for case in profiles for c in case['integer_comparisons']),
              'scope': 'Exact selected-entry and complete integer-profile comparisons; '
              'separate numerical pressure, positive-scale, and derivative controls. '
              'No stronger factoring bound is claimed.',
              'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest()
                                for name in sources}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass', 'integer_orders_compared': report['integer_orders_compared'],
                      'entry_cases': len(entries), 'gradient_cases': len(profiles)}))


if __name__ == '__main__':
    main()
