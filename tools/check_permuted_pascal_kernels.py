#!/usr/bin/env python3
"""Exact Hadamard bound, all small permutations, and actual Pascal cuts."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import random

from check_weighted_pascal_flags import exact_weighted_profile, original_fock_profile, subset_sums
from filtered_pascal_exponents import upset_orders
from permuted_pascal_kernels import (
    explicit_infinity_gauge, permutation_upset_orders, permuted_pascal_kernel,
)
from tensor_pascal_box_crossings import exterior_tensor, matrix_product, polynomial_power
from weighted_pascal_flags import clifford_data, inverse_kernel, natural_heights


def order_at_one(polynomial, length):
    if not polynomial:
        return None
    translated = 0
    for j in range(length):
        if polynomial >> j & 1:
            translated ^= polynomial_power(3, j)
    return (translated & -translated).bit_length() - 1


def hadamard_controls():
    result = []
    for length in (2, 4, 8):
        orders = [order_at_one(value, length) for value in range(1 << length)]
        count = 0
        for first in range(1 << length):
            for second in range(1 << length):
                product = first & second
                if product and first and second:
                    assert orders[product] >= orders[first] + orders[second] - (length - 1)
                count += 1
        result.append({'length': length, 'all_binary_polynomial_pairs': count})
    return result


def permutation_controls():
    counts = {'all_permutations': 0, 'inequalities_checked': 0,
              'nonmonotone_order_functions': 0}
    for images in permutations(range(8)):
        orders = permutation_upset_orders(images)
        assert orders == upset_orders(images)
        assert orders[0] == 7 and orders[7] == 0
        for i in range(8):
            for j in range(8):
                assert orders[i] + orders[j] - orders[i | j] <= 7
                counts['inequalities_checked'] += 1
        if any(orders[i] < orders[j] for i in range(8) for j in range(8) if i & j == i):
            counts['nonmonotone_order_functions'] += 1
        counts['all_permutations'] += 1
    assert counts['all_permutations'] == 40320
    assert counts['nonmonotone_order_functions'] == 31104
    return counts


def check_kernel(images, beta, full_fock=False):
    seed = permuted_pascal_kernel(images)
    length, matrix, pascal = seed['n'], seed['matrix'], seed['pascal']
    first = [row[:length] for row in matrix[:length]]
    second = [row[length:] for row in matrix[length:]]
    expected = [[polynomial_power(3, 2 * length) if value else 0 for value in row] for row in pascal]
    assert matrix_product(first, second) == expected
    assert matrix_product(matrix, inverse_kernel(seed)) == [
        [int(i == j) for j in range(2 * length)] for i in range(2 * length)]
    for i, row in enumerate(second):
        for j, value in enumerate(row):
            actual_order = (value & -value).bit_length() - 1
            assert actual_order == 2 * seed['upset_orders'][i | j]
            assert actual_order + seed['zero_energies'][i] + seed['zero_energies'][j] >= 2 * length
    zero_heights = natural_heights(seed)[0]
    values, shifts, row_indices, column_indices = clifford_data(seed, 0)
    rows = [zero_heights[i] for i in row_indices]
    columns = [zero_heights[i] for i in column_indices]
    zero_orders = exact_weighted_profile(values, shifts, rows + [-x for x in rows],
                                         columns + [-x for x in columns])
    assert zero_orders == sorted((seed['zero_energies'] + [-x for x in seed['zero_energies']]) * 2)
    gauge = explicit_infinity_gauge(seed, beta)
    h = gauge['heights']
    infinity = exact_weighted_profile(matrix, [[0] * (2 * length) for _ in matrix],
                                     h[2 * length:], h[:2 * length], True)
    assert infinity == gauge['one_particle_orders']
    record = {'length': length, 'images': list(images), 'upset_orders': seed['upset_orders'],
              'zero_energies': seed['zero_energies'], 'zero_clifford_orders': zero_orders,
              'infinity_gauge': gauge, 'original_fock_orders_checked': 0}
    if full_fock:
        tensor = exterior_tensor(seed)
        actual_zero = original_fock_profile(seed, tensor, zero_heights, 0)
        actual_infinity = original_fock_profile(seed, tensor, h, 1)
        assert actual_zero == subset_sums(seed['zero_energies'] * 2)
        assert actual_infinity == subset_sums(infinity)
        record['original_fock_orders_checked'] = len(actual_zero) + len(actual_infinity)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hadamard = hadamard_controls()
    permutations_checked = permutation_controls()
    samples = [([0, 1], [1, 3], True), ([1, 0], [2, 7], True)]
    samples += [(list(images), [1, 3, 5, 7], True) for images in
                ((0, 1, 2, 3), (0, 1, 3, 2), (0, 3, 1, 2), (2, 0, 3, 1))]
    samples += [([0, 1, 2, 3, 4, 6, 5, 7], list(range(1, 16, 2)), False),
                ([0, 1, 2, 3, 4, 5, 7, 6], [1, 3, 5, 7, 11, 13, 15, 17], False)]
    rng = random.Random(202611010552)
    for length in (8, 16, 16, 32):
        images = list(range(length))
        rng.shuffle(images)
        weights = [rng.randrange(1, 2 * length) for _ in range(length.bit_length() - 1)]
        beta = [sum(weight for j, weight in enumerate(weights) if i >> j & 1) - length
                for i in range(length)]
        samples.append((images, beta, False))
    kernels = [check_kernel(images, beta, fock) for images, beta, fock in samples]
    assert kernels[6]['upset_orders'] == [7, 1, 1, 4, 3, 1, 2, 0]
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__), Path(__file__).with_name('permuted_pascal_kernels.py'),
               Path(__file__).with_name('filtered_pascal_exponents.py'),
               Path(__file__).with_name('weighted_pascal_flags.py'),
               Path(__file__).with_name('check_weighted_pascal_flags.py'),
               Path(__file__).with_name('tensor_pascal_box_crossings.py'),
               Path(__file__).with_name('cube_symmetric_exterior_blocks.py'),
               root / 'research/permuted-pascal-zero-spectra.md']
    report = {'status': 'pass', 'hadamard_controls': hadamard,
              'permutation_controls': permutations_checked, 'kernel_cases': kernels,
              'factoring_score': '1/2',
              'source_sha256': {str(path.relative_to(root)): sha256(path.read_bytes()).hexdigest()
                                for path in sources}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass', 'hadamard_pairs': sum(
        row['all_binary_polynomial_pairs'] for row in hadamard),
        'permutations': permutations_checked['all_permutations'],
        'kernel_cases': len(kernels), 'original_fock_orders_checked': sum(
            row['original_fock_orders_checked'] for row in kernels)}, indent=2))


if __name__ == '__main__':
    main()
