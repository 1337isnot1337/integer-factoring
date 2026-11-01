#!/usr/bin/env python3
"""Literal rational-gauge replay and derivative checks for shared gauges.

The controls include negative letter potentials. Their ramified Smith
profiles are calculated independently from literal ambient Pfaffians with
integer binary-polynomial arithmetic, then compared with the weighted
compound oracle. The finite-difference check is numerical and separate.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_block_smith_crossings import smith_profile
from check_cancellation_module_automata import pmul
from coherent_graded_quotients import parameter_images, quotient_maps
from coherent_hot_compounds import height, support
from coherent_inverse_systems import inverse_system
from coherent_weighted_gauges import SharedGaugeOracle
from eight_letter_reflection_retraction import halve_exponents, shift_polynomial
from reflection_dual_rail_gaussian import gauged_kernel
from search_coherent_core_products import core_maps
from search_coherent_hot_retractions import _atomic_write


def selected_kernel(modes, place='zero'):
    interval = 4
    while interval < 4 * modes:
        interval *= 2
    original = gauged_kernel(interval, place)
    half = interval // 2
    indices = [half * leg + half - 1 - a for leg in range(4) for a in range(modes)]
    return [[halve_exponents(original[i][j]) for j in indices] for i in indices]


def kernel_controls():
    count = 0
    for modes in range(1, 9):
        zero, infinity = selected_kernel(modes), selected_kernel(modes, 'infinity')
        for i, j in product(range(4 * modes), repeat=2):
            leg_i, a = divmod(i, modes)
            leg_j, b = divmod(j, modes)
            swapped_i = modes * (0, 2, 1, 3)[leg_i] + a
            swapped_j = modes * (0, 2, 1, 3)[leg_j] + b
            assert infinity[i][j] == zero[swapped_i][swapped_j]
            if leg_i == leg_j:
                assert zero[i][j] == 0
            elif {leg_i, leg_j} in ({0, 1}, {2, 3}):
                assert zero[i][j] == 3 << (2 * (a & b))
            count += 1
    return count


def literal_profiles(embedding, retraction, modes, numerators):
    kernel = selected_kernel(modes)
    q = len(embedding)
    ps, ls = tuple(map(support, embedding)), tuple(map(support, retraction))

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return 1
        if mask.bit_count() % 2:
            return 0
        first = mask & -mask
        i = first.bit_length() - 1
        remaining = rest = mask ^ first
        value = 0
        while remaining:
            chosen = remaining & -remaining
            remaining ^= chosen
            j = chosen.bit_length() - 1
            value ^= pmul(kernel[i][j], pfaffian(rest ^ chosen))
        return value

    @lru_cache(None)
    def ramify(poly):
        return sum(1 << (16 * degree) for degree in support(poly))

    matrices = [[[0] * (q * q) for _ in range(q * q)] for _ in range(2)]
    for x, y, u, v in product(range(q), repeat=4):
        values = [0, 0]
        for a, b, c, d in product(ls[x], ls[y], ps[u], ps[v]):
            values[0] ^= pfaffian(a | (b << modes) | (c << (2 * modes)) | (d << (3 * modes)))
            values[1] ^= pfaffian(a | (c << modes) | (b << (2 * modes)) | (d << (3 * modes)))
        shift = sum(numerators[j] for j in (x, y, u, v))
        matrices[0][q * x + y][q * u + v] = shift_polynomial(ramify(values[0]), shift)
        matrices[1][q * x + u][q * y + v] = shift_polynomial(ramify(values[1]), shift)
    heights = [height(support(v)[0]) for v in embedding]
    parities = [(heights[a] + heights[b]) % 2 for a, b in product(range(q), repeat=2)]
    complete = []
    for matrix in matrices:
        assert all(not matrix[i][j] for i, j in product(range(q * q), repeat=2)
                   if parities[i] != parities[j])
        orders = []
        for parity in (0, 1):
            indices = [i for i, value in enumerate(parities) if value == parity]
            block = [[matrix[i][j] for j in indices] for i in indices]
            degree = max(value.bit_length() - 1 for row in block for value in row)
            result, _ = smith_profile(block, degree)
            orders.extend(result)
        complete.append(sorted(orders))
    return complete, pfaffian.cache_info().currsize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    counts = Counter(kernel_entries=kernel_controls())
    inverse = inverse_system([(1 << 6) ^ (1 << 9)], 4)
    quotient = quotient_maps(parameter_images(3, 5, 1), 3)
    cases = [('inverse_q6', inverse['embedding'], inverse['retraction'], 4),
             ('coherent_q8', *core_maps(8, True), 4),
             ('nonlinear_q8', quotient['embedding'], quotient['retraction'], 5)]
    records = []
    random = np.random.default_rng(104729)
    for name, embedding, retraction, modes in cases:
        oracle = SharedGaugeOracle(embedding, retraction, modes)
        numerators = []
        for p, l, h in zip(embedding, retraction, oracle.heights.astype(int)):
            values = {sum(1 if a % 2 == 0 else -1 for a in range(modes) if word >> a & 1)
                      for word in support(p | l)}
            assert len(values) == 1
            numerators.append(int(h + 8 * values.pop()))
        parameters = np.asarray([.5] + [numerators[j] / 32 for j in oracle.free])
        result = oracle.evaluate(parameters, details=True)
        assert min(result['potentials']) < 0
        literal, pfaffians = literal_profiles(embedding, retraction, modes, numerators)
        for exact, key in zip(literal, ('zero_orders', 'infinity_orders')):
            scaled = [32 * value for value in result[key]]
            assert all(value == round(value) for value in scaled)
            assert exact == sorted(map(round, scaled))
            counts['complete_rational_literal_cut_comparisons'] += 1
            counts['complete_smith_orders'] += len(exact)
        counts['literal_pfaffians'] += pfaffians
        counts['literal_tensor_entries'] += 2 * len(embedding) ** 4
        interior = parameters.copy()
        interior[1:] += random.uniform(-.0001, .0001, oracle.q - 1)
        at = oracle.evaluate(interior)
        assert at['minimum_cone_slack'] > 0
        gradient_errors = []
        for _ in range(6):
            direction = random.normal(size=oracle.q)
            direction /= np.linalg.norm(direction)
            step = 1e-6
            plus = oracle.evaluate(interior + step * direction)['rho']
            minus = oracle.evaluate(interior - step * direction)['rho']
            numerical = (plus - minus) / (2 * step)
            expected = float(at['gradient'] @ direction)
            error = abs(numerical - expected)
            assert error < 1e-6, (name, error, numerical, expected)
            gradient_errors.append(error)
            counts['directional_gradient_checks'] += 1
        records.append({'name': name, 'modes': modes, 'embedding': embedding,
                        'retraction': retraction, 'letter_numerators_over_32': numerators,
                        'scale_t': .5, 'cone_rows': len(oracle.cone),
                        'minimum_cone_slack': result['minimum_cone_slack'],
                        'literal_orders_after_ramification_32': literal,
                        'numerical_pressure_rho': result['rho'],
                        'gradient_absolute_errors': gradient_errors})
    names = (Path(__file__).name, 'coherent_weighted_gauges.py', 'coherent_hot_mixed.py',
             'coherent_hot_mixed.cpp', 'coherent_hot_compounds_fast.cpp',
             'coherent_hot_compounds_fast.py', 'coherent_hot_compounds.py',
             'coherent_graded_quotients.py', 'coherent_inverse_systems.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_rectangular_smith_f4.cpp', 'selected_gaussian_coefficients.cpp',
             'reflection_dual_rail_gaussian.py', 'check_block_smith_crossings.py',
             'check_cancellation_module_automata.py', 'eight_letter_reflection_retraction.py')
    report = {'status': 'pass', 'scope': __doc__, 'counts': dict(counts), 'records': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'factoring_exponent_improved': False}
    _atomic_write(args.output, report)
    print(json.dumps({'status': 'pass', 'counts': dict(counts)}, indent=2))


if __name__ == '__main__':
    main()
