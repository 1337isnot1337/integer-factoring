#!/usr/bin/env python3
"""Independent regular-support and original-tensor checks for hot-core defects."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_rectangular_defect_support import ZERO, ONE, Z, apply_shear, pairing
from check_rectangular_f4_repair import add as field_add, scale
from native_matrix_smith_f4 import MatrixOracle
from reflection_dual_rail_gaussian import tensor, logical_heights
from relabelled_gaussian_sums import flatten, substitute_power
import reflection_hot_core_profiles as single
import reflection_hot_block_profiles as block


def full_images(interval, ramification, core):
    m, p, k = interval // 2, ramification, core
    c, s, u, large = ONE + Z ** (2 * p), Z ** (4 * p), Z ** (2 * p), Z ** ((2 * interval + 2) * p)
    answer = {}
    def vector():
        return [ZERO] * (2 * m)
    for j in range(k):
        creation, annihilation = vector(), vector()
        creation[m + j] = ONE
        answer['left_first', j] = (creation, annihilation)
        creation, annihilation = vector(), vector()
        creation[j] = ONE
        for r in range(m):
            creation[m + r] = c * s ** (r - (j & r))
        answer['left_second', j] = (creation, annihilation)
        creation, annihilation = vector(), vector()
        annihilation[j] = ONE / (u * s ** j)
        for r in range(m):
            annihilation[m + r] = c * s ** (-(j & r)) / large
        answer['right_first', j] = (creation, annihilation)
        creation, annihilation = vector(), vector()
        annihilation[m + j] = ONE / (large * s ** j)
        answer['right_second', j] = (creation, annihilation)
    return answer


def support_control(interval, ramification, core):
    m, p, k = interval // 2, ramification, core
    before = full_images(interval, ramification, core)
    images = {name: (list(first), list(second)) for name, (first, second) in before.items()}
    orders = [p * (2 + 4 * i) for i in range(interval)]
    trace = []
    groups = m // k
    s = Z ** (4 * p * k)
    geometric = sum((s ** h for h in range(groups - 1)), ZERO)
    if groups >= 2:
        for t in range(k):
            pivot = m + t + k
            for q in range(2, groups):
                apply_shear(images, m + t + k * q, pivot, s ** (q - 1), orders, trace)
            for q in range(2, groups):
                apply_shear(images, pivot, m + t + k * q, ONE / geometric, orders, trace)
    expected, indices, active_orders = block.active_images(interval, ramification, core)
    assert tuple(orders[i] for i in indices) == active_orders
    for name, (creation, annihilation) in images.items():
        assert all(not creation[i] and not annihilation[i] for i in range(interval) if i not in indices)
        assert [creation[i] for i in indices] == expected[name][0]
        assert [annihilation[i] for i in indices] == expected[name][1]
    for first, second in product(before, repeat=2):
        assert pairing(before[first], before[second]) == pairing(images[first], images[second])
    # Pull input operators back through the diagonal before checking integrality.
    for (side, j), (creation, annihilation) in images.items():
        if side.startswith('right'):
            creation = [value * Z ** (-order) for value, order in zip(creation, orders)]
            annihilation = [value * Z ** order for value, order in zip(annihilation, orders)]
        assert all(value.order() >= 0 for value in creation + annihilation if value)
    return {'interval': interval, 'ramification': ramification, 'core': core,
            'active_modes': len(indices), 'spectator_modes': interval - len(indices),
            'regular_shears': len(trace), 'all_pairings': len(images) ** 2}


def polynomial_seed(interval, ramification, factor, core):
    factor = block.normalize_factor(factor, core)
    clearing = max(0, -min(degree for scalar in factor.values() for degree in scalar))
    base = tensor(interval)['values']
    m = interval // 2
    source = [substitute_power(value, ramification) for value, high in base]
    assert all(not high for value, high in base)
    answer = [(0, 0)] * len(base)
    for mask, scalar in factor.items():
        physical = sum(1 << (leg * m + m - 1 - j)
                       for leg in range(4) for j in range(core) if mask >> (leg * core + j) & 1)
        for original, value in enumerate(source):
            if original & physical or not value:
                continue
            target = original | physical
            for degree, coefficient in scalar.items():
                answer[target] = field_add(answer[target], scale((value << (degree + clearing), 0), coefficient))
    for actual, expected in zip(answer, base):
        endpoint = (actual[0].bit_count() & 1) + 2 * (actual[1].bit_count() & 1)
        assert endpoint == (expected[0].bit_count() & 1)
    bound = 2 * ramification * interval * interval + clearing + max(degree for scalar in factor.values() for degree in scalar)
    assert max(max(value).bit_length() - 1 for value in answer) <= bound
    return answer, clearing


def original_profiles(interval, ramification, values):
    m, q = interval // 2, 1 << (interval // 2)
    bits = max(1, max(value.bit_length() for pair in values for value in pair))
    precision = max(512, 1 << (bits - 1).bit_length())
    answer = []
    for side, place in enumerate(('zero', 'infinity')):
        h = logical_heights(interval, place)
        heights = [[ramification * sum(h[leg * m + j] for j in range(m) if state >> j & 1)
                    for state in range(q)] for leg in range(4)]
        first, second = (((0, 1), (2, 3)) if side == 0 else ((0, 2), (1, 3)))
        rows = [heights[first[0]][a] + heights[first[1]][b] for a, b in product(range(q), repeat=2)]
        columns = [heights[second[0]][a] + heights[second[1]][b] for a, b in product(range(q), repeat=2)]
        answer.append(MatrixOracle(flatten(values, m, bool(side)), q * q, precision=precision).profile(
            rows, columns, infinity=bool(side))['orders'])
    return answer


def literal_control(interval, ramification, core, factor):
    values, clearing = polynomial_seed(interval, ramification, factor, core)
    expected = original_profiles(interval, ramification, values)
    data = block.profile_data(interval, ramification, factor, core)
    assert clearing == data['common_laurent_clearing']
    for side in range(2):
        actual = np.asarray(data['profiles'][side]) + data['scalar_orders'][side]
        for spectator in data['spectator_orders']:
            actual = np.r_[actual, actual + spectator]
        assert np.array_equal(np.sort(actual), np.sort(expected[side]))
    return {'interval': interval, 'ramification': ramification, 'core': core,
            'masks': sorted(factor), 'clearing': clearing, 'complete_cut_order': len(expected[0]),
            'active_cut_order': data['active_matrix_order'], 'endpoint_coefficients': len(values)}


def verify(extended=False):
    supports = [support_control(M, P, k) for M in (4, 8, 16, 32, 64)
                for P in (1, 3) for k in (1, 2, 4) if k <= M // 2]
    literal = []
    for M, P in ((4, 1), (8, 1), (4, 3)):
        for mask in single.EVEN_MASKS:
            factor = {0: {0: 1}, mask: {-1: 2, 1: 2}} if mask else {0: {0: 1, 2: 2, 3: 2}}
            literal.append(literal_control(M, P, 1, factor))
            data = single.profile_data(M, P, factor)
            general = block.profile_data(M, P, factor, 1)
            assert data['profiles'] == general['profiles']
    for mask in (3, 5, 15, 17, 85, 170, 240, 255):
        literal.append(literal_control(8, 1, 2, {0: {0: 1}, mask: {-2: 3, 1: 3}}))
    literal.append(literal_control(8, 3, 2,
        {0: {0: 1}, 5: {0: 2, 1: 2}, 85: {-1: 3, 2: 3}, 255: {2: 1, 5: 1}}))
    if extended:
        for mask in block.even_masks(2):
            factor = {0: {0: 1}, mask: {0: 2, 1: 2}} if mask else {0: {0: 1}}
            literal.append(literal_control(4, 1, 2, factor))
    backgrounds = []
    for M in (4, 8, 16, 32, 64):
        data = single.profile_data(M, 1, {0: {0: 1}})
        expected = [0]
        for order in data['active_one_particle_orders']:
            expected += [value + order for value in expected]
        assert tuple(sorted(expected)) == data['profiles'][0] == data['profiles'][1]
        backgrounds.append({'interval': M, **single.optimize_scale(data)})
    root = Path(__file__).resolve().parent.parent
    files = ['tools/reflection_hot_core_profiles.py', 'tools/reflection_hot_block_profiles.py',
             'tools/check_reflection_hot_core_profiles.py', 'research/reflection-hot-core-defects.md']
    return {'status': 'pass', 'uniform_regular_support_controls': supports,
            'original_polynomial_controls': literal,
            'original_complete_cut_profiles': 2 * len(literal),
            'endpoint_coefficient_equalities': sum(record['endpoint_coefficients'] for record in literal),
            'unperturbed_floating_controls': backgrounds, 'extended': extended,
            'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest() for name in files},
            'scope': 'Exact polynomial and profile controls for a proved bounded-core family; no improved exponent claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify(args.extended)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'],
        'regular_support_cases': len(result['uniform_regular_support_controls']),
        'original_complete_cut_profiles': result['original_complete_cut_profiles'],
        'endpoint_coefficient_equalities': result['endpoint_coefficient_equalities']}))


if __name__ == '__main__':
    main()
