#!/usr/bin/env python3
"""Exact controls for the conditional three-axis cube crossing transfer.

This checks the geometry and cost interface. It does not construct an
eight-leg family meeting the required width and polynomial-piece bounds.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import comb
from pathlib import Path
from random import Random


def edge(first, second):
    return tuple(sorted((first, second)))


def cube_geometry(g):
    vertices = list(product(range(g), repeat=3))
    generic = {edge(a, b) for a, b in combinations(vertices, 2)
               if all(x != y for x, y in zip(a, b))}
    covered, cubes = {}, []
    for coordinate_pairs in product(combinations(range(g), 2), repeat=3):
        corners = [tuple(coordinate_pairs[axis][word >> axis & 1]
                         for axis in range(3)) for word in range(8)]
        for word in range(4):
            item = edge(corners[word], corners[word ^ 7])
            assert item not in covered
            covered[item] = coordinate_pairs
            assert all(corners[word][axis] != corners[word ^ 7][axis] for axis in range(3))
        for axis in range(3):
            assert Counter(corner[axis] for corner in corners) == Counter({
                coordinate_pairs[axis][0]: 4, coordinate_pairs[axis][1]: 4})
        cubes.append(coordinate_pairs)
    assert set(covered) == generic
    assert len(cubes) == comb(g, 2) ** 3
    assert len(generic) == 4 * len(cubes) == (g * (g - 1)) ** 3 // 2
    residual = [edge(a, b) for a, b in combinations(vertices, 2)
                if edge(a, b) not in generic]
    assert len(residual) == 3 * g ** 4 * (g - 1) // 2
    residual_counts = [Counter() for _ in range(3)]
    for a, b in residual:
        for axis in range(3):
            if a[axis] != b[axis]:
                residual_counts[axis][tuple(sorted((a[axis], b[axis])))] += 1
    for axis in range(3):
        assert len(residual_counts[axis]) == comb(g, 2)
        assert set(residual_counts[axis].values()) == {2 * g ** 3 - g ** 2}
    return cubes, residual_counts, {
        'side': g, 'vertices': len(vertices), 'cubes': len(cubes),
        'generic_edges': len(generic), 'residual_edges': len(residual),
        'axis_residual_edges_per_plane_pair': 2 * g ** 3 - g ** 2,
    }


def width_controls(g, cubes, residual_counts):
    random = Random(31003 + g)
    records = []
    for label_bits in (1, 2, 3):
        volume = g ** 6 * label_bits
        weights = [{key: value * label_bits for key, value in counts.items()}
                   for counts in residual_counts]
        cube_sums = [0, 0, 0]
        for coordinate_pairs in cubes:
            while True:
                widths = [random.randrange(4 * label_bits + 1) for _ in range(3)]
                if sum(widths) <= 9 * label_bits:
                    break
            for axis, width in enumerate(widths):
                weights[axis][coordinate_pairs[axis]] += width
                cube_sums[axis] += width
        selected = min(range(3), key=cube_sums.__getitem__)
        observed = 2 * sum(weights[selected].values())
        bound = (Fraction(3, 4) * (g - 1) ** 3 * g ** 3 * label_bits
                 + g ** 3 * (g - 1) * (2 * g - 1) * label_bits)
        assert observed <= bound
        assert observed <= (Fraction(3, 4) + Fraction(2, g)) * volume
        assert all(max(axis.values()) <= g ** 4 * label_bits for axis in weights)
        assert all(2 * sum(axis.values()) <= Fraction(g - 1, g) * volume for axis in weights)
        records.append({'label_bits': label_bits, 'selected_axis': selected,
                        'oriented_volume': observed, 'volume_bound_numerator': bound.numerator,
                        'volume_bound_denominator': bound.denominator,
                        'maximum_observed_edge_width': max(max(axis.values()) for axis in weights)})
    return records


def cut_weight(matrix, order, m, first=0, second=1):
    return sum(matrix[a][b] for a in order[first * m:(first + 1) * m]
               for b in order[second * m:(second + 1) * m])


def equitable_partition_controls():
    matrices = []
    pairs = list(combinations(range(4), 2))
    for bits in range(1 << len(pairs)):
        matrix = [[0] * 4 for _ in range(4)]
        for bit, (i, j) in enumerate(pairs):
            matrix[i][j] = matrix[j][i] = bits >> bit & 1
        matrices.append((matrix, 2))
    random = Random(731104)
    for _ in range(8):
        matrix = [[0] * 6 for _ in range(6)]
        for i, j in combinations(range(6), 2):
            matrix[i][j] = matrix[j][i] = random.randrange(14)
        matrices.append((matrix, 3))
    means = swaps = 0
    for matrix, groups in matrices:
        g, m = len(matrix), len(matrix) // groups
        total = sum(sum(row) for row in matrix)
        maximum = max(max(row) for row in matrix)
        orders = list(permutations(range(g)))
        for first, second in combinations(range(groups), 2):
            observed = sum(cut_weight(matrix, order, m, first, second) for order in orders)
            assert Fraction(observed, len(orders)) == Fraction(total * g, (g - 1) * groups ** 2)
            means += 1
        for order in orders:
            before = cut_weight(matrix, order, m)
            for a, b in combinations(range(g), 2):
                changed = list(order)
                changed[a], changed[b] = changed[b], changed[a]
                after = cut_weight(matrix, changed, m)
                assert abs(after - before) <= 2 * m * maximum
                swaps += 1
    return {'graphs': len(matrices), 'exact_cut_means': means,
            'permutation_swap_bounds': swaps}


def pure_tensor_negative_control():
    endpoint_checks = factor_checks = 0
    for assignment in range(256):
        labels = [assignment >> word & 1 for word in range(8)]
        endpoint = int(all(labels[word] == labels[word ^ 7] for word in range(4)))
        total = 0
        for piece in range(16):
            wanted = [0] * 8
            for word in range(4):
                wanted[word] = wanted[word ^ 7] = piece >> word & 1
            value = int(labels == wanted)
            total += value
            for axis in range(3):
                left = int(all(labels[word] == wanted[word]
                               for word in range(8) if not (word >> axis & 1)))
                right = int(all(labels[word] == wanted[word]
                                for word in range(8) if word >> axis & 1))
                assert left * right == value
                factor_checks += 1
        assert total == endpoint
        endpoint_checks += 1
    for g in range(2, 10):
        for label_bits in range(1, 6):
            exponent = 4 * label_bits * comb(g, 2) ** 3
            assert exponent == Fraction(1, 2) * (g - 1) ** 3 * g ** 3 * label_bits
    return {'endpoint_entries': endpoint_checks, 'all_cut_factor_entries': factor_checks,
            'pieces_at_one_bit': 16, 'piece_log2_at_lambda_bits': '4*lambda',
            'qualifying_polynomial_piece_family': False}


def schedule_controls():
    checks = balance = 0
    for x0 in (2, 5, 18, 32):
        level_log = x0
        for level in range(13):
            formula = ((5 * x0 + 9) * 3 ** level - 18 * level - 3 * x0 - 9) // 2
            assert level_log == formula
            next_log = 3 * (level_log + 6 * level + x0)
            assert next_log <= 6 * level_log
            if x0 == 32:
                # ln(K) <= log2(K), so this is a rational sufficient form of (8).
                ratio = 1 << (6 * level + x0)
                zeta_inverse_square = 1 << (2 * level + 16)
                assert ratio >= 16 * zeta_inverse_square * (2 * level_log + 1)
                balance += 1
            level_log = next_log
            checks += 1
    return {'exact_logarithm_and_growth_checks': checks,
            'equitable_schedule_inequalities': balance,
            'cutoff_label_lower_exponent': [12, 25],
            'piece_enumeration_upper_exponent': [13, 25]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    geometry = []
    for side in range(2, 7):
        cubes, residuals, record = cube_geometry(side)
        record['width_controls'] = width_controls(side, cubes, residuals)
        geometry.append(record)
    root = Path(__file__).resolve().parent.parent
    paths = [Path(__file__), root / 'research/cube-crossing-transfer.md']
    report = {'status': 'pass', 'scope': __doc__, 'geometry': geometry,
              'equitable_partition': equitable_partition_controls(),
              'pure_tensor_negative_control': pure_tensor_negative_control(),
              'schedule': schedule_controls(),
              'new_local_crossing_family_constructed': False,
              'new_factoring_exponent_established': False,
              'source_sha256': {str(path.relative_to(root)): sha256(path.read_bytes()).hexdigest()
                                for path in paths}}
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
