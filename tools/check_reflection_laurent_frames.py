#!/usr/bin/env python3
"""Exact paired-frame construction and full-profile precontrols."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random

import numpy as np

from f16_identity_operator_profiles import MUL
from f16_weighted_local_oracle import F16WeightedLocalOracle
from fused_six_vertex_crossings import power
from native_rectangular_smith_f4 import RectangularOracle
from reflection_laurent_physical_frames import catalogue, digest, transformed


def binary_rank(rows):
    pivots = {}
    for value in rows:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def third_cut_rank(seed, exponent):
    q = seed['q_x']
    rows = []
    for x in range(q):
        for v in range(q):
            row = 0
            for y in range(q):
                for u in range(q):
                    a, b = seed['values'][x + q * (y + q * (u + q * v))]
                    assert b == 0
                    row |= ((a >> exponent) & 1) << (q * y + u)
            rows.append(row)
    return binary_rank(rows)


def evaluate(value, point):
    result = 0
    for component, multiplier in zip(value, (1, 2)):
        while component:
            bit = component & -component
            result ^= MUL[multiplier][power(point, bit.bit_length() - 1)]
            component ^= bit
    return result


def evaluated_action(values, q, leg, target, source, coefficient):
    """Independent sequential physical-leg actions over F16."""
    current = list(values)
    for acting_leg, row, column in ((leg, target, source), (3 - leg, source, target)):
        stride = q ** acting_leg
        output = current.copy()
        for index in range(len(current)):
            if (index // stride) % q == row:
                source_index = index + (column - row) * stride
                output[index] ^= MUL[coefficient][current[source_index]]
        current = output
    return current


def cut_matrix(seed, cut):
    q = seed['q_x']
    matrix = [[{} for _ in range(q * q)] for _ in range(q * q)]
    for index, (a, b) in enumerate(seed['values']):
        x, y, u, v = ((index // q ** leg) % q for leg in range(4))
        row, column = (q * x + y, q * u + v) if cut == 0 else (q * x + u, q * y + v)
        entry = {}
        for exponent in range(max(a.bit_length(), b.bit_length())):
            coefficient = ((a >> exponent) & 1) + 2 * ((b >> exponent) & 1)
            if coefficient:
                entry[exponent] = coefficient
        matrix[row][column] = entry
    return matrix


def exact_orders(seed, cut, point, heights, scale):
    q = seed['q_x']
    row_legs, column_legs = ((0, 1), (2, 3)) if cut == 0 else ((0, 2), (1, 3))
    rows = [heights[row_legs[0]][a] + heights[row_legs[1]][b] for a in range(q) for b in range(q)]
    columns = [heights[column_legs[0]][a] + heights[column_legs[1]][b] for a in range(q) for b in range(q)]
    oracle = F16WeightedLocalOracle(cut_matrix(seed, cut), seed['degree'])
    return oracle.profile(rows, columns, point, scale)['orders']


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    records, controls = [], []
    generator = random.Random(312109)
    for interval in (4, 8):
        seed, allowed = catalogue(interval)
        q = seed['q_x']
        fields = ['leg', 'target', 'source', 'place', 'scalar', 'minimum_pole', 'maximum_pole']
        compact = []
        for row in allowed:
            orders = row['pole_orders']
            assert orders == list(range(min(orders), max(orders) + 1))
            compact.append([row[key] for key in fields[:5]] + [min(orders), max(orders)])
        candidates = [row for row in allowed if row['scalar'] == 1 and row['leg'] == 0]
        selected = [candidates[i] for i in sorted(set((0, len(candidates) // 3, 2 * len(candidates) // 3, len(candidates) - 1)))]
        entries, inversion_entries, endpoint_entries = 0, 0, 0
        for row in selected:
            for constant in (False, True):
                exponent = max(row['pole_orders']) * (1 if row['place'] == 'infinity' else -1)
                scalar = 2 if constant else 1
                frame = transformed(seed, row['leg'], row['target'], row['source'], exponent, scalar, constant)
                assert frame['degree'] <= seed['degree']
                restored = transformed(frame, row['leg'], row['target'], row['source'], exponent, scalar, constant)
                assert restored['values'] == seed['values']
                inversion_entries += q ** 4
                endpoint_entries += q ** 4
                # F4 evaluations alone can hide polynomial discrepancies;
                # use a genuine F16 point as well as an F4 point.
                for point in (2, 4):
                    before = [evaluate(value, point) for value in seed['values']]
                    monomial = power(point, exponent % 15)
                    coefficient = MUL[scalar][monomial ^ int(constant)]
                    expected = evaluated_action(before, q, row['leg'], row['target'], row['source'], coefficient)
                    actual = [evaluate(value, point) for value in frame['values']]
                    assert actual == expected
                    entries += q ** 4
                if interval == 4:
                    native = RectangularOracle(frame['values'], q, q)
                    scale = Fraction(5, 4)
                    height_sets = [[[Fraction(0)] + [Fraction(generator.randrange(-4, 5), 3) for _ in range(q - 1)]
                                    for _ in range(4)] for _ in range(2)]
                    params = [float(value) for heights in height_sets for leg_heights in heights for value in leg_heights[1:]] + [float(scale)]
                    orders = np.asarray(native.prefixes(params, 1)['orders']).reshape(2, q * q)
                    maximum_error = 0.0
                    for side, cut, point in ((0, 1, None), (1, 0, (0, 1))):
                        exact = exact_orders(frame, cut, point, height_sets[side], scale)
                        # The native two-end objective uses F(1/t) directly;
                        # the independent oracle reverses homogeneously by D.
                        if point is None:
                            exact = [value - frame['degree'] * scale for value in exact]
                        error = max(abs(float(a) - b) for a, b in zip(exact, orders[side]))
                        assert error < 2e-8
                        maximum_error = max(maximum_error, error)
                    controls.append({'frame': frame['frame'], 'degree': frame['degree'],
                                     'matrix_sha256': digest(frame['values']),
                                     'complete_exact_profile_comparisons': 2,
                                     'maximum_profile_error': maximum_error})
        zero_rank = third_cut_rank(seed, 0)
        assert zero_rank == q * q
        assert not any(row['place'] == 'zero' for row in allowed)
        try:
            transformed(seed, 0, 0, 1, -1)
        except ArithmeticError:
            negative_rejected = True
        else:
            raise AssertionError('An actual pole was incorrectly accepted.')
        records.append({'interval': interval, 'alphabet': q, 'degree': seed['degree'],
                        'seed_sha256': digest(seed['values']), 'catalogue_fields': fields,
                        'catalogue': compact,
                        'admissible_monomial_frames': sum(len(row['pole_orders']) for row in allowed),
                        'zero_third_cut_rank': zero_rank,
                        'top_third_cut_rank': third_cut_rank(seed, seed['degree']),
                        'independent_evaluation_entries': entries,
                        'inversion_entries': inversion_entries, 'endpoint_entries': endpoint_entries,
                        'pole_negative_control_rejected': negative_rejected})
    directory = Path(__file__).parent
    names = (Path(__file__).name, 'reflection_laurent_physical_frames.py',
             'reflection_dual_rail_gaussian.py', 'f16_weighted_local_oracle.py',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    result = {'scope': 'Exact finite admissibility catalogues and constructor/profile controls; no pressure or factoring improvement certified.',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
              'cases': records, 'profile_controls': controls,
              'totals': {'complete_exact_profiles': sum(row['complete_exact_profile_comparisons'] for row in controls),
                         'independent_evaluation_entries': sum(row['independent_evaluation_entries'] for row in records),
                         'inversion_entries': sum(row['inversion_entries'] for row in records),
                         'endpoint_entries': sum(row['endpoint_entries'] for row in records)}}
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps(result['totals']))


if __name__ == '__main__':
    main()
