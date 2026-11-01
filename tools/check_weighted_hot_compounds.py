#!/usr/bin/env python3
"""Full reflection profiles versus independently assembled weighted compounds."""
from __future__ import annotations
import argparse
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import random

from affine_place_profile_oracle import primitive_seed
from check_reflection_natural_cut_duality import determinant
from hot_code_charge_profiles import hot_energy, natural_profiles, one_particle
from hot_compound_crossing import compound_crossing
from native_matrix_smith_f4 import MatrixOracle
from native_rectangular_smith_f4 import RectangularOracle
from restrict_implicit_reflection import restricted_tensor


def logical(mask, modes):
    return sum(1 << (modes - 1 - bit) for bit in range(modes) if mask >> bit & 1)


def allowed(codes, gauges, stabilized=True, modes=None):
    for first, second in ((0, 1), (2, 3)):
        for i, x in enumerate(codes[first]):
            for j, y in enumerate(codes[second]):
                for a in range(max(codes[first]).bit_length()):
                    if x >> a & 1 or x | 1 << a not in codes[first]:
                        continue
                    dx = gauges[first][codes[first].index(x | 1 << a)] - gauges[first][i]
                    for b in range(max(codes[second]).bit_length()):
                        if y >> b & 1 or y | 1 << b not in codes[second]:
                            continue
                        dy = gauges[second][codes[second].index(y | 1 << b)] - gauges[second][j]
                        if dx + dy + 4 * (a & b) < 0:
                            return False
    if not stabilized:
        for leg, code in enumerate(codes):
            for i, x in enumerate(code):
                for a in range(modes):
                    b = modes - 1 - a
                    if a >= b or x >> a & 1 or x >> b & 1:
                        continue
                    end = x | 1 << a | 1 << b
                    if end in code and gauges[leg][code.index(end)] < gauges[leg][i]:
                        return False
    return True


def block_orders(data, gauges, swap=False):
    if swap:
        gauges = [gauges[0], gauges[2], gauges[1], gauges[3]]
    xs, ys = data['x_letters'], data['y_letters']
    pairs = [(i, j) for i, x in enumerate(xs) for j, y in enumerate(ys)]
    orders = []
    for block in data['blocks']:
        chosen = [(i, j) for i, j in pairs
                  if xs[i].bit_count() + ys[j].bit_count() == block['charge']]
        rows = [gauges[0][i] + gauges[1][j] for i, j in chosen]
        columns = [gauges[2][j] + gauges[3][i] for i, j in chosen]
        matrix = [[(value, 0) for value in row] for row in block['matrix']]
        orders.extend(MatrixOracle(matrix, block['size']).profile(rows, columns, 2)['orders'])
    return sorted(orders)


def monotone(code, rng):
    values = {0: Fraction(0)}
    for mask in sorted(code, key=lambda x: (x.bit_count(), x)):
        if mask:
            prior = [value for smaller, value in values.items() if smaller & ~mask == 0]
            values[mask] = max(prior) + Fraction(rng.randrange(1, 7), 4)
    return [values[mask] for mask in code]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(202611010051)
    cases = ((4, [0, 1, 2, 3], [0, 1, 2, 3]),
             (8, [0, 1, 2, 3], [0, 1, 2, 3]),
             (8, [0, 1, 2], [0, 1]),
             (16, [0, 1, 2, 3, 4, 5], [0, 1, 2, 3]))
    profiles = channels = endpoints = determinant_entries = costs = 0
    records = []
    for interval, xs, ys in cases:
        ambient = interval // 2
        raw = restricted_tensor(interval, [logical(x, ambient) for x in xs],
                                [logical(y, ambient) for y in ys])
        seed, common = primitive_seed(raw)
        assert len(common) == 1 and next(iter(common.values())) == 1
        shift = interval * interval // 2 - next(iter(common))
        data = natural_profiles(xs, ys, keep_blocks=True)
        oracle = RectangularOracle(seed['values'], len(xs), len(ys))
        codes = [xs, ys, ys, xs]
        stabilized = max(map(hot_energy, xs + ys)) < 4 * interval
        for trial in range(10):
            zero = [monotone(code, rng) for code in codes]
            infinity = [monotone(code, rng) for code in codes]
            if stabilized and trial < 3:
                c, d = Fraction(trial + 1, 4), Fraction(3 - trial, 4)
                zero = [[factor * mask.bit_count() for mask in code]
                        for code, factor in zip(codes, (c, -c, d, -d))]
                infinity = [[factor * mask.bit_count() for mask in code]
                            for code, factor in zip(codes, (c, d, -c, -d))]
            assert allowed(codes, zero, stabilized, ambient)
            swapped = [infinity[0], infinity[2], infinity[1], infinity[3]]
            assert allowed(codes, swapped, stabilized, ambient)
            params = []
            for side, gauges in enumerate((infinity, zero)):
                for leg, (code, weights) in enumerate(zip(codes, gauges)):
                    for index in range(1, len(code)):
                        natural = hot_energy(code[index]) // 2 if leg in ((2, 3) if side == 0 else (0, 1)) else 0
                        params.append(float(weights[index] + natural))
            actual = oracle.prefixes(params + [1.])['orders']
            n = len(xs) * len(ys)
            expected = [[value - shift for value in block_orders(data, infinity, True)],
                        [value + shift for value in block_orders(data, zero)]]
            sides = [sorted(actual[:n]), sorted(actual[n:])]
            assert sides == expected
            px = [Fraction(i + 1, len(xs) * (len(xs) + 1) // 2) for i in range(len(xs))]
            py = [Fraction(len(ys) - i, len(ys) * (len(ys) + 1) // 2) for i in range(len(ys))]
            mean = sum(p * hot_energy(x) for p, x in zip(px, xs)) + sum(p * hot_energy(y) for p, y in zip(py, ys))
            for gauges in (zero, infinity):
                for leg, weights in enumerate(gauges):
                    mean += sum(p * g for p, g in zip(px if leg in (0, 3) else py, weights))
            def logsum(side):
                low = min(side)
                return -low + math.log(sum(math.exp(low - v) for v in side))
            predicted_cost = logsum(block_orders(data, infinity, True)) + logsum(block_orders(data, zero)) + float(mean)
            actual_cost = oracle.evaluate(params + [1.], list(map(float, px)), list(map(float, py)))[0]
            assert abs(predicted_cost - actual_cost) < 3e-12
            costs += 1
            profiles += 2
            channels += 2 * n
            records.append({'interval': interval, 'x_letters': xs, 'y_letters': ys,
                            'additional_zero_weights': [[str(v) for v in leg] for leg in zero],
                            'additional_infinity_weights': [[str(v) for v in leg] for leg in infinity],
                            'complete_orders_times_four': [[int(4 * v) for v in side] for side in sides]})
        # Endpoint and complete small minor controls for the separate strict tensor.
        if max(data['modes'], 0) <= 2:
            crossing = compound_crossing(xs, ys)
            matrix = one_particle(data['modes'])
            qx, qy = len(xs), len(ys)
            for index, value in enumerate(crossing['values']):
                x, y = index % qx, index // qx % qy
                u, v = index // (qx * qy) % qy, index // (qx * qy * qy)
                rr = xs[x] | ys[y] << data['modes']
                cc = ys[u] | xs[v] << data['modes']
                rows = [i for i in range(len(matrix)) if rr >> i & 1]
                columns = [i for i in range(len(matrix)) if cc >> i & 1]
                independent = determinant([[matrix[i][j] for j in columns] for i in rows]) if len(rows) == len(columns) else 0
                assert value == (independent, 0)
                assert value[0].bit_count() % 2 == int(x == v and y == u)
                determinant_entries += 1
                endpoints += 1
    # A genuine pole in a triangular factor can change the complete profile.
    xs = ys = [0, 1]
    raw = restricted_tensor(8, [0, 8])
    seed, common = primitive_seed(raw)
    shift = 32 - next(iter(common))
    data = natural_profiles(xs, ys, keep_blocks=True)
    bad = [[Fraction(0), Fraction(-1)], [Fraction(0), Fraction(-1)],
           [Fraction(0), Fraction(0)], [Fraction(0), Fraction(0)]]
    assert not allowed([xs, ys, ys, xs], bad)
    params = [0, 0, 2, 2, 1, 1, 0, 0, 1]
    actual = sorted(RectangularOracle(seed['values'], 2, 2).prefixes(params)['orders'][4:])
    predicted = [v + shift for v in block_orders(data, bad)]
    assert actual != predicted
    names = (Path(__file__).name, 'hot_compound_crossing.py', 'hot_code_charge_profiles.py',
             'hot_code_charge_profiles.cpp', 'restrict_implicit_reflection.py',
             'affine_place_profile_oracle.py', 'check_reflection_natural_cut_duality.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_rectangular_smith_f4.py', 'native_rectangular_smith_f4.cpp')
    report = {'scope': 'Complete weighted profiles and physical means inside a sufficient Gaussian integrality cone.',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'counts': {'complete_profiles': profiles, 'profile_channels': channels,
                         'full_cost_comparisons': costs, 'independent_compound_entries': determinant_entries,
                         'compound_endpoints': endpoints, 'outside_cone_countercontrols': 1},
              'records': records,
              'outside_cone_countercontrol': {'true_orders': actual, 'invalid_shortcut_orders': predicted},
              'passed': True}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'passed': True, 'counts': report['counts']}))


if __name__ == '__main__':
    main()
