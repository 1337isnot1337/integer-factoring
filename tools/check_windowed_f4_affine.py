#!/usr/bin/env python3
"""Independent F4 profiles and binary-fiber equivalence for coupled frames."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from f4_lattice_linear import eliminate, encode
from full_alphabet_two_cut_lattices import FullAlphabetTwoCutLattice, cut_location
from graded_adjugate_caps import inequalities
from joint_energy_conjugations import Space
from weighted_smith_oracle_f4 import WeightedSmithOracle
from windowed_f4_affine_lattices import F4AffineLattice, CoupledF4Fiber
from windowed_f4_smith import matrix_product, pack, planes, regular_smith_frames


def independent_profiles(values, model):
    result = []
    for side in range(2):
        matrix = [[0]*model.size for _ in range(model.size)]
        for entry, pair in enumerate(values):
            i, j = cut_location(model.labels[entry], model.q, side)
            matrix[i][j] = pack(pair)
        h = model.heights[4*side:4*side+4]
        rw = [h[0][a]+h[2 if side == 0 else 1][b] for a, b in product(range(model.q), repeat=2)]
        cw = [h[1 if side == 0 else 2][a]+h[3][b] for a, b in product(range(model.q), repeat=2)]
        orders, _ = WeightedSmithOracle(matrix, infinity=side == 0).profile(rw, cw)
        result.append(tuple(orders))
    return tuple(result)


def augmented_rows(model, targets, original_binary=False):
    out = []
    for target, equations in zip(targets, model.equations):
        for power, row, rhs in equations:
            if power >= target:
                break
            if original_binary:
                n = model.variables
                out += [row | (rhs << (2*n)), row << n]
            else:
                out.append(row | (rhs << model.variables))
    return out


def equal_row_spaces(first, second):
    a, b = Space(first), Space(second)
    assert len(a.basis) == len(b.basis)
    for v in first:
        b.binary_coordinates(v)
    for v in second:
        a.binary_coordinates(v)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = Random(0x01A195C2F4)
    binary_controls, compared_systems, endpoint_entries, independent_orders = [], 0, 0, 0
    for trial in range(8):
        q, degree = 2, 3
        values = []
        for index in range(q**4):
            x, y, u, v = ((index//q**j) % q for j in range(4))
            a = rng.getrandbits(degree)
            values.append((int(x == v and y == u) ^ a ^ (a << 1), 0))
        heights = [[rng.randint(-2, 2) for _ in range(q)] for _ in range(8)]
        original = FullAlphabetTwoCutLattice(values, q, degree, heights, raise_cap=3)
        coupled = F4AffineLattice(values, q, degree, heights, raise_cap=3)
        assert original.orders == coupled.orders and original.caps == coupled.caps
        targets = [original.orders, tuple(0 for _ in original.orders)]
        targets += [tuple(rng.randrange(c+1) for c in original.caps) for _ in range(4)]
        for target in targets:
            equal_row_spaces(augmented_rows(original, target, True), augmented_rows(coupled, target))
            old, new = original.solve(target), coupled.solve(target)
            assert (old is None) == (new is None)
            if old is not None:
                assert len(new) == 2*len(old)
                witness = coupled.sample(new, rng)
                profiles = coupled.check_witness(witness, target)
                assert original.check_witness(witness, target) == profiles
                assert independent_profiles(witness, coupled) == profiles
                independent_orders += 2*q*q
                endpoint_entries += q**4
            compared_systems += 1
        binary_controls.append({'trial': trial, 'source': values, 'heights': heights,
                                'orders': original.orders, 'field_coefficient_variables': original.variables,
                                'coupled_binary_variables': coupled.variables})
    nonbinary_controls = []
    for q, trials in ((2, 12), (3, 4)):
        for trial in range(trials):
            degree, values = 3, []
            for index in range(q**4):
                x, y, u, v = ((index//q**j) % q for j in range(4))
                a, b = rng.getrandbits(degree), rng.getrandbits(degree)
                values.append((int(x == v and y == u) ^ a ^ (a << 1), b ^ (b << 1)))
            heights = [[rng.randint(-2, 2) for _ in range(q)] for _ in range(8)]
            model = F4AffineLattice(values, q, degree, heights, raise_cap=3)
            source_profiles = model.profiles(values)
            assert independent_profiles(values, model) == source_profiles
            independent_orders += 2*q*q
            for side in range(2):
                regular = tuple(v+model.scalar_shifts[side] for v in source_profiles[side])
                assert regular == model.orders[side*q*q:(side+1)*q*q]
                constants = [[entry & 3 for entry in row] for row in model.left_frames[side]]
                assert len(eliminate([encode(row) for row in constants])) == q*q
            floors = tuple(max(0, v-1) for v in model.orders)
            fiber = CoupledF4Fiber(model, floors)
            assert fiber.variables % 2 == 0
            full, compressed = model.solve(floors), fiber.solve(floors)
            assert full is not None and compressed is not None
            witness = fiber.sample(compressed, rng)
            actual = model.check_witness(witness, floors)
            assert independent_profiles(witness, model) == actual
            independent_orders += 2*q*q
            endpoint_entries += q**4
            hi, h0 = heights[:4], heights[4:]
            ex = [hi[0][j]+h0[0][j]+hi[3][j]+h0[3][j] for j in range(q)]
            ey = [hi[1][j]+h0[1][j]+hi[2][j]+h0[2][j] for j in range(q)]
            d = min(ex)+min(ey)-actual[0][0]-actual[1][0]
            ex, ey = [v-min(ex) for v in ex], [v-min(ey) for v in ey]
            r = [v-actual[1][0] for v in actual[1]]
            s = [v-actual[0][0] for v in actual[0]]
            assert all(cap['slack'] >= 0 for cap in inequalities(d, r, s, ex, ey))
            nonbinary_controls.append({'q': q, 'trial': trial, 'source': values, 'heights': heights,
                'source_profiles': source_profiles, 'floors': floors, 'binary_floor_fiber_variables': fiber.variables,
                'witness': witness, 'witness_profiles': actual, 'frame_windows': model.frame_windows})
    # Known nonzero orders after dense genuinely F4 polynomial units.
    frame_controls = []
    from weighted_smith_oracle_f4 import polynomial_multiply
    for n in (2, 3, 4):
        unit = [[int(i == j) for j in range(n)] for i in range(n)]
        for _ in range(8):
            i, j = rng.sample(range(n), 2)
            p = rng.getrandbits(8)
            unit[i] = [a ^ polynomial_multiply(p, b) for a, b in zip(unit[i], unit[j])]
        expected = tuple(1+5*i for i in range(n))
        diagonal = [[(1 << (2*expected[i])) if i == j else 0 for j in range(n)] for i in range(n)]
        matrix = matrix_product(matrix_product(unit, diagonal), list(map(list, zip(*unit))))
        result = regular_smith_frames(matrix, 4)
        assert result[0] == expected
        oracle = WeightedSmithOracle(matrix)
        independent, _ = oracle.profile([0]*n, [0]*n)
        assert tuple(independent) == expected
        independent_orders += n
        frame_controls.append({'source': matrix, 'orders': expected, 'window': result[4]})
    retry = regular_smith_frames([[1, 0], [0, 2 << 200]], 4)
    assert retry[0] == (0, 100) and len(retry[4]['attempts']) > 1
    root = Path(__file__).resolve().parents[1]
    names = ['tools/check_windowed_f4_affine.py', 'tools/windowed_f4_smith.py',
             'tools/windowed_f4_affine_lattices.py', 'tools/full_alphabet_two_cut_lattices.py',
             'tools/weighted_smith_oracle_f4.py', 'tools/native_matrix_smith_f4.py',
             'tools/graded_adjugate_caps.py', 'tools/f4_lattice_linear.py']
    report = {'status': 'pass', 'binary_augmented_row_space_comparisons': compared_systems,
              'independent_complete_fraction_free_orders': independent_orders,
              'literal_endpoint_entries': endpoint_entries, 'binary_controls': binary_controls,
              'nonbinary_controls': nonbinary_controls, 'field_frame_controls': frame_controls,
              'forced_field_precision_retry': retry[4],
              'source_sha256': {n: sha256((root/n).read_bytes()).hexdigest() for n in names},
              'scope': 'Coupled F4 equations and genuinely F4 reference frames; the old solver '
              'already allowed F4 coefficients in binary reference frames. No new contraction '
              'or improved factoring exponent is asserted.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'binary_augmented_row_space_comparisons',
                       'independent_complete_fraction_free_orders', 'literal_endpoint_entries')}))


if __name__ == '__main__':
    main()
