#!/usr/bin/env python3
"""Exact endpoint, profile, phase-projection, and Laurent-clearing controls."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

from check_block_smith_crossings import smith_profile
from check_localized_quartic_profiles import gauged_matrix
from check_rectangular_f4_repair import ZERO, add, multiply, scale
from native_matrix_smith_f4 import MatrixOracle
from relabelled_gaussian_sums import (
    RelabelingPressure, affine_relabeling, charge_projection, check_endpoint,
    flatten, gaussian_tensor, natural_heights,
)
from tie_preserving_crossing_exponents import construct


def subset_sums(energies):
    values = [0]
    for energy in energies:
        values += [value + energy for value in values]
    return sorted(values)


def natural_weights(count, place):
    q = 1 << count
    h = natural_heights(tuple(range(count)), count)[place]
    first, second = ((0, 1), (2, 3)) if place == 0 else ((0, 2), (1, 3))
    rows = [h[first[0]][a] + h[first[1]][b] for a, b in product(range(q), repeat=2)]
    columns = [h[second[0]][a] + h[second[1]][b] for a, b in product(range(q), repeat=2)]
    return rows, columns


def phase_check(count, charges):
    q = 1 << count
    base = gaussian_tensor(tuple(range(count)), count)
    projected = charge_projection(base, charges)
    check_endpoint(projected, count)
    powers = [1, 2, 3]
    for mask, original in enumerate(base):
        x, y, u, v = ((mask >> (count * leg)) & (q - 1) for leg in range(4))
        difference = (charges[x] + charges[u] - charges[y] - charges[v]) % 3
        literal = ZERO
        for j in range(3):
            literal = add(literal, scale((original, 0), powers[(j * difference) % 3]))
        assert literal == projected[mask]
    block_records = []
    for place in range(2):
        matrix = flatten(projected, count, bool(place))
        row_weights, column_weights = natural_weights(count, place)
        if place == 0:
            row_charge = [(charges[a] - charges[b]) % 3 for a, b in product(range(q), repeat=2)]
            column_charge = [(charges[b] - charges[a]) % 3 for a, b in product(range(q), repeat=2)]
        else:
            row_charge = column_charge = [(charges[a] + charges[b]) % 3
                                          for a, b in product(range(q), repeat=2)]
        whole = MatrixOracle(matrix, q*q).profile(row_weights, column_weights,
                                                  infinity=bool(place))['orders']
        pieces = []
        sizes = []
        for charge in range(3):
            rows = [i for i, value in enumerate(row_charge) if value == charge]
            columns = [j for j, value in enumerate(column_charge) if value == charge]
            assert len(rows) == len(columns)
            sizes.append(len(rows))
            for i in rows:
                for j, value in enumerate(matrix[i]):
                    if column_charge[j] != charge:
                        assert value == ZERO
            if rows:
                block = [[matrix[i][j] for j in columns] for i in rows]
                pieces += MatrixOracle(block, len(rows)).profile(
                    [row_weights[i] for i in rows], [column_weights[j] for j in columns],
                    infinity=bool(place))['orders']
        assert sorted(pieces) == sorted(whole)
        block_records.append({'place': place, 'block_sizes': sizes,
                              'complete_orders': sorted(int(value) for value in whole)})
    return {'modes': count, 'charges': charges, 'literal_phase_entries': len(base),
            'blocks': block_records}


def verify():
    endpoint_cases = []
    for count, ambient in ((1, 1), (2, 2), (3, 4), (4, 4)):
        values = gaussian_tensor(tuple(range(count)), ambient)
        checks = check_endpoint(tuple((value, 0) for value in values), count)
        endpoint_cases.append({'modes': count, 'ambient': ambient, 'entries': checks})
    profile_cases = []
    for count in (1, 2, 4):
        strict = gaussian_tensor(tuple(range(count)), count)
        raw = gaussian_tensor(tuple(range(count)), count, False)
        for place in range(2):
            row_weights, column_weights = natural_weights(count, place)
            orders = []
            for values in (raw, strict):
                matrix = flatten(tuple((value, 0) for value in values), count, bool(place))
                result = MatrixOracle(matrix, 1 << (2 * count)).profile(
                    row_weights, column_weights, infinity=bool(place))['orders']
                orders.append(sorted(int(value) for value in result))
            energies = ([j for j in range(1, 2 * count, 2) for _ in range(2)] if place == 0
                        else list(range(1 - 2 * count, 2 * count, 2)))
            assert orders[0] == orders[1] == subset_sums(energies)
            independent = False
            if count <= 2:
                matrix, minimum = gauged_matrix(strict, count, bool(place))
                degree = max(value.bit_length() - 1 for row in matrix for value in row)
                python_orders, _ = smith_profile(matrix, degree)
                assert sorted(value + minimum for value in python_orders) == orders[1]
                independent = True
            profile_cases.append({'modes': count, 'place': place, 'orders': orders[1],
                                  'independent_binary_elimination': independent})

    base = gaussian_tensor((0, 1, 2, 3), 4)
    permutation = list(range(16))
    permutation[6], permutation[9] = permutation[9], permutation[6]
    mixture = affine_relabeling(base, permutation, scalar=2)
    check_endpoint(mixture, 4)
    four = [1, 2, 4 + 2, 12 + 1]
    mask = sum(1 << bit for bit in four)
    residual = multiply(mixture[0], mixture[mask])
    for pairs in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        left, right = [sum(1 << four[i] for i in pair) for pair in pairs]
        residual = add(residual, multiply(mixture[left], mixture[right]))
    assert residual == (0, (1 << 10) | (1 << 11))

    phase_records = []
    for charges in ([0, 1, 1, 2], [0, 1, 1, 0], [0, 0, 0, 1]):
        phase_records.append(phase_check(2, list(charges)))
    phase_records.append(phase_check(4, [int(label == 12) for label in range(16)]))

    small = gaussian_tensor((0, 1), 2)
    values = tuple((value, 0) for value in small)
    clearing_cases = []
    for kind in ('common', 'quadratic', 'charge'):
        shifts = []
        for mask in range(len(small)):
            x, y, u, v = ((mask >> (2 * leg)) & 3 for leg in range(4))
            shifts.append(-7 if kind == 'common' else
                          ((x & 1) - (v & 1)) * ((u >> 1) - (y >> 1)) if kind == 'quadratic' else
                          -2 * int((int(x == 3) + int(u == 3) - int(y == 3) - int(v == 3)) % 3 != 0))
        model = RelabelingPressure(values, (0, 1), 2, entry_shifts=shifts)
        cleared = tuple((value << (shift + model.clearing_order), 0)
                        if shift + model.clearing_order >= 0 else
                        (value >> (-shift - model.clearing_order), 0)
                        for value, shift in zip(small, shifts))
        for original, shift, actual in zip(small, shifts, cleared):
            if shift + model.clearing_order < 0:
                assert original % (1 << (-shift - model.clearing_order)) == 0
        check_endpoint(cleared, 2)
        for place in range(2):
            rw, cw = natural_weights(2, place)
            laurent = model.oracles[place].profile(rw, cw, infinity=bool(place))['orders']
            polynomial = MatrixOracle(flatten(cleared, 2, bool(place)), 16).profile(
                rw, cw, infinity=bool(place))['orders']
            shift = model.clearing_order * (-1 if place else 1)
            assert sorted(polynomial) == sorted(value + shift for value in laurent)
        clearing_cases.append({'kind': kind, 'clearing_order': model.clearing_order,
                               'polynomial_degree': model.maximum_degree})

    ties = construct(2)
    assert len(ties['support']) == 62
    assert ties['constraints_by_place'] == [229, 229]
    assert (ties['constraint_rank'], ties['nullity'], ties['physical_leg_gauge_rank']) == (49, 13, 12)
    root = Path(__file__).resolve().parent.parent
    sources = ['tools/check_strict_gaussian_crossings.py', 'tools/relabelled_gaussian_sums.py',
               'tools/tie_preserving_crossing_exponents.py', 'research/strict-gaussian-relabeling-crossings.md']
    return {'status': 'pass', 'endpoint_cases': endpoint_cases,
            'inherited_profile_cases': profile_cases, 'quartic_wick_residual': {'low_bits': '0x0', 'high_bits': hex(residual[1])},
            'charge_projection_cases': phase_records, 'Laurent_clearing_cases': clearing_cases,
            'tie_preserving_directions': ties,
            'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest() for name in sources},
            'scope': 'Construction and finite exact profile controls. No improved pressure or factoring score is asserted.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'],
                      'strict_endpoint_entries': sum(record['entries'] for record in result['endpoint_cases']),
                      'inherited_profile_cases': len(result['inherited_profile_cases']),
                      'charge_projection_cases': len(result['charge_projection_cases']),
                      'Laurent_clearing_cases': len(result['Laurent_clearing_cases']),
                      'nonseparable_tie_directions': len(result['tie_preserving_directions']['additional_integer_directions'])}, indent=2))


if __name__ == '__main__':
    main()
