#!/usr/bin/env python3
"""Transport Ward operators through the exact stabilized local splitting."""

import argparse
import hashlib
import json
from pathlib import Path

from check_modular_crossing import rank_two
from check_pfaffian_charge_gauges import E, N, exterior, kronecker
from complement_shell_tensor import multiply
from pfaffian_particle_hole import local_exterior_table, one_particle_matrix, particle_hole_mask
from stable_pfaffian_bounds import HG_INTERTWINER_COLUMNS


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def apply(columns, vector):
    result = 0
    for j in bits(vector):
        result ^= columns[j]
    return result


def inverse(columns):
    pivots = {}
    for j, column in enumerate(columns):
        value, coordinate = column, 1 << j
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value, coordinate
                break
            old, old_coordinate = pivots[pivot]
            value ^= old
            coordinate ^= old_coordinate
    assert len(pivots) == len(columns)
    result = []
    for j in range(len(columns)):
        value, coordinate = 1 << j, 0
        while value:
            old, old_coordinate = pivots[value.bit_length() - 1]
            value ^= old
            coordinate ^= old_coordinate
        result.append(coordinate)
    return result


def conjugate(operator, columns):
    reverse = inverse(columns)
    return [apply(reverse, apply(operator, column)) for column in columns]


def tensor_sum(first, second):
    a, b = len(first), len(second)
    return [sum(1 << (i * b + y) for i in bits(first[x])) ^
            sum(1 << (x * b + j) for j in bits(second[y]))
            for x in range(a) for y in range(b)]


def tensor_product(first, second):
    size = len(second)
    return [sum(1 << (i * size + j) for i in bits(a) for j in bits(b))
            for a in first for b in second]


def flatten(columns):
    size = len(columns)
    return sum(((column >> i) & 1) << (size * i + j)
               for j, column in enumerate(columns) for i in range(size))


def unflatten(value, size):
    return [sum(((value >> (size * i + j)) & 1) << i for i in range(size)) for j in range(size)]


def block_ranks(columns, intervals):
    return [[rank_two([(columns[j] >> a) & ((1 << (b - a)) - 1) for j in range(c, d)])
             for c, d in intervals] for a, b in intervals]


def local_operators():
    gauge = exterior(kronecker([[1, 0], [N, E]], [[E, N], [0, 1]]), scalar_parity=True)
    incidence = [sum(int(i != j and gauge[i][j] == N) << i for i in range(16)) for j in range(16)]
    creation = [0] * 16
    for y in range(4):
        creation[particle_hole_mask(0, y)] = 1 << particle_hole_mask(3, y)
    one, middle = (1, 2, 4, 8), (3, 5, 9, 6, 10, 12)
    def restrict(operator, indices):
        return [sum(((operator[j] >> i) & 1) << k for k, i in enumerate(indices)) for j in indices]
    return incidence, creation, restrict(incidence, one), restrict(creation, one), restrict(incidence, middle), restrict(creation, middle)


def hg_certificate():
    _, _, rg, og, rh, oh = local_operators()
    columns = list(HG_INTERTWINER_COLUMNS)
    physical_o = conjugate(tensor_sum(oh, [0] * 4), columns)
    ancillary_o = conjugate(tensor_sum([0] * 6, og), columns)
    total_o = [a ^ b for a, b in zip(physical_o, ancillary_o)]
    total_r = conjugate(tensor_sum(rh, rg), columns)
    expected_r = [0] * 24
    expected_o = [0] * 24
    for offset in (0, 4):
        for j in range(4):
            expected_r[offset + j] = rg[j] << offset
            expected_o[offset + j] = og[j] << offset
    for offset in (8, 16):
        for j in range(4):
            for bit in range(2):
                expected_r[offset + 2 * j + bit] = sum(1 << (offset + 2 * i + bit) for i in bits(rg[j]))
                expected_o[offset + 2 * j + bit] = sum(1 << (offset + 2 * i + bit) for i in bits(og[j]))
    for j in range(4):
        expected_o[16 + 2 * j] ^= 1 << (8 + 2 * j)
        expected_o[8 + 2 * j + 1] ^= 1 << (16 + 2 * j + 1)
    assert total_r == expected_r
    assert total_o == expected_o
    intervals = ((0, 4), (4, 8), (8, 16), (16, 24))
    return {'physical_creation_columns': physical_o,
            'ancillary_creation_columns': ancillary_o,
            'total_creation_columns': total_o, 'total_incidence_columns': total_r,
            'physical_creation_block_ranks': block_ranks(physical_o, intervals),
            'total_creation_block_ranks': block_ranks(total_o, intervals),
            'total_incidence_block_ranks': block_ranks(total_r, intervals)}


def stable_columns():
    """G tensor (g+h)^2 basis columns in the S tensor G basis."""
    left = [0] * 64
    for copy, mask in enumerate((0, 15)):
        for j in range(4):
            left[4 * copy + j] = 1 << (4 * mask + j)
    middle = (3, 5, 9, 6, 10, 12)
    for j, column in enumerate(HG_INTERTWINER_COLUMNS):
        target = 8 + j if j < 8 else 40 + j
        left[target] = sum(1 << (4 * middle[i // 4] + i % 4) for i in bits(column))
    for copy, masks in enumerate(((1, 2, 4, 8), (7, 11, 13, 14))):
        for a, mask in enumerate(masks):
            for b in range(4):
                left[16 + 16 * copy + 4 * a + b] = 1 << (4 * mask + b)
    cube = (1, 22, 104, 128, 6, 96, 20, 40)
    right = [0] * 64
    for which in (0, 1):
        for g in range(4):
            natural = (g & 1) if which == 0 else (g >> 1)
            spectator = (g >> 1) if which == 0 else (g & 1)
            for copy in range(2):
                target = (2 * which + copy) * 4 + g
                for value in bits(cube[4 + 2 * copy + natural]):
                    a, b, c = value & 1, value >> 1 & 1, value >> 2 & 1
                    anc = a + 2 * spectator if which == 0 else 2 * a + spectator
                    right[target] ^= 1 << ((anc * 4 + (b + 2 * which)) * 4 + (c + 2 * which))
            for twist in range(2):
                target = 48 + 8 * which + 2 * g + twist
                for value in bits(cube[natural + 2 * twist]):
                    a, b, c = value & 1, value >> 1 & 1, value >> 2 & 1
                    anc = a + 2 * spectator if which == 0 else 2 * a + spectator
                    right[target] ^= 1 << ((anc * 4 + (b + 2 * which)) * 4 + (c + 2 * which))
    for copy in (0, 1):
        for extra in range(4):
            for anc in range(4):
                g, h = extra & 1, extra >> 1
                first, second = (g, h + 2) if copy == 0 else (h + 2, g)
                right[16 + 16 * copy + 4 * extra + anc] = 1 << ((anc * 4 + first) * 4 + second)
    reverse = inverse(right)
    return [apply(left, column) for column in reverse]


def full_stable_certificate():
    columns = stable_columns()
    assert rank_two(columns) == 64
    g_grades = ((1, 0), (0, 1), (2, 1), (1, 2))
    f_grades = ((1, 0), (0, 1), (0, 0), (1, 1))
    source_grades = []
    for mask in range(16):
        first = 1 + (mask >> 2 & 1) - (mask >> 1 & 1)
        second = 1 + (mask >> 3 & 1) - (mask & 1)
        source_grades.extend((first + a, second + b) for a, b in g_grades)
    target_grades = [(a + c + e, b + d + f) for a, b in g_grades
                     for c, d in f_grades for e, f in f_grades]
    for j, column in enumerate(columns):
        assert all(source_grades[i] == target_grades[j] for i in bits(column))
    stride = 8
    p, q, r = 2, 1 << stride, 1 << (stride * stride)
    g, h = ((r, p), (1, 1)), ((1, 1), (q, p ^ q ^ r))
    big = one_particle_matrix(p, q, r)
    f = [[0] * 4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            f[i][j], f[i + 2][j + 2] = g[i][j], h[i][j]
    table = local_exterior_table(p, q, r)
    local = [[0] * 16 for _ in range(16)]
    for a in range(4):
        for y in range(4):
            for v in range(4):
                for u in range(4):
                    local[particle_hole_mask(a, y)][particle_hole_mask(v, u)] = table[((a * 4 + y) * 4 + u) * 4 + v]
    source = [[multiply(local[i // 4][j // 4], big[i % 4][j % 4]) for j in range(64)] for i in range(64)]
    target = [[multiply(big[i // 16][j // 16], multiply(f[i // 4 % 4][j // 4 % 4], f[i % 4][j % 4])) for j in range(64)] for i in range(64)]
    for i in range(64):
        for j in range(64):
            first = second = 0
            for k in bits(columns[j]):
                first ^= source[i][k]
            for k in range(64):
                if columns[k] >> i & 1:
                    second ^= target[k][j]
            assert first == second
    rs, os, rg, og, _, _ = local_operators()
    total_r = conjugate(tensor_sum(rs, rg), columns)
    total_o = conjugate(tensor_sum(os, og), columns)
    rf = (0, 1, 8, 0)
    expected_r = tensor_sum(rg, tensor_sum(rf, rf))
    assert total_r == expected_r
    correction = [a ^ b for a, b in zip(total_o, tensor_sum(og, [0] * 16))]
    ancillary_blocks = [sum(((correction[16 * c + j] >> (16 * a + i)) & 1) << (16 * i + j)
                            for i in range(16) for j in range(16)) for a in range(4) for c in range(4)]
    return {'intertwiner_columns': columns, 'generic_polynomial_equations': 4096,
            'grade_preserving': True,
            'total_incidence_columns': total_r, 'total_creation_columns': total_o,
            'creation_correction_ancillary_Schmidt_rank': rank_two(ancillary_blocks),
            'creation_correction_nonzero_entries': sum(value.bit_count() for value in correction)}


def closed_operator_interface():
    columns = stable_columns()
    reverse = inverse(columns)
    rs, os, _, _, _, _ = local_operators()
    def transport(operator):
        return [apply(reverse, apply(operator, column)) for column in columns]
    def coefficients(operator):
        return [sum(((operator[16 * b + j] >> (16 * a + i)) & 1) << (4 * a + b)
                    for a in range(4) for b in range(4)) for i in range(16) for j in range(16)]
    physical_r = transport(tensor_sum(rs, [0] * 4))
    physical_o = transport(tensor_sum(os, [0] * 4))
    pivots, basis = {}, []
    def add(value):
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot][0]
            else:
                pivots[pivot] = value, 1 << len(basis)
                basis.append(value)
                break
    add(flatten([1 << i for i in range(4)]))
    for operator in (physical_r, physical_o):
        for coefficient in coefficients(operator):
            add(coefficient)
    position = 0
    while position < len(basis):
        operator = transport(tensor_sum([0] * 16, unflatten(basis[position], 4)))
        position += 1
        for coefficient in coefficients(operator):
            add(coefficient)
    assert len(basis) == 8
    def coordinates(value):
        result = 0
        while value:
            pivot = value.bit_length() - 1
            old, coordinate = pivots[pivot]
            value ^= old
            result ^= coordinate
        return result
    def decompose(operator):
        result = [0] * len(basis)
        for index, coefficient in enumerate(coefficients(operator)):
            for which in bits(coordinates(coefficient)):
                result[which] ^= 1 << index
        reconstructed = [0] * 64
        for a, b in zip(basis, result):
            for j, value in enumerate(tensor_product(unflatten(a, 4), unflatten(b, 16))):
                reconstructed[j] ^= value
        assert reconstructed == operator
        return result
    transitions = [decompose(transport(tensor_sum([0] * 16, unflatten(a, 4)))) for a in basis]
    emission_r, emission_o = decompose(physical_r), decompose(physical_o)
    multiplication = []
    for a in basis:
        row = []
        for b in basis:
            value = flatten([apply(unflatten(a, 4), col) for col in unflatten(b, 4)])
            row.append(coordinates(value))
        multiplication.append(row)
    grades = ((1, 0), (0, 1), (2, 1), (1, 2))
    shifts = []
    for value in basis:
        found = {(grades[i][0] - grades[j][0], grades[i][1] - grades[j][1])
                 for i in range(4) for j in range(4) if value >> (4 * i + j) & 1}
        assert len(found) == 1
        shifts.append(list(found.pop()))
    # Two-site direct conjugation checks the recurrence independently.
    doubled = []
    for original in range(1024):
        old, new = original // 16, original % 16
        result = 0
        for row in bits(columns[old]):
            physical, ancillary = row // 4, row % 4
            for last in bits(columns[16 * ancillary + new]):
                result ^= 1 << (64 * physical + last)
        doubled.append(result)
    reverse_doubled = inverse(doubled)
    once = [transport(tensor_sum([0] * 16, unflatten(a, 4))) for a in basis]
    checked = 0
    for local, initial, emission in ((rs, physical_r, emission_r), (os, physical_o, emission_o)):
        expected = tensor_product(initial, [1 << i for i in range(16)])
        for operator, coefficient in zip(once, emission):
            for j, value in enumerate(tensor_product(operator, unflatten(coefficient, 16))):
                expected[j] ^= value
        original = tensor_sum(tensor_sum(local, local), [0] * 4)
        actual = [apply(reverse_doubled, apply(original, column)) for column in doubled]
        assert actual == expected
        checked += len(actual)
    return {'ancillary_algebra_dimension': len(basis), 'homogeneous_basis_matrices': basis,
            'basis_grade_shifts': shifts, 'multiplication_table': multiplication,
            'transport_transitions': transitions, 'incidence_emission': emission_r,
            'creation_emission': emission_o, 'two_site_direct_operator_columns_checked': checked,
            'global_single_insertion_bond_bound': len(basis) + 1}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'middle_stable_operators': hg_certificate(), 'full_stable_operators': full_stable_certificate(),
              'closed_operator_interface': closed_operator_interface()}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
