#!/usr/bin/env python3
"""Independent exact Clifford, scalar, and complete Fock Smith controls."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from check_rectangular_defect_support import Rat, ZERO, ONE, Z, identity, mm
from reflection_dual_rail_gaussian import laurent_kernel, logical_heights, tensor
from reflection_unrestricted_gaussian import bogoliubov_data, fock_spectrum
from relabelled_gaussian_sums import flatten
from native_matrix_smith_f4 import MatrixOracle


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def add(first, second):
    return [[x+y for x, y in zip(a, b)] for a, b in zip(first, second)]


def inverse(matrix):
    size = len(matrix)
    work = [list(row)+unit for row, unit in zip(matrix, identity(size))]
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        value = work[column][column]
        work[column] = [entry/value for entry in work[column]]
        for row in range(size):
            if row != column and work[row][column]:
                value = work[row][column]
                work[row] = [entry+value*other for entry, other in zip(work[row], work[column])]
    return [row[size:] for row in work]


def rational_smith(matrix):
    """Exact minimum-valuation DVR elimination, with no series window."""
    size = len(matrix)
    work = [list(row) for row in matrix]
    orders = []
    for position in range(size):
        order, row, column = min((work[i][j].order(), i, j)
                                 for i in range(position, size)
                                 for j in range(position, size) if work[i][j])
        work[position], work[row] = work[row], work[position]
        for line in work:
            line[position], line[column] = line[column], line[position]
        pivot = work[position][position]
        for i in range(position+1, size):
            if work[i][position]:
                coefficient = work[i][position]/pivot
                assert coefficient.order() >= 0
                for j in range(position+1, size):
                    work[i][j] = work[i][j]+coefficient*work[position][j]
                work[i][position] = ZERO
        # Clearing the pivot row would change no remaining submatrix
        # entry, since the pivot column below it is already zero.
        assert all(not work[position][j] or work[position][j].order() >= order
                   for j in range(position+1, size))
        orders.append(order)
    assert orders == sorted(orders)
    return orders


def from_polynomial(polynomial, shift=0, infinity=False):
    if not infinity:
        return Rat(polynomial)*Z**shift
    value = ZERO
    while polynomial:
        bit = polynomial & -polynomial
        polynomial ^= bit
        value = value+Z**(-shift-(bit.bit_length()-1))
    return value


def closed_matrix(interval, side, infinity=False):
    matrix, shifts, rows, columns = bogoliubov_data(interval, side)
    assert all(not high for row in matrix for low, high in row)
    return [[from_polynomial(value[0], shift, infinity)
             for value, shift in zip(line, offsets)]
            for line, offsets in zip(matrix, shifts)], rows, columns


def kernel_matrix(interval, side):
    m = interval//2
    raw = laurent_kernel(interval)
    legs = (0, 1) if side == 0 else (0, 2)
    rows = [leg*m+i for leg in legs for i in range(m)]
    columns = [i for i in range(2*interval) if i not in rows]
    def block(first, second):
        return [[sum((Z**power for power in raw[i][j]), ZERO) for j in second] for i in first]
    p, k, q = block(rows, rows), block(rows, columns), block(columns, columns)
    kt = inverse(transpose(k))
    pk, kq = mm(p, kt), mm(kt, q)
    upper = add(k, mm(pk, q))
    return [upper[i]+pk[i] for i in range(interval)] + [kq[i]+kt[i] for i in range(interval)]


def quadratic_control(matrix):
    modes = len(matrix)//2
    columns = transpose(matrix)
    for column in columns:
        assert sum((column[i]*column[modes+i] for i in range(modes)), ZERO) == ZERO
    for i, first in enumerate(columns):
        for j, second in enumerate(columns):
            value = sum((first[k]*second[modes+k]+first[modes+k]*second[k]
                         for k in range(modes)), ZERO)
            assert value == int((i+modes) % (2*modes) == j)


def subset_prediction(signed, mean):
    modes = len(signed)//2
    assert signed == [-value for value in reversed(signed)]
    positive = signed[modes:]
    minimum = mean-Fraction(sum(positive), 2)
    assert minimum.denominator == 1
    result = [int(minimum)]
    for order in positive:
        result += [value+order for value in result]
    return sorted(result)


def mode_operators(modes):
    size = 1 << modes
    creations, annihilations = [], []
    for mode in range(modes):
        create = [[ZERO]*size for _ in range(size)]
        annihilate = [[ZERO]*size for _ in range(size)]
        for state in range(size):
            if state >> mode & 1:
                annihilate[state ^ (1 << mode)][state] = ONE
            else:
                create[state | (1 << mode)][state] = ONE
        creations.append(create)
        annihilations.append(annihilate)
    return creations, annihilations


def conjugation_matrix(fock, modes):
    creations, annihilations = mode_operators(modes)
    operators = creations+annihilations
    fi = inverse(fock)
    columns = []
    for operator in operators:
        image = mm(mm(fock, operator), fi)
        coefficients = [image[1 << i][0] for i in range(modes)] + [image[0][1 << i] for i in range(modes)]
        reconstructed = [[sum((coefficient*op[i][j] for coefficient, op in zip(coefficients, operators)), ZERO)
                          for j in range(1 << modes)] for i in range(1 << modes)]
        assert reconstructed == image
        columns.append(coefficients)
    return transpose(columns)


def scaled(matrix, scalar):
    return [[value*scalar for value in row] for row in matrix]


def toy_controls():
    cases = []
    # Half-integral mean with an odd Cartan order; the minimum is integral.
    cases.append(('odd determinant order', [[ONE, ZERO], [ZERO, Z**3]], 1, Fraction(3, 2)))
    c, a = mode_operators(2)
    cases.append(('creation pole', add(identity(4), scaled(mm(c[0], c[1]), Z**-3)), 2, Fraction(0)))
    cases.append(('particle-hole exchange', add(c[0], a[0]), 2, Fraction(0)))
    # v=e_2+f_2 has q(v)=1; the e_1 term in the image of f_1 is required.
    eichler = add(identity(4), mm(c[0], add(c[1], a[1])))
    cases.append(('quadratic Eichler correction', eichler, 2, Fraction(0)))
    records = []
    for name, fock, modes, mean in cases:
        b = conjugation_matrix(fock, modes)
        quadratic_control(b)
        signed = rational_smith(b)
        actual = rational_smith(fock)
        assert sum(actual) == mean*(1 << modes)
        assert actual == subset_prediction(signed, mean)
        if name == 'quadratic Eichler correction':
            assert b[0][modes] == ONE
            omitted = [row[:] for row in b]
            omitted[0][modes] = ZERO
            column = [row[modes] for row in omitted]
            assert sum((column[i]*column[modes+i] for i in range(modes)), ZERO) == ONE
        records.append({'case': name, 'modes': modes, 'twice_mean': int(2*mean),
                        'bogoliubov_orders': signed, 'fock_orders': actual})
    return records


def complete_weights(interval, side, heights):
    m, q = interval//2, 1 << (interval//2)
    first = (0, 1) if side == 0 else (0, 2)
    second = tuple(leg for leg in range(4) if leg not in first)
    def weights(legs):
        return [sum(heights[leg*m+j] for leg, state in zip(legs, (a, b))
                    for j in range(m) if state >> j & 1)
                for a, b in product(range(q), repeat=2)]
    return weights(first), weights(second)


def exact_margin(interval, side, heights):
    raw = laurent_kernel(interval)
    _, _, rows, columns = bogoliubov_data(interval, side)
    values = []
    for group in (rows, columns):
        for position, i in enumerate(group):
            for j in group[position+1:]:
                if raw[i][j]:
                    order = min(raw[i][j]) if side == 0 else -max(raw[i][j])
                    values.append(order+heights[i]+heights[j])
    return min(values)


def reflection_control(interval, side, heights, exact_fock):
    b, rows, columns = closed_matrix(interval, side, infinity=bool(side))
    rw = [heights[i] for i in rows]
    cw = [heights[i] for i in columns]
    rw += [-value for value in rw]
    cw += [-value for value in cw]
    gauged = [[entry*Z**(rw[i]+cw[j]) for j, entry in enumerate(row)] for i, row in enumerate(b)]
    signed = rational_smith(gauged)
    mean = Fraction(interval*interval if side == 0 else 0)+Fraction(sum(heights), 2)
    predicted = subset_prediction(signed, mean)
    native = fock_spectrum(interval, side, heights)
    assert native['bogoliubov_profile']['orders'] == signed
    assert native['mean'] == mean
    assert native['minimum'] == predicted[0]
    values = tensor(interval)['values']
    matrix = flatten(values, interval//2, bool(side))
    rows, columns = complete_weights(interval, side, heights)
    if exact_fock:
        literal = [[from_polynomial(value[0], infinity=bool(side))*Z**(rows[i]+columns[j])
                    for j, value in enumerate(row)] for i, row in enumerate(matrix)]
        actual = rational_smith(literal)
    else:
        bits = max(max(value).bit_length() for value in values)
        oracle = MatrixOracle(matrix, 1 << interval, precision=max(512, 1 << (bits-1).bit_length()))
        actual = oracle.profile(rows, columns, infinity=bool(side))['orders']
    assert actual == predicted
    assert sum(actual) == mean*(1 << interval)
    return {'interval': interval, 'side': side, 'heights': heights,
            'internal_regularity_margin': exact_margin(interval, side, heights),
            'twice_mean': int(2*mean), 'bogoliubov_orders': signed,
            'fock_orders': predicted, 'exact_rational_fock': exact_fock}


def verify():
    rng = random.Random(77231)
    structural = []
    for interval in (4, 8, 16):
        for side in (0, 1):
            matrix, _, _ = closed_matrix(interval, side)
            assert matrix == kernel_matrix(interval, side)
            quadratic_control(matrix)
            structural.append({'interval': interval, 'side': side,
                               'quadratic_columns': 2*interval,
                               'polar_pairings': 4*interval*interval})
    complete = []
    for interval in (4, 8):
        for side in (0, 1):
            controls = [[0]*(2*interval), list(logical_heights(interval, 'zero' if side == 0 else 'infinity'))]
            odd = [0]*(2*interval)
            odd[0] = 1
            controls.append(odd)
            for trial in range(5 if interval == 4 else 2):
                controls.append([rng.randrange(-4*interval, 4*interval+1) for _ in range(2*interval)])
            for heights in controls:
                complete.append(reflection_control(interval, side, heights, exact_fock=interval == 4))
    assert any(record['internal_regularity_margin'] < 0 for record in complete)
    assert any(record['twice_mean'] % 2 for record in complete)
    root = Path(__file__).resolve().parent.parent
    files = ['research/clifford-gaussian-smith-profiles.md',
             'tools/reflection_unrestricted_gaussian.py',
             'tools/check_reflection_unrestricted_gaussian.py',
             'tools/reflection_dual_rail_gaussian.py']
    return {'status': 'pass', 'toy_controls': toy_controls(),
            'closed_matrix_controls': structural, 'complete_profile_controls': complete,
            'nonregular_complete_profiles': sum(record['internal_regularity_margin'] < 0 for record in complete),
            'exact_rational_complete_profiles': sum(record['exact_rational_fock'] for record in complete),
            'complete_fock_orders_checked': sum(len(record['fock_orders']) for record in complete),
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in files},
            'scope': 'Exact Clifford and complete Fock controls for the uniform unrestricted Gaussian profile theorem; no improved exponent claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'],
                      'toy_controls': len(result['toy_controls']),
                      'complete_profiles': len(result['complete_profile_controls']),
                      'nonregular_profiles': result['nonregular_complete_profiles'],
                      'exact_rational_complete_profiles': result['exact_rational_complete_profiles'],
                      'complete_fock_orders_checked': result['complete_fock_orders_checked']}))


if __name__ == '__main__':
    main()
