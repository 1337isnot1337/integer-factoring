#!/usr/bin/env python3
"""Eight-k-mode literal Gaussian transport for homogeneous physical codes.

This transports the internal Gaussian factors as well as the compound
matrix. Deletion closure is not assumed. The explicit finite constructor
below computes the k=1 occupied/vacant controls; larger k operator images
are provided without claiming that their exponentially large matrices
have been evaluated.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from check_tensor_hot_compounds_f16 import binary_product as pmul


def dot(a, b):
    value = 0
    for x, y in zip(a, b):
        value ^= pmul(x, y)
    return value


def pairing(a, b):
    return dot(a[0], b[1]) ^ dot(a[1], b[0])


def images(M, k):
    if type(M) is not int or type(k) is not int or k < 1 or k & (k-1) or M < 4*k or M & (M-1):
        raise ValueError('Use powers of two k and M, with M>=4k.')
    N, n = M//k, 8*k
    S = sum(1 << (2*k*j) for j in range(N-2))
    tail_power = 1 << (2*M-2*k)
    reflected = (1 << (2*M)) ^ 1
    index = lambda group, position, residue: group*4*k+position*k+residue

    def unit(creation, place):
        out = [[0]*n, [0]*n]
        out[0 if creation else 1][place] = 1
        return out

    Rx, Ry, Cv, Cu = [], [], [], []
    for a in range(k):
        rcx, rax = unit(True, index(0, 0, a)), unit(False, index(0, 0, a))
        rcy, ray = unit(True, index(1, 0, a)), unit(False, index(1, 0, a))
        ccv, cav = unit(True, index(0, 0, a)), unit(False, index(0, 0, a))
        ccu, cau = unit(True, index(1, 0, a)), unit(False, index(1, 0, a))
        for c in range(k):
            fca, fac, gac = 3 << (2*(c & ~a)), 3 << (2*(a & ~c)), 3 << (2*(a & c))
            for position, weight in ((0, 1), (1, S << (2*k)), (3, tail_power)):
                rcx[0][index(1, position, c)] ^= pmul(fca, weight)
            for position, weight in ((0, 1), (2, S), (3, 1)):
                rax[0][index(1, position, c)] ^= pmul(gac, weight)
                ray[0][index(0, position, c)] ^= pmul(gac, weight)
            for position in (0, 1, 3):
                ray[1][index(0, position, c)] ^= fac
            for position in (0, 2, 3):
                ccv[1][index(1, position, c)] ^= gac
                ccu[1][index(0, position, c)] ^= gac
            for position, weight in ((0, 1), (1, 1 << (2*k)), (3, tail_power)):
                cav[1][index(1, position, c)] ^= pmul(fca, weight)
            for position, weight in ((0, 1), (1, S), (3, 1)):
                ccu[0][index(0, position, c)] ^= pmul(fac, weight)
        ray[0][index(1, 3, k-1-a)] ^= reflected
        ccu[1][index(1, 3, k-1-a)] ^= reflected
        Rx.append((rcx, rax))
        Ry.append((rcy, ray))
        Cv.append((ccv, cav))
        Cu.append((ccu, cau))
    for groups in ((Rx, Ry), (Cv, Cu)):
        operators = groups[0]+groups[1]
        for i, (ci, ai) in enumerate(operators):
            for j, (cj, aj) in enumerate(operators):
                assert pairing(ci, cj) == pairing(ai, aj) == 0
                assert pairing(ci, aj) == int(i == j)
    result = (Rx, Ry, Cv, Cu)
    assert max(value.bit_length()-1 for group in result for pair in group for operator in pair
               for vector in operator for value in vector) <= 2*M
    physical_active = [c+k*B for B in (0, 1, N-2, N-1) for c in range(k)]
    active = physical_active + [M+j for j in physical_active]
    orders = [2*j+1 for j in physical_active] + [2*M+2*j+1 for j in physical_active]
    assert sum(orders) == 16*k*M
    return result, orders, active


def left_action(matrix, operator):
    n = len(operator[0])
    size = 1 << n
    out = [[0]*size for _ in range(size)]
    for source, row in enumerate(matrix):
        nonzero = [(j, value) for j, value in enumerate(row) if value]
        if not nonzero:
            continue
        for creation, vector in ((True, operator[0]), (False, operator[1])):
            for mode, coefficient in enumerate(vector):
                if not coefficient or bool(source >> mode & 1) == creation:
                    continue
                target = source ^ (1 << mode)
                for column, value in nonzero:
                    out[target][column] ^= pmul(coefficient, value)
    return out


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def right_action(matrix, operator):
    return transpose(left_action(transpose(matrix), (operator[1], operator[0])))


def number_projector(matrix, operator_pair, occupied, right=False):
    creator, annihilator = operator_pair
    if right:
        product = right_action(right_action(matrix, creator), annihilator)
    else:
        product = left_action(left_action(matrix, annihilator), creator)
    if occupied:
        return product
    return [[a ^ b for a, b in zip(row, other)] for row, other in zip(matrix, product)]


def k1_matrix(M, occupied):
    if occupied not in (0, 1):
        raise ValueError('Choose the vacant or occupied lowest mode.')
    (rx, ry, cv, cu), orders, active = images(M, 1)
    size = 1 << len(orders)
    matrix = [[0]*size for _ in range(size)]
    for state in range(size):
        matrix[state][state] = 1 << sum(order for i, order in enumerate(orders) if state >> i & 1)
    for operators, right in ((rx, False), (ry, False), (cv, True), (cu, True)):
        matrix = number_projector(matrix, operators[0], occupied, right)
    degree = max(value.bit_length()-1 for row in matrix for value in row)
    assert degree <= 32*M
    return matrix, degree, active


def binary_rank(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length()-1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def k1_profile(M, occupied):
    from native_matrix_smith_f4 import MatrixOracle
    matrix, degree, active = k1_matrix(M, occupied)
    endpoint = [sum((value.bit_count() % 2) << j for j, value in enumerate(row)) for row in matrix]
    assert binary_rank(endpoint) == 64
    converted = [[(value, 0) for value in row] for row in matrix]
    precision = 128
    while precision <= degree:
        precision *= 2
    attempts = []
    while precision <= 8192:
        try:
            result = MatrixOracle(converted, 64, precision=precision).profile()['orders']
            attempts.append([precision, True])
            break
        except ArithmeticError:
            attempts.append([precision, False])
            precision *= 2
    else:
        raise ArithmeticError('The complete k=1 core exceeded the precision budget.')
    assert len(result) == 64 and all(x >= 0 and x == round(x) for x in result)
    return {'M': M, 'k': 1, 'occupied': occupied, 'active_modes': 8,
            'active_ambient_order': 256, 'active_generic_rank': 64,
            'orders': sorted(map(round, result)), 'actual_coefficient_degree': degree,
            'coefficient_degree_bound': 32*M, 'precision_attempts': attempts,
            'active_full_mode_indices': active,
            'spectator_orders': [2*j+1 for j in range(2, M-2)]
                               + [2*M+2*j+1 for j in range(2, M-2)],
            'physical_hot_heights': [occupied],
            'physical_tail_heights': [2*j+1 for j in range(1, M)],
            'physical_alphabet_dimension': 1 << (M-1)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, nargs='+', default=[4, 8, 16, 32])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    names = [Path(__file__).name, 'check_tensor_hot_compounds_f16.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_weighted_smith_f4.py', 'native_rectangular_smith_f4.cpp']
    report = {'scope': __doc__, 'completed': False, 'records': [],
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for M in args.modes:
        for occupied in (0, 1):
            started = time.monotonic()
            record = k1_profile(M, occupied)
            record['seconds'] = time.monotonic()-started
            report['records'].append(record)
            args.output.write_text(json.dumps(report, indent=2)+'\n')
            print(json.dumps({k: record[k] for k in ('M', 'occupied', 'active_generic_rank',
                                                   'actual_coefficient_degree', 'seconds')}), flush=True)
    report['completed'] = True
    args.output.write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
