#!/usr/bin/env python3
"""Full-mode conjugation checks for the literal eight-k active-core lemma."""
import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

import tensor_hot_compounds_f16 as field
import tensor_literal_hot_pfaffians_f16 as literal
from tensor_literal_full_tail_core import dot, images, k1_profile, pairing, pmul


def multiply(a, b):
    columns = list(zip(*b))
    return [[dot(row, column) for column in columns] for row in a]


def shift(value, power):
    if power >= 0:
        return value << power
    assert not value & ((1 << -power)-1)
    return value >> -power


def concentration(M, k):
    """Independent rational-basis numerators, with common denominator S*T."""
    N, length = M//k, M//k-2
    S = sum(1 << (2*k*j) for j in range(length))
    T = 1 ^ (1 << (2*k*(length-1)))
    denominator = pmul(S, T)
    assert denominator & 1
    C, U = [[0]*M for _ in range(M)], [[0]*M for _ in range(M)]
    for c in range(k):
        for B in (0, N-1):
            C[c+k*B][c+k*B] = U[c+k*B][c+k*B] = denominator
        first, last = c+k, c+k*(N-2)
        for j in range(length):
            place = c+k*(j+1)
            C[place][first] = T << (2*k*j)
            C[place][last] = T
            U[first][place] = denominator
            U[last][place] = denominator << (2*k*(length-1-j))
        for j in range(1, length-1):
            place = c+k*(j+1)
            C[place][place] = denominator
            C[first][place] = pmul(S, 1 ^ (1 << (2*k*(length-1-j))))
            C[last][place] = pmul(S, (1 << (2*k*(length-1-j))) ^ (1 << (2*k*(length-1))))
            for h in range(length):
                U[place][c+k*(h+1)] = (T << (2*k*j)) ^ (T << (2*k*(length-1-h)))
            U[place][place] ^= denominator
    return C, U, denominator


def block_double(matrix):
    n = len(matrix)
    return [row+[0]*n for row in matrix] + [[0]*n+row for row in matrix]


def full_operators(M, k):
    """Compose both internal Gaussian shears before any tail concentration."""
    F = [[3 << (2*(a & ~b)) for b in range(M)] for a in range(M)]
    A = [[3 << (2*(a & b)) for b in range(M)] for a in range(M)]
    FA = multiply(F, A)
    AFt = multiply(A, [list(row) for row in zip(*F)])
    reflected = (1 << (2*M)) ^ 1
    for a, b in product(range(M), repeat=2):
        assert FA[a][b] == AFt[a][b] == (reflected if a+b == M-1 else 0)

    def unit(creation, place):
        op = [[0]*(2*M), [0]*(2*M)]
        op[0 if creation else 1][place] = 1
        return op

    result = [[], [], [], []]
    for a in range(k):
        rcx, rax = unit(True, a), unit(False, a)
        rcy, ray = unit(True, M+a), unit(False, M+a)
        ccv, cav = unit(True, a), unit(False, a)
        ccu, cau = unit(True, M+a), unit(False, M+a)
        for b in range(M):
            rcx[0][M+b] ^= F[b][a]
            rax[0][M+b] ^= A[a][b]
            ray[1][b] ^= F[a][b]
            ray[0][b] ^= A[a][b]
            ray[0][M+b] ^= FA[b][a]
            ccv[1][M+b] ^= A[a][b]
            cav[1][M+b] ^= F[b][a]
            ccu[0][b] ^= F[a][b]
            ccu[1][b] ^= A[a][b]
            ccu[1][M+b] ^= AFt[a][b]
        for side, value in enumerate(((rcx, rax), (rcy, ray), (ccv, cav), (ccu, cau))):
            result[side].append(value)
    return result, 2*M*M


def operator_controls():
    counts = {'matrix_inverse_entries': 0, 'weighted_regular_entries': 0,
              'reflection_product_entries': 0, 'transported_operator_entries': 0,
              'canonical_pairings': 0, 'endpoint_operator_entries': 0}
    records = []
    for k, M in ((1, 4), (1, 8), (1, 16), (2, 8), (2, 16), (2, 32), (4, 16), (4, 32)):
        C, U, denominator = concentration(M, k)
        squared = pmul(denominator, denominator)
        for matrix in (multiply(C, U), multiply(U, C)):
            for i, j in product(range(M), repeat=2):
                assert matrix[i][j] == (squared if i == j else 0)
                counts['matrix_inverse_entries'] += 1
        V = [[shift(C[i][j], 2*(j-i)) for j in range(M)] for i in range(M)]
        Vinv = [[shift(U[i][j], 2*(j-i)) for j in range(M)] for i in range(M)]
        counts['weighted_regular_entries'] += 2*M*M
        fullC, fullU, fullV, fullVinv = map(block_double, (C, U, V, Vinv))
        expected, orders, active = images(M, k)
        operators, checked = full_operators(M, k)
        counts['reflection_product_entries'] += checked
        for side, group in enumerate(operators):
            creator_map, annihilator_map = (fullU, fullC) if side < 2 else (fullVinv, fullV)
            for a, pair in enumerate(group):
                for which, operator in enumerate(pair):
                    creator = [dot(row, operator[0]) for row in creator_map]
                    annihilator = [dot(operator[1], column) for column in zip(*annihilator_map)]
                    for component, actual in enumerate((creator, annihilator)):
                        target = [0]*(2*M)
                        for index, value in zip(active, expected[side][a][which][component]):
                            target[index] = pmul(denominator, value)
                        assert actual == target
                        counts['transported_operator_entries'] += 2*M
        for groups in (expected[:2], expected[2:]):
            modes = groups[0]+groups[1]
            for i, (ci, ai) in enumerate(modes):
                for j, (cj, aj) in enumerate(modes):
                    assert pairing(ci, cj) == pairing(ai, aj) == 0
                    assert pairing(ci, aj) == int(i == j)
                    counts['canonical_pairings'] += 3
        for side, group in enumerate(expected):
            group_index = 0 if side in (0, 2) else 1
            for a, pair in enumerate(group):
                for creation, operator in enumerate(pair):
                    for component, vector in enumerate(operator):
                        endpoint = [value.bit_count() % 2 for value in vector]
                        target = [0]*(8*k)
                        if component == creation:
                            target[group_index*4*k+a] = 1
                        assert endpoint == target
                        counts['endpoint_operator_entries'] += len(vector)
        records.append({'k': k, 'M': M, 'active_indices': active,
                        'active_diagonal_orders': orders, 'concentration_denominator': denominator})
    return counts, records


def literal_controls():
    records = []
    for occupied in (0, 1):
        maps = [{w: 1} for w in range(16) if w & 1 == occupied]
        q = len(maps)
        orders = []
        for parity in (0, 1):
            pairs = [(i, j) for i in range(q) for j in range(q)
                     if (next(iter(maps[i])).bit_count()+next(iter(maps[j])).bit_count()) % 2 == parity]
            rows = [field.tensor_vector(maps[i], maps[j], 4) for i, j in pairs]
            columns = [field.tensor_vector(maps[j], maps[i], 4) for i, j in pairs]
            orders.extend(literal.compound_block(4, rows, columns)['orders_t'])
        record = k1_profile(4, occupied)
        assert sorted(orders) == record['orders']
        record['literal_orders_compared'] = len(orders)
        record['literal_code_is_deletion_closed'] = not occupied
        records.append(record)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/tensor-literal-full-tail-controls.json'))
    args = parser.parse_args()
    started = time.monotonic()
    counts, records = operator_controls()
    literal_records = literal_controls()
    names = [Path(__file__).name, 'tensor_literal_full_tail_core.py',
             'tensor_literal_hot_pfaffians_f16.py', 'tensor_literal_hot_pfaffians_f16.cpp',
             'tensor_hot_compounds_f16_fast.cpp', 'tensor_hot_compounds_f16.py',
             'tensor_hot_compounds_f16.cpp', 'tensor_orthogonal_scaling.py',
             'check_tensor_hot_compounds_f16.py', 'coherent_hot_compounds_f4.cpp',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp', 'native_weighted_smith_f4.py']
    report = {'passed': True, 'counts': counts, 'operator_controls': records,
              'literal_profile_controls': literal_records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}, 'seconds': time.monotonic()-started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'passed': True, 'counts': counts,
                      'literal_profile_orders': sum(r['literal_orders_compared'] for r in literal_records),
                      'seconds': report['seconds']}), flush=True)


if __name__ == '__main__':
    main()
