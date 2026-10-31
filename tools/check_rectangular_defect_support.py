#!/usr/bin/env python3
"""Exact DVR-compatible concentration of interval crossing defect operators.

This verifies rational one-particle identities and finite CAR supports.
It does not compute a pressure bound or claim an improved factoring score.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


def pmul(a, b):
    answer = 0
    while b:
        bit = b & -b
        answer ^= a << (bit.bit_length()-1)
        b ^= bit
    return answer


def pdivmod(a, b):
    if not b:
        raise ZeroDivisionError
    q = 0
    while a and a.bit_length() >= b.bit_length():
        shift = a.bit_length()-b.bit_length()
        q ^= 1 << shift
        a ^= b << shift
    return q, a


def pgcd(a, b):
    while b:
        a, b = b, pdivmod(a, b)[1]
    return a


class Rat:
    __slots__ = ('n', 'd')
    def __init__(self, n=0, d=1):
        if not d:
            raise ZeroDivisionError
        if not n:
            self.n, self.d = 0, 1
        else:
            common = pgcd(n, d)
            self.n, rn = pdivmod(n, common)
            self.d, rd = pdivmod(d, common)
            assert not rn and not rd
    def __bool__(self):
        return bool(self.n)
    def __add__(self, other):
        if not isinstance(other, Rat):
            other = Rat(other)
        return Rat(pmul(self.n, other.d)^pmul(other.n, self.d), pmul(self.d, other.d))
    __radd__ = __add__
    def __mul__(self, other):
        if not isinstance(other, Rat):
            other = Rat(other)
        return Rat(pmul(self.n, other.n), pmul(self.d, other.d))
    __rmul__ = __mul__
    def __truediv__(self, other):
        if not isinstance(other, Rat):
            other = Rat(other)
        return Rat(pmul(self.n, other.d), pmul(self.d, other.n))
    def __pow__(self, exponent):
        if exponent < 0:
            return Rat(self.d, self.n)**(-exponent)
        answer, base = Rat(1), self
        while exponent:
            if exponent & 1:
                answer = answer*base
            base = base*base
            exponent >>= 1
        return answer
    def __eq__(self, other):
        if not isinstance(other, Rat):
            other = Rat(other)
        return self.n == other.n and self.d == other.d
    def __repr__(self):
        return f'Rat({self.n},{self.d})'
    def order(self):
        if not self.n:
            return None
        return (self.n & -self.n).bit_length()-(self.d & -self.d).bit_length()
    def pair(self):
        return [self.n, self.d]


ZERO, ONE, Z, UNIT = Rat(), Rat(1), Rat(2), Rat(3)
T, Q = Z**2, UNIT**2


def vec(n, entries=()):
    answer = [ZERO]*n
    for index, value in entries:
        answer[index] = value
    return answer


def dot(a, b):
    answer = ZERO
    for x, y in zip(a, b):
        if x and y:
            answer = answer+x*y
    return answer


def mm(a, b):
    bt = list(zip(*b))
    return [[dot(row, column) for column in bt] for row in a]


def identity(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def infinity_matrices(M):
    """Direct Gaussian matrix and the explicit regular Smith operations."""
    n = 2*M
    matrix = [[ZERO]*n for _ in range(n)]
    U = [[ZERO]*n for _ in range(n)]
    Ui = [[ZERO]*n for _ in range(n)]
    V = [[ZERO]*n for _ in range(n)]
    Vi = [[ZERO]*n for _ in range(n)]
    sigma = [Z**(1-2*M+2*a) for a in range(n)]
    for a in range(M):
        U[a][M+a] = ONE
        U[M+a][a] = ONE/UNIT
        Ui[a][M+a] = UNIT
        Ui[M+a][a] = ONE
        V[a][a] = Vi[a][a] = ONE
        for b in range(M):
            row_v, row_y = 2*a-M+1, -M
            col_u, col_x = 2*b-M+1, M
            if a & b == a:
                matrix[a][b] = (UNIT/Z)**(2*a+1) * Z**(row_v+col_u)
            if a & b == b:
                matrix[a][M+b] = Z**(-2*M+row_v+col_x)
            if a == b:
                matrix[M+a][b] = Z**(row_y+col_u)
            matrix[M+a][M+b] = UNIT*Z**(-1-2*(M-1-(a | b))+row_y+col_x)
            if a & b == a:
                U[M+a][M+b] = Q**a
                Ui[a][b] = UNIT*Q**a
                V[a][M+b] = T**(b-a)*Q**(M-b)
            if a & b == b:
                V[M+a][M+b] = UNIT
                Vi[M+a][M+b] = ONE/UNIT
            Vi[a][M+b] = UNIT*T**((a | b)-a)
    return matrix, U, Ui, V, Vi, sigma


def infinity_images(M):
    n, top = 2*M, M-1
    images = {}
    for i in (0, top):
        images[f'out_v{i}_create'] = (vec(n, [(M+i, ONE/UNIT)]), vec(n))
        images[f'out_y{i}_create'] = (vec(n, [(i, ONE)]+[(M+a, Q**a) for a in range(M) if a & i == a]), vec(n))
        images[f'out_v{i}_annihilate'] = (vec(n), vec(n, [(a, UNIT*Q**i) for a in range(M) if i & a == i]+[(M+i, UNIT)]))
        images[f'out_y{i}_annihilate'] = (vec(n), vec(n, [(i, ONE)]))
        images[f'in_u{i}_create'] = (vec(n, [(i, Z**(1-2*M+2*i))]), vec(n))
        images[f'in_u{i}_annihilate'] = (vec(n), vec(n, [(i, Z**(2*M-1-2*i))]+[(M+a, Z**(-1-2*i)*Q**(M-a)) for a in range(M) if i & a == i]))
    images['in_x0_annihilate'] = (vec(n), vec(n, [(M, UNIT/Z)]))
    return images


def zero_images(M):
    n, top = 2*M, M-1
    images = {}
    sigma = [Z**(2*a+1) for a in range(M)]
    for i in (0, top):
        images[f'left_y{i}_create'] = (vec(n, [(M+i, ONE)]), vec(n))
        images[f'left_y{i}_annihilate'] = (vec(n, [(a, UNIT*T**(a & i)) for a in range(M)]), vec(n, [(M+i, ONE)]))
        images[f'right_u{i}_create'] = (vec(n, [(M+i, sigma[i])]), vec(n, [(a, UNIT*T**(i & a)/sigma[a]) for a in range(M)]))
        images[f'right_u{i}_annihilate'] = (vec(n), vec(n, [(M+i, ONE/sigma[i])]))
        images[f'right_v{i}_create'] = (vec(n, [(a, sigma[a]) for a in range(M) if a & i == a]), vec(n, [(M+a, UNIT*Q**i/sigma[a]) for a in range(M) if i & a == i]))
        images[f'right_v{i}_annihilate'] = (vec(n), vec(n, [(a, ONE/sigma[a]) for a in range(M) if i & a == i]))
    images['left_x0_create'] = (vec(n, [(0, ONE)]), vec(n))
    return images


def pairing(first, second):
    return dot(first[0], second[1])+dot(first[1], second[0])


def car_checks(images, prefixes):
    count = 0
    for first in prefixes:
        for second in prefixes:
            for kind1 in ('create', 'annihilate'):
                for kind2 in ('create', 'annihilate'):
                    expected = int(first == second and kind1 != kind2)
                    assert pairing(images[first+'_'+kind1], images[second+'_'+kind2]) == expected
                    count += 1
    return count


def apply_shear(images, row, column, coefficient, orders, trace):
    assert row != column and coefficient
    order = coefficient.order()
    companion_order = order+orders[column]-orders[row]
    assert order >= 0 and companion_order >= 0
    for creation, annihilation in images.values():
        creation[row] = creation[row]+coefficient*creation[column]
        annihilation[column] = annihilation[column]+coefficient*annihilation[row]
    trace.append({'row': row, 'column': column, 'coefficient': coefficient.pair(),
                  'coefficient_order': order, 'companion_order': companion_order})


def concentrate_infinity(M):
    images = infinity_images(M)
    orders = [1-2*M+2*a for a in range(2*M)]
    trace = []
    if M < 4:
        return images, list(range(2*M)), orders, trace
    active = sorted({0, 1, M-1, M, M+1, 2*M-2, 2*M-1})
    pivot = M-2
    for j in range(2, M-2):
        apply_shear(images, M+j, M+pivot, Q**(j-pivot), orders, trace)
    for j in range(2, M-1):
        apply_shear(images, 1, j, ONE, orders, trace)
    for j in range(2, M-2):
        apply_shear(images, M+1, M+j, Q**(1-j), orders, trace)
    assert all(not a[i] and not b[i] for a, b in images.values() for i in range(2*M) if i not in active)
    return images, active, orders, trace


def concentrate_zero(M):
    images = zero_images(M)
    orders = [2*a+1 for _ in range(2) for a in range(M)]
    trace = []
    if M < 8:
        return images, list(range(2*M)), orders, trace
    active_first = {0, 1, 2, M-3, M-2, M-1}
    spectators = [j for j in range(M) if j not in active_first]
    high = M-2
    for j in spectators+[high]:
        apply_shear(images, j, 1, T**(j-1), orders, trace)
    for j in spectators:
        apply_shear(images, j, high, (ONE+T**(j-1))/(ONE+T**(high-1)), orders, trace)
    high_b, low_b = M-3, 2
    for j in spectators+[low_b]:
        apply_shear(images, high_b, j, T**(high_b-j), orders, trace)
    for j in spectators:
        apply_shear(images, low_b, j, (ONE+T**(high_b-j))/(ONE+T**(high_b-low_b)), orders, trace)
    for j in range(1, M-2):
        apply_shear(images, M+(M-2), M+j, T**((M-2)-j), orders, trace)
    active = sorted(active_first | {M, 2*M-2, 2*M-1})
    assert all(not a[i] and not b[i] for a, b in images.values() for i in range(2*M) if i not in active)
    return images, active, orders, trace


def check_infinity_diagonalization(M):
    matrix, U, Ui, V, Vi, sigma = infinity_matrices(M)
    assert mm(U, Ui) == mm(Ui, U) == identity(2*M)
    assert mm(V, Vi) == mm(Vi, V) == identity(2*M)
    expected = [[sigma[i] if i == j else ZERO for j in range(2*M)] for i in range(2*M)]
    assert mm(mm(U, matrix), V) == expected
    for transformation in (U, Ui, V, Vi):
        assert all(not x or x.order() >= 0 for row in transformation for x in row)
    images = infinity_images(M)
    for i in (0, M-1):
        for prefix, index in ((f'out_v{i}', i), (f'out_y{i}', M+i)):
            assert images[prefix+'_create'][0] == [U[j][index] for j in range(2*M)]
            assert images[prefix+'_annihilate'][1] == Ui[index]
        assert images[f'in_u{i}_create'][0] == [sigma[j]*Vi[j][i] for j in range(2*M)]
        assert images[f'in_u{i}_annihilate'][1] == [V[i][j]/sigma[j] for j in range(2*M)]
    assert images['in_x0_annihilate'][1] == [V[M][j]/sigma[j] for j in range(2*M)]
    return {'M': M, 'matrix_size': 2*M, 'two_sided_inverse_and_diagonalization': True,
            'regular_transformations_and_inverses': True, 'transported_core_generators_match': True}


def check_zero_transport(M):
    images = zero_images(M)
    for i in (0, M-1):
        for j in range(M):
            original = ZERO
            for a in range(M):
                if a & i == a and a & j == a:
                    original = original+UNIT*Q**a
            assert original == UNIT*T**(i & j)
    return {'M': M, 'pascal_transport_identity_entries': 2*M}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--artifact-directory', type=Path)
    args = parser.parse_args()
    diagonalizations = [check_infinity_diagonalization(M) for M in (2, 4, 8, 16)]
    records = []
    for M in (2, 4, 8, 16, 32, 64, 128):
        check_zero_transport(M)
        for place, initial, procedure, groups in (
            ('infinity', infinity_images, concentrate_infinity,
             [[f'out_{letter}{i}' for letter in ('v', 'y') for i in (0, M-1)],
              [f'in_u{i}' for i in (0, M-1)]]),
            ('zero', zero_images, concentrate_zero,
             [[f'left_y{i}' for i in (0, M-1)],
              [f'right_{letter}{i}' for letter in ('u', 'v') for i in (0, M-1)]])):
            before = initial(M)
            before_checks = sum(car_checks(before, group) for group in groups)
            images, active, orders, trace = procedure(M)
            after_checks = sum(car_checks(images, group) for group in groups)
            # All cross-pairings are preserved too; they need not be canonical.
            for name, first in before.items():
                for other_name, second in before.items():
                    assert pairing(first, second) == pairing(images[name], images[other_name])
            active_images = {name: {'creation': [a[i].pair() for i in active],
                                    'annihilation': [b[i].pair() for i in active]}
                             for name, (a, b) in images.items()}
            data = {'place': place, 'M': M, 'active_indices': active,
                    'active_orders': [orders[i] for i in active],
                    'spectator_orders': [orders[i] for i in range(2*M) if i not in active],
                    'active_images': active_images, 'regular_shears': trace}
            record = {'place': place, 'M': M, 'active_modes': len(active),
                      'spectator_modes': 2*M-len(active), 'regular_shears': len(trace),
                      'canonical_pairings_checked': before_checks+after_checks,
                      'all_cross_pairings_checked': len(images)**2,
                      'coefficient_sha256': sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()}
            records.append(record)
            if args.artifact_directory:
                args.artifact_directory.mkdir(parents=True, exist_ok=True)
                (args.artifact_directory/f'defect-{place}-M{M}.json').write_text(json.dumps(data, separators=(',', ':'))+'\n')
    result = {'status': 'exact rational Gaussian transforms and bounded defect supports; no Smith pressure claim',
              'infinity_diagonalizations': diagonalizations, 'cases': records}
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
