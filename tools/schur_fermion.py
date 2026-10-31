#!/usr/bin/env python3
"""Direct Schur modes and crossing coefficients from charged-fermion determinants.

Research implementation.  Source modes are Appendix B at anthropics/formal-math
16a5ce41304fcdaea4e0be7c132e1c361dc26da5; Schur normalization is from PR #30.
The independent comparison uses rational power sums only in the audit.
"""

from __future__ import annotations

from functools import lru_cache
from math import comb
from fractions import Fraction


@lru_cache(None)
def frobenius(partition):
    """Return decreasing arm and leg coordinates, accepting increasing parts."""
    lam = tuple(sorted(partition, reverse=True))
    rank = sum(row >= i + 1 for i, row in enumerate(lam))
    arms = tuple(lam[i] - i - 1 for i in range(rank))
    legs = tuple(sum(row > i for row in lam) - i - 1 for i in range(rank))
    return arms, legs


@lru_cache(None)
def state_fields(partition):
    """Paired creator list and sign in the Schur--semi-infinite-wedge map."""
    arms, legs = frobenius(partition)
    fields = tuple(field for a, b in zip(arms, legs)
                   for field in ((1, a), (-1, b)))
    return fields, (-1) ** sum(legs)


@lru_cache(None)
def conjugate(partition):
    return tuple(sorted(sum(x >= k for x in partition)
                        for k in range(1, max(partition, default=0) + 1)))


def wick_determinant(fields, contraction, sign=1, prime=None):
    plus = [i for i, field in enumerate(fields) if field[1] == 1]
    minus = [i for i, field in enumerate(fields) if field[1] == -1]
    assert len(plus) == len(minus)
    matrix = [[contraction(i, j) for j in minus] for i in plus]
    m = len(plus)
    inversions = sum(i > j for i in plus for j in minus)
    sign *= (-1) ** (inversions + m * (m - 1) // 2)
    value = sign * determinant(matrix, prime)
    return value if prime is None else value % prime


def crossing_value(x, y, u, v, z, prime=None):
    """Value of z^(2h) sum_s A_s(x,y,u,v) z^s for z != 0,1."""
    h = sum(x)
    assert h == sum(y) == sum(u) == sum(v)
    assert z not in (0, 1)
    q = z - 1
    xf, sx = state_fields(conjugate(x))
    vf, sv = state_fields(v)
    uf, su = state_fields(u)
    yf, sy = state_fields(y)
    fields = (tuple((0, -species, k) for species, k in reversed(uf))
              + tuple((1, species, k) for species, k in xf)
              + tuple((2, species, k) for species, k in vf)
              + tuple((3, species, k) for species, k in yf))
    power = (lambda a, b: Fraction(a) ** b) if prime is None else (
        lambda a, b: pow(a % prime, b, prime))

    def contraction(i, j):
        if i > j:
            return -contraction(j, i)
        bi, specie_i, ki = fields[i]
        bj, specie_j, kj = fields[j]
        if specie_i == specie_j or bi == bj:
            return 0
        if (bi, bj) == (0, 1):
            return sum(comb(kj, a) * power(q, a) * comb(ki, a)
                       * power(z, ki - a) for a in range(min(ki, kj) + 1))
        if (bi, bj) == (0, 2):
            return comb(ki, kj) if ki >= kj else 0
        if (bi, bj) == (0, 3):
            return int(ki == kj)
        if (bi, bj) == (1, 2):
            return (-1) ** ki * (comb(kj, ki) if kj >= ki else 0) * power(q, -kj - 1)
        if (bi, bj) == (1, 3):
            return sum(comb(ki, a) * power(q, a) * (-1) ** a
                       * comb(kj + a, a) * power(z, -kj - a - 1)
                       for a in range(ki + 1))
        assert (bi, bj) == (2, 3)
        return (-1) ** ki * comb(ki + kj, ki)

    value = wick_determinant(fields, contraction, sx * sv * su * sy, prime)
    value *= power(z, 2 * h) * power(q, h + len(xf) // 2)
    return value if prime is None else value % prime


def determinant(matrix, prime=None):
    n = len(matrix)
    if not n:
        return 1
    a = [list(row) for row in matrix]
    sign = 1
    if prime is not None:
        a = [[x % prime for x in row] for row in a]
        value = 1
        for j in range(n):
            pivot = next((i for i in range(j, n) if a[i][j]), None)
            if pivot is None:
                return 0
            if pivot != j:
                a[j], a[pivot] = a[pivot], a[j]
                value = -value
            pivot_value = a[j][j]
            value = value * pivot_value % prime
            inverse = pow(pivot_value, -1, prime)
            for i in range(j + 1, n):
                scale = a[i][j] * inverse % prime
                if scale:
                    for k in range(j + 1, n):
                        a[i][k] = (a[i][k] - scale * a[j][k]) % prime
        return value % prime
    previous = 1
    for j in range(n - 1):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            sign = -sign
        pivot_value = a[j][j]
        for i in range(j + 1, n):
            for k in range(j + 1, n):
                numerator = pivot_value * a[i][k] - a[i][j] * a[j][k]
                if isinstance(numerator, Fraction) or isinstance(previous, Fraction):
                    a[i][k] = numerator / previous
                else:
                    assert numerator % previous == 0
                    a[i][k] = numerator // previous
        previous = pivot_value
    return sign * a[-1][-1]


def schur_mode(parameter, output, input_, prime=None):
    """Coefficient of s_output in s_parameter[|output|-|input_|] s_input_."""
    param, sp = state_fields(parameter)
    out, so = state_fields(output)
    inp, si = state_fields(input_)
    # A field has (block, species, nonnegative mode or derivative index).
    fields = (tuple((0, -species, k) for species, k in reversed(out))
              + tuple((1, species, k) for species, k in param)
              + tuple((2, species, k) for species, k in inp))
    if prime == 2:
        plus = [i for i, field in enumerate(fields) if field[1] == 1]
        minus = [i for i, field in enumerate(fields) if field[1] == -1]
        pivots = {}
        for i in plus:
            row = 0
            for column, j in enumerate(minus):
                left, right = (fields[i], fields[j]) if i < j else (fields[j], fields[i])
                bi, _, ki = left
                bj, _, kj = right
                if bi == bj:
                    value = 0
                elif (bi, bj) == (0, 1):
                    value = int(kj & ~ki == 0)
                elif (bi, bj) == (1, 2):
                    value = int(ki & kj == 0)
                else:
                    assert (bi, bj) == (0, 2)
                    value = int(ki == kj)
                row |= value << column
            while row:
                lead = row.bit_length() - 1
                if lead not in pivots:
                    pivots[lead] = row
                    break
                row ^= pivots[lead]
            if not row:
                return 0
        return 1

    def contraction(i, j):
        if i > j:
            return -contraction(j, i)
        bi, specie_i, ki = fields[i]
        bj, specie_j, kj = fields[j]
        if specie_i == specie_j or bi == bj:
            return 0
        if (bi, bj) == (0, 1):
            return comb(ki, kj) if ki >= kj else 0
        if (bi, bj) == (1, 2):
            return (-1) ** ki * comb(ki + kj, ki)
        assert (bi, bj) == (0, 2)
        return int(ki == kj)

    return wick_determinant(fields, contraction, sp * so * si, prime)


def binary_polynomial_mod(value, divisor):
    while value.bit_length() >= divisor.bit_length():
        value ^= divisor << (value.bit_length() - divisor.bit_length())
    return value


def binary_polynomial_gcd(a, b):
    while b:
        a, b = b, binary_polynomial_mod(a, b)
    return a


class BinaryField:
    """Deterministically constructed F_(2^degree), with literal polynomial bits."""

    def __init__(self, degree):
        assert degree >= 1
        self.degree = degree
        self.order = 1 << degree
        for modulus in range(self.order | 1, 2 * self.order, 2):
            self.modulus = modulus
            power = 2
            good = True
            for i in range(1, degree + 1):
                power = self._multiply(power, power)
                if i <= degree // 2 and binary_polynomial_gcd(power ^ 2, modulus) != 1:
                    good = False
                    break
            if good and power == 2:
                break
        else:
            if degree == 1:
                self.modulus = 3
            else:
                raise AssertionError('irreducible polynomial search failed')
        # Small fields dominate the intended exact finite experiments.
        self.table = ([[self._multiply(a, b) for b in range(self.order)]
                       for a in range(self.order)] if degree <= 8 else None)
        self.inverse = [0] + [self.power(a, self.order - 2) for a in range(1, self.order)]
        assert all(self.multiply(a, self.inverse[a]) == 1 for a in range(1, self.order))

    def _multiply(self, a, b):
        value = 0
        while b:
            if b & 1:
                value ^= a
            b >>= 1
            a <<= 1
            if a & self.order:
                a ^= self.modulus
        return value

    def multiply(self, a, b):
        return self.table[a][b] if self.table is not None else self._multiply(a, b)

    def power(self, a, exponent):
        if exponent < 0:
            assert a
            exponent %= self.order - 1
        answer = 1
        while exponent:
            if exponent & 1:
                answer = self.multiply(answer, a)
            exponent >>= 1
            a = self.multiply(a, a)
        return answer

    def determinant(self, matrix):
        a = [list(row) for row in matrix]
        n = len(a)
        value = 1
        mul = self.multiply
        for j in range(n):
            pivot = next((i for i in range(j, n) if a[i][j]), None)
            if pivot is None:
                return 0
            a[j], a[pivot] = a[pivot], a[j]
            pivot_value = a[j][j]
            value = mul(value, pivot_value)
            inverse = self.inverse[pivot_value]
            for i in range(j + 1, n):
                scale = mul(a[i][j], inverse)
                if scale:
                    for k in range(j + 1, n):
                        a[i][k] ^= mul(scale, a[j][k])
        return value


def crossing_value_binary(x, y, u, v, z, field):
    """Characteristic-two folded kernel: every nonzero entry is a monomial."""
    h = sum(x)
    assert h == sum(y) == sum(u) == sum(v)
    assert z not in (0, 1)
    q = z ^ 1
    xf, _ = state_fields(conjugate(x))
    vf, _ = state_fields(v)
    uf, _ = state_fields(u)
    yf, _ = state_fields(y)
    fields = (tuple((0, -species, k) for species, k in reversed(uf))
              + tuple((1, species, k) for species, k in xf)
              + tuple((2, species, k) for species, k in vf)
              + tuple((3, species, k) for species, k in yf))
    plus = [i for i, f in enumerate(fields) if f[1] == 1]
    minus = [i for i, f in enumerate(fields) if f[1] == -1]
    powers_z = {a: field.power(z, a) for a in range(-2 * h - 1, 2 * h + 1)}
    powers_q = {a: field.power(q, a) for a in range(-h, 2 * h + 1)}

    def contraction(i, j):
        if i > j:
            i, j = j, i
        bi, si, ki = fields[i]
        bj, sj, kj = fields[j]
        if si == sj or bi == bj:
            return 0
        if (bi, bj) == (0, 1):
            return powers_z[ki & ~kj]
        if (bi, bj) == (0, 2):
            return int(kj & ~ki == 0)
        if (bi, bj) == (0, 3):
            return int(ki == kj)
        if (bi, bj) == (1, 2):
            return powers_q[-kj - 1] if ki & ~kj == 0 else 0
        if (bi, bj) == (1, 3):
            return powers_z[-kj - 1 - (ki & ~kj)]
        assert (bi, bj) == (2, 3)
        return int(ki & kj == 0)

    value = field.determinant([[contraction(i, j) for j in minus] for i in plus])
    value = field.multiply(value, powers_z[2 * h])
    return field.multiply(value, powers_q[h + len(xf) // 2])


@lru_cache(None)
def binary_interpolation_data(h):
    count = 4 * h + 1
    degree = (count + 2).bit_length()
    field = BinaryField(degree)
    points = tuple(range(2, count + 2))
    master = [1]
    for point in points:
        new = [0] * (len(master) + 1)
        for k, value in enumerate(master):
            new[k] ^= field.multiply(value, point)
            new[k + 1] ^= value
        master = new
    basis = []
    for point in points:
        # Synthetic division by Z+point; all roots in master are distinct.
        coefficients = [0] * count
        coefficients[-1] = master[-1]
        for k in range(count - 1, 0, -1):
            coefficients[k - 1] = master[k] ^ field.multiply(point, coefficients[k])
        assert master[0] == field.multiply(point, coefficients[0])
        denominator = 0
        for scalar in reversed(coefficients):
            denominator = field.multiply(denominator, point) ^ scalar
        assert denominator
        inverse = field.inverse[denominator]
        basis.append(tuple(field.multiply(c, inverse) for c in coefficients))
    return field, points, tuple(basis)


def crossing_coefficients_binary(x, y, u, v):
    h = sum(x)
    field, points, basis = binary_interpolation_data(h)
    coefficients = [0] * (4 * h + 1)
    for point, lagrange in zip(points, basis):
        value = crossing_value_binary(x, y, u, v, point, field)
        if value:
            for i, scalar in enumerate(lagrange):
                coefficients[i] ^= field.multiply(value, scalar)
    assert all(value in (0, 1) for value in coefficients), 'coefficients did not descend to F2'
    return tuple(coefficients)
