"""Polynomial-time values of the normalized characteristic-two crossing pieces.

A polynomial over F_2 is encoded by the bits of a nonnegative Python integer.
For four distinct-odd-part labels x, v, u, y of weight h, final_piece_profile
returns the polynomial whose coefficient of Z^(s+2h) is A_s(x,y,u,v).

The formula folds the source dual polynomial into a polynomial Pfaffian.
It enumerates neither that dual polynomial nor intermediate partitions.
See research/modular-pfaffian.md for the ordinary proof and exact scope.
The underlying operators and crossing identity are Claude (Anthropic)'s
Appendix B, pinned in that note. This module is specific to characteristic 2.
"""

from functools import lru_cache


@lru_cache(maxsize=65536)
def polynomial_product(left, right):
    """Carryless multiplication in F_2[Z], with no numerical approximation."""
    assert left >= 0 and right >= 0
    if left.bit_count() > right.bit_count():
        left, right = right, left
    answer = 0
    while left:
        bit = left & -left
        answer ^= right << (bit.bit_length() - 1)
        left ^= bit
    return answer


def polynomial_exact_quotient(numerator, denominator):
    """The exact divisions required by fraction-free Pfaffian elimination."""
    assert numerator >= 0 and denominator > 0
    quotient = 0
    degree = denominator.bit_length()
    while numerator.bit_length() >= degree:
        shift = numerator.bit_length() - degree
        quotient ^= 1 << shift
        numerator ^= denominator << shift
    assert numerator == 0, "non-exact polynomial division"
    return quotient


@lru_cache(maxsize=4096)
def one_plus_z_power(exponent):
    """(1+Z)^exponent, using Frobenius for its binary digits."""
    assert exponent >= 0
    answer = 1
    bit_index = 0
    while exponent:
        if exponent & 1:
            answer ^= answer << (1 << bit_index)
        exponent >>= 1
        bit_index += 1
    return answer


def polynomial_pfaffian(matrix):
    """Pfaffian of an alternating matrix over F_2[Z].

    Simultaneous row/column swaps introduce no sign in characteristic two.
    The Pfaffian condensation identity makes every displayed division exact.
    A zero pivot row means the remaining Pfaffian is zero.
    """
    a = [row[:] for row in matrix]
    size = len(a)
    assert size % 2 == 0
    assert all(len(row) == size for row in a)
    assert all(a[i][i] == 0 for i in range(size))
    assert all(a[i][j] == a[j][i]
               for i in range(size) for j in range(i))
    if not size:
        return 1
    previous = 1
    for first in range(0, size - 2, 2):
        pivot = next((j for j in range(first + 1, size)
                      if a[first][j]), None)
        if pivot is None:
            return 0
        second = first + 1
        if pivot != second:
            a[pivot], a[second] = a[second], a[pivot]
            for row in a:
                row[pivot], row[second] = row[second], row[pivot]
        current = a[first][second]
        for j in range(first + 2, size):
            for k in range(j + 1, size):
                numerator = (
                    polynomial_product(current, a[j][k]) ^
                    polynomial_product(a[first][j], a[second][k]) ^
                    polynomial_product(a[first][k], a[second][j]))
                a[j][k] = a[k][j] = polynomial_exact_quotient(
                    numerator, previous)
        previous = current
    return a[size - 2][size - 1]


def folded_matrix(x, v, u, y):
    """The six folded blocks, in vertex order (U, A, B, Y)."""
    h = sum(x)
    assert sum(v) == sum(u) == sum(y) == h
    assert all(part > 0 and part % 2
               for label in (x, v, u, y) for part in label)
    assert all(len(set(label)) == len(label) for label in (x, v, u, y))
    nu, na, nb, ny = map(len, (u, x, v, y))
    size = nu + na + nb + ny
    assert size % 2 == 0
    matrix = [[0] * size for _ in range(size)]

    def put(i, j, polynomial):
        matrix[i][j] = matrix[j][i] = polynomial

    for ui, t in enumerate(u):
        for ai, n in enumerate(x):
            put(ui, nu + ai, 3 << (n | t))
        for bi, m in enumerate(v):
            if ((m - 1) & ~(t - 1)) == 0:
                put(ui, nu + na + bi, one_plus_z_power(m))
        for yi, t2 in enumerate(y):
            if t == t2:
                put(ui, nu + na + nb + yi, 1 << t2)
    for ai, n in enumerate(x):
        for bi, m in enumerate(v):
            if ((n - 1) & ~(m - 1)) == 0:
                put(nu + ai, nu + na + bi, 1 << n)
        for yi, t in enumerate(y):
            put(nu + ai, nu + na + nb + yi, 3 << ((n & t) - 1))
    for bi, m in enumerate(v):
        for yi, t in enumerate(y):
            if (m & t) == 1:
                put(nu + na + bi, nu + na + nb + yi,
                    one_plus_z_power(m) << t)
    return matrix


def final_piece_profile(x, v, u, y):
    """All A_s(x,y,u,v), stored as coefficient bits for -2h <= s <= 2h."""
    h = sum(x)
    answer = polynomial_pfaffian(folded_matrix(x, v, u, y))
    assert answer.bit_length() <= 4 * h + 1
    return answer
