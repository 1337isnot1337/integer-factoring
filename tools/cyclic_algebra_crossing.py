#!/usr/bin/env python3
"""Crossings from regular polynomial transforms of a binary full cycle.

Polynomials over F2 are integer bit masks. The constructor returns the two
polynomial cut matrices, with the complete degree n-1 reversal at infinity.
"""
from __future__ import annotations

from check_controlled_permutation_smith import binomial_polynomial, two_control_table


def cyclic_multiply(a: int, b: int, n: int) -> int:
    """Multiply modulo X**n + 1; both inputs have degree less than n."""
    mask, answer = (1 << n) - 1, 0
    while b:
        bit = b & -b
        shift = bit.bit_length() - 1
        answer ^= ((a << shift) & mask) ^ (a >> (n - shift) if shift else 0)
        b ^= bit
    return answer


def reverse(value: int, length: int) -> int:
    return int(f'{value:0{length}b}'[::-1], 2)


def default_cycle(q: int) -> tuple[int, ...]:
    if q < 2 or q & (q - 1):
        raise ValueError('the alphabet must be a power of two, at least two')
    change = list(range(q))
    change[0], change[1] = change[1], change[0]
    control = sum(1 << i for i in range(1, q // 2))
    table = two_control_table(q, control, change)
    return tuple(q * j + table[j][i] for i in range(q) for j in range(q))


def construct(q: int, g: int, cycle=None):
    """Return row/column K(z) and row H(w)=w**(n-1)K(1/w).

    Here U0 is the supplied full cycle, S0=I+U0, U=I+g(S0), and
    K(z)=Phi(sum_j z**(n-1-j) g(S0)**j U). Require g(0)=0, g'(0)=1.
    """
    default = default_cycle(q)
    n = q * q
    if not isinstance(g, int) or g < 0 or g.bit_length() > n or g & 3 != 2:
        raise ValueError('g must have degree below q^2 and begin with S')
    cycle = tuple(default if cycle is None else cycle)
    if len(cycle) != n or sorted(cycle) != list(range(n)):
        raise ValueError('supply a permutation of q^2 coordinates')
    orbit, position, current = [], {}, 0
    while current not in position:
        position[current] = len(orbit)
        orbit.append(current)
        current = cycle[current]
    if len(orbit) != n:
        raise ValueError('the permutation must be one full cycle')

    # Translate the truncated S0 basis into F2[U0]/(U0**n+1).
    s = 0
    for degree in range(n):
        if g >> degree & 1:
            s ^= binomial_polynomial(degree)
    power = s ^ 1
    entries = [0] * n
    for degree in range(n):
        bits = power
        while bits:
            bit = bits & -bits
            entries[bit.bit_length() - 1] ^= 1 << degree
            bits ^= bit
        power = cyclic_multiply(power, s, n)
    if power:
        raise ArithmeticError('nilpotency failed')

    infinity = [[entries[(position[q*x+u] - position[q*v+y]) % n]
                 for u in range(q) for v in range(q)]
                for x in range(q) for y in range(q)]
    row = [[reverse(value, n) for value in line] for line in infinity]
    column = [[row[q*x+y][q*u+v] for y in range(q) for v in range(q)]
              for x in range(q) for u in range(q)]
    return row, column, infinity


def coefficient_rows(matrix, degree: int):
    return [sum(((value >> degree) & 1) << j for j, value in enumerate(row))
            for row in matrix]
