#!/usr/bin/env python3
"""Integer charge and area normal forms for growing qutrit Weyl crossings."""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from check_fixed_field_descent import Field
from check_steiner_cubic_crossings import curve, table
from weighted_smith_oracle import polynomial_multiply


def remainder(a, b):
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length()-b.bit_length())
    return a


def gcd(a, b):
    while b:
        a, b = b, remainder(a, b)
    return a


def irreducible(modulus):
    degree = modulus.bit_length()-1
    field = Field(modulus)
    x = 2
    for i in range(1, degree+1):
        x = field.mul(x, x)
        if i <= degree//2 and gcd(x ^ 2, modulus) != 1:
            return False
    return x == 2


@lru_cache(None)
def clock(length):
    if length < 1:
        raise ValueError('Use a positive physical word length.')
    degree = (2*length+1).bit_length()
    modulus = next((1 << degree) | low for low in range(1, 1 << degree, 2)
                   if irreducible((1 << degree) | low))
    field = Field(modulus)
    order = field.order-1
    n, divisor, primes = order, 2, []
    while divisor*divisor <= n:
        if n % divisor == 0:
            primes.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        primes.append(n)
    root = next(a for a in range(2, field.order)
                if all(field.power(a, order//p) != 1 for p in primes))
    powers = tuple(field.power(root, k) for k in range(order))
    assert len(set(powers)) == order > 2*length
    return field, powers


def statistics(physical):
    charge = area = 0
    for x, y, u, v in physical:
        a, c = x-v, x+y-u-v
        area += charge*a-v*c
        charge += c
    return charge, area


def local_rows():
    scalar = curve(table(3, 'affine'))
    return [[(3*u+v, scalar[3*x+y][3*u+v], x-v, x+y-u-v, v)
             for u, v in product(range(3), repeat=2) if scalar[3*x+y][3*u+v]]
            for x, y in product(range(3), repeat=2)]


polynomial_product = lru_cache(None)(polynomial_multiply)


def area_zero_cuts(length):
    if length < 1:
        raise ValueError('Use a positive physical word length.')
    rows = [[(0, 1, 0, 0)]]
    for _ in range(length):
        rows = [[(9*c+d, polynomial_product(p, q), charge+dc, area+charge*a-v*dc)
                 for (c, p, charge, area), (d, q, a, dc, v) in product(left, right)]
                for left, right in product(rows, local_rows())]
    size = len(rows)
    cuts = [[[0]*size for _ in range(3*length+1)] for _ in range(2)]
    digits = [[(r//9**(length-1-i)) % 9 for i in range(length)] for r in range(size)]
    support = endpoints = 0
    for r, row in enumerate(rows):
        for c, polynomial, charge, area in row:
            if charge or area:
                continue
            new_row = new_column = swap = 0
            for first, second in zip(digits[r], digits[c]):
                x, y = divmod(first, 3)
                u, v = divmod(second, 3)
                new_row = 9*new_row+3*x+u
                new_column = 9*new_column+3*y+v
                swap = 9*swap+3*y+x
            endpoint = polynomial.bit_count() % 2
            assert endpoint == int(c == swap)
            endpoints += endpoint
            support += 1
            while polynomial:
                term = polynomial & -polynomial
                k = term.bit_length()-1
                cuts[0][k][r] ^= 1 << c
                cuts[1][k][new_row] ^= 1 << new_column
                polynomial ^= term
    assert endpoints == size
    assert cuts[0][0] == [1 << r for r in range(size)]
    return cuts, support
