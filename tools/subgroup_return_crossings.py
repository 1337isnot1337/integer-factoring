#!/usr/bin/env python3
"""Ordered subgroup-return filters on the three-letter Steiner crossing."""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from check_steiner_cubic_crossings import curve, table
from weighted_smith_oracle import polynomial_multiply

FIRST = (3, 0, 1, 2)
SECOND = (3, 2, 0, 1)
PERMUTATIONS = {1: FIRST, 2: SECOND,
                -1: tuple(FIRST.index(i) for i in range(4)),
                -2: tuple(SECOND.index(i) for i in range(4))}
FAMILIES = ('free-commutator', 'four-state')


@lru_cache(None)
def reduce_word(word):
    stack = []
    for letter in word:
        if stack and stack[-1] == -letter:
            stack.pop()
        else:
            stack.append(letter)
    return tuple(stack)


def subgroup_return(word, family):
    if family == 'free-commutator':
        word = reduce_word(word)
        if len(word) % 4:
            return False
        power = len(word)//4
        return word == (1, 2, -1, -2)*power or word == (2, 1, -2, -1)*power
    if family == 'four-state':
        state = 0
        for letter in reversed(word):
            state = PERMUTATIONS[letter][state]
        return state == 0
    raise ValueError('Choose one of the two stated return filters.')


def local_rows():
    scalar = curve(table(3, 'affine'))
    rows = []
    for x, y in product(range(3), repeat=2):
        row = []
        for u, v in product(range(3), repeat=2):
            polynomial = scalar[3*x+y][3*u+v]
            if polynomial:
                word = reduce_word(tuple(g for g in (x, y, -u, -v) if g))
                row.append((3*u+v, polynomial, word))
        rows.append(row)
    return rows


polynomial_product = lru_cache(None)(polynomial_multiply)


def coefficient_rows(length, family):
    if length < 1 or family not in FAMILIES:
        raise ValueError('Use positive length and a stated return-filter family.')
    rows = [[(0, 1, ())]]
    for _ in range(length):
        rows = [[(9*c+d, polynomial_product(p, q), reduce_word(a+b))
                 for (c, p, a), (d, q, b) in product(left, right)]
                for left, right in product(rows, local_rows())]
    size = len(rows)
    cuts = [[[0]*size for _ in range(3*length+1)] for _ in range(2)]
    digits = [[(r//9**(length-1-i)) % 9 for i in range(length)] for r in range(size)]
    support = endpoints = 0
    for r, row in enumerate(rows):
        for c, polynomial, word in row:
            if not subgroup_return(word, family):
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
