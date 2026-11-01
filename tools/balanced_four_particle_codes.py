#!/usr/bin/env python3
"""Exact strict Gaussian codes on arbitrary four-particle states of eight modes.

The five-state reflection alphabet is only one choice among the seventy
four-particle states.  This constructor permits all of them on every boundary.
It retains the complete polynomial tensor, removes its exact common factor,
and checks the strict endpoint after normalization.
"""
from __future__ import annotations

from itertools import combinations

from saturated_dual_kernel_directions import divide, polynomial_gcd
from selected_gaussian_coefficients import selected_tensor


REFERENCE_CODE = (15, 77, 43, 23, 51)
FOUR_PARTICLE_STATES = tuple(sum(1 << i for i in positions)
                             for positions in combinations(range(8), 4))


def grade(mask):
    return 2*sum(i for i in range(8) if (mask >> i) & 1)+4


def canonical_code(code):
    code = tuple(code)
    if len(code) != 5 or len(set(code)) != 5 or any(mask not in FOUR_PARTICLE_STATES for mask in code):
        raise ValueError('Choose five distinct four-particle masks on eight modes.')
    return tuple(sorted(code, key=lambda mask: (grade(mask), mask)))


def natural_energies(code):
    grades = tuple(grade(mask) for mask in code)
    return tuple(value-min(grades) for value in grades)


def compress_square(value):
    """Return f with input f(w^2); reject any odd exponent."""
    result = 0
    while value:
        bit = value & -value
        value ^= bit
        exponent = bit.bit_length()-1
        if exponent & 1:
            raise ArithmeticError('A balanced four-particle coefficient has an odd exponent.')
        result ^= 1 << (exponent//2)
    return result


def construct(code):
    code = canonical_code(code)
    original = selected_tensor(tuple(range(8)), 8, code)
    compressed = []
    for low, high in original['values']:
        assert not high
        compressed.append(compress_square(low))
    common = 0
    for value in compressed:
        common = polynomial_gcd(common, value)
    assert common and common.bit_count() % 2 == 1
    normalized = []
    for index, value in enumerate(compressed):
        quotient, remainder = divide(value, common)
        assert not remainder
        x, y, u, v = ((index//5**i) % 5 for i in range(4))
        assert quotient.bit_count() % 2 == int(x == v and y == u)
        normalized.append((quotient, 0))
    return {'values': tuple(normalized), 'code_masks': code, 'dimensions': [5]*4,
            'paired_energies': natural_energies(code), 'entry_gcd_in_z': common,
            'entry_gcd_degree': common.bit_length()-1,
            'degree': max(value.bit_length()-1 for value, high in normalized),
            'raw_degree_in_z': max(value.bit_length()-1 for value in compressed),
            'strict_endpoint': True, 'ramification_w_equals_square_root_z': 2,
            'native_complete_degree_bound_in_w': original['coefficient_degree_bound'],
            'native_exact_precision': original['exact_precision']}

