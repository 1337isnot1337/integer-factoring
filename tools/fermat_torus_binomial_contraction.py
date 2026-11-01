#!/usr/bin/env python3
"""Division-free contraction of Frobenius binomials on a Fermat torus."""
from __future__ import annotations


def grouped_choices(m, factors, multiply, add, zero, one):
    groups = {}
    negative_shift = span = 0
    for a, b, exponent in factors:
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Use nonnegative integer Frobenius exponents.")
        reduced = exponent % (2 * m)
        if reduced >= m:
            reduced -= m
            a, b = b, a
            negative_shift += 1 << reduced
        span += 1 << reduced
        previous = groups.get(reduced, [one])
        current = [zero] * (len(previous) + 1)
        for count, value in enumerate(previous):
            current[count] = add(current[count], multiply(a, value))
            current[count + 1] = add(current[count + 1], multiply(b, value))
        groups[reduced] = current
    return groups, negative_shift, span


def carry_coefficient(groups, target, factor_count, bits, multiply, add, zero, one):
    """Extract one coefficient with a fixed, value-independent carry circuit."""
    if target < 0:
        return zero
    state = [zero] * (factor_count + 1)
    state[0] = one
    for bit in range(bits):
        local = groups.get(bit, [one])
        updated = [zero] * (factor_count + 1)
        requested = (target >> bit) & 1
        for carry, weight in enumerate(state):
            for count, coefficient in enumerate(local):
                total = carry + count
                if (total & 1) == requested:
                    next_carry = total >> 1
                    assert next_carry <= factor_count
                    updated[next_carry] = add(
                        updated[next_carry], multiply(weight, coefficient))
        state = updated
    return state[0]


def torus_binomial_sum(m, factors, *, multiply, add=lambda a, b: a ^ b,
                      zero=0, one=1, monomial=0, scalar=None):
    """Sum c*v^h*product(a+b*v^(2^e)) over the scalar Mth roots.

    Coefficients may belong to any commutative algebra of characteristic
    two. Its coefficient field supplies the M=2^m+1 scalar roots. This routine
    performs only the supplied additions and multiplications.
    """
    if not isinstance(m, int) or m < 1:
        raise ValueError("Use a positive integer torus parameter.")
    if not isinstance(monomial, int):
        raise ValueError("Use an integer monomial exponent.")
    factors = list(factors)
    if scalar is None:
        scalar = one
    modulus = (1 << m) + 1
    groups, negative_shift, span = grouped_choices(
        m, factors, multiply, add, zero, one)
    shift = (monomial - negative_shift) % modulus
    first = (-shift) % modulus
    bits = max(1, m, span.bit_length())
    result = zero
    queries = 0
    for target in range(first, span + 1, modulus):
        result = add(result, carry_coefficient(
            groups, target, len(factors), bits, multiply, add, zero, one))
        queries += 1
    assert queries <= len(factors) + 1
    return multiply(scalar, result)


def dyadic_box_star_sum(m, boxes, *, multiply, add=lambda a, b: a ^ b,
                       zero=0, one=1, monomial=0, scalar=None):
    """Contract product of S_(2^width)(gamma*v^(2^exponent))."""
    factors = []
    for gamma, exponent, width in boxes:
        if not isinstance(width, int) or width < 0:
            raise ValueError("Use nonnegative integer box widths.")
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Use nonnegative integer Frobenius exponents.")
        coefficient = gamma
        for bit in range(width):
            factors.append((one, coefficient, exponent + bit))
            coefficient = multiply(coefficient, coefficient)
    return torus_binomial_sum(
        m, factors, multiply=multiply, add=add, zero=zero, one=one,
        monomial=monomial, scalar=scalar)
