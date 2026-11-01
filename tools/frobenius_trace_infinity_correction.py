#!/usr/bin/env python3
"""Polynomial-cost Frobenius rational traces with arbitrary formal degree."""
from __future__ import annotations

from check_frobenius_rational_trace import coefficient_carries, rational_value
from fermat_torus_binomial_contraction import carry_coefficient


def zero_safe_binomial_trace(w, degree, scalar, numerators, denominators, *,
                             multiply, add=lambda a, b: a ^ b, zero=0, one=1):
    """Fixed division-free trace circuit, with field-valued iota(0)=0.

    Over a coefficient algebra this evaluates the polynomial power
    x^(2^(degree+1)-3), rather than asserting an inverse for zero divisors.
    """
    if (not isinstance(w, int) or w < 1 or not isinstance(degree, int)
            or degree < 1 or degree % w):
        raise ValueError("Use a positive field degree divisible by the subfield degree.")
    numerators, denominators = list(numerators), list(denominators)
    if any(not (0 <= exponent < w)
           for _, _, exponent in numerators + denominators):
        raise ValueError("Use Frobenius exponents in the promised range.")
    factors = list(numerators)
    for constant, slope, exponent in denominators:
        a, b = constant, slope
        for bit in range(degree + 1):
            if bit == 0 or bit >= 2:
                factors.append((a, b, (exponent + bit) % w))
            a, b = multiply(a, a), multiply(b, b)
    assert len(factors) == len(numerators) + degree * len(denominators)
    groups, span, constant = {}, 0, scalar
    for a, b, exponent in factors:
        span += 1 << exponent
        constant = multiply(constant, a)
        previous = groups.get(exponent, [one])
        current = [zero] * (len(previous) + 1)
        for index, value in enumerate(previous):
            current[index] = add(current[index], multiply(a, value))
            current[index + 1] = add(current[index + 1], multiply(b, value))
        groups[exponent] = current
    full = zero
    bits = max(1, w, span.bit_length())
    for target in range((1 << w) - 1, span + 1, (1 << w) - 1):
        full = add(full, carry_coefficient(
            groups, target, len(factors), bits, multiply, add, zero, one))
    full = multiply(scalar, full)
    return full, add(full, constant)


def degree_free_trace(field, w, scalar, numerators, denominators, statistics=None):
    """Return full-field and nonzero-field sums, with zero-free denominators."""
    if not isinstance(w, int) or w < 1 or field.degree % w:
        raise ValueError("The explicit field must contain the requested subfield.")
    if any(not (0 <= exponent < w)
           for _, _, exponent in numerators + denominators):
        raise ValueError("Use Frobenius exponents in the promised range.")
    counts = statistics if statistics is not None else {}
    def count(name):
        counts[name] = counts.get(name, 0) + 1

    base = 1 << w
    roots, grouped = [], {}
    for constant, slope, exponent in denominators:
        if not slope:
            raise ValueError("Denominator slopes must be nonzero.")
        ratio = field.mul(constant, field.inverse(slope))
        root = field.power(ratio, 1 << ((-exponent) % field.degree))
        if field.power(root, base) == root:
            raise ValueError("A denominator root lies in the summation field.")
        roots.append(root)
        grouped[root] = grouped.get(root, 0) + (1 << exponent)

    finite = 0
    for root, multiplicity in grouped.items():
        shifted_numerators = [
            (a ^ field.mul(b, field.power(root, 1 << exponent)), b, exponent)
            for a, b, exponent in numerators
        ]
        shifted_denominators = []
        shifted_scalar = scalar
        for (constant, slope, exponent), other_root in zip(denominators, roots):
            if other_root == root:
                shifted_scalar = field.mul(shifted_scalar, field.inverse(slope))
            else:
                value = constant ^ field.mul(slope, field.power(root, 1 << exponent))
                inverse = field.inverse(value)
                shifted_scalar = field.mul(shifted_scalar, inverse)
                shifted_denominators.append((field.mul(slope, inverse), exponent))
        delta = field.power(root, base) ^ root
        assert delta
        for offset in range((multiplicity - 1) // base + 1):
            factors = list(shifted_denominators)
            term_scalar = shifted_scalar
            for bit in range((offset + 1).bit_length()):
                if (offset + 1) >> bit & 1:
                    value = field.power(delta, 1 << bit)
                    inverse = field.inverse(value)
                    term_scalar = field.mul(term_scalar, inverse)
                    factors.append((inverse, bit))
            target = multiplicity - 1 - base * offset
            finite ^= field.mul(term_scalar, coefficient_carries(
                field, shifted_numerators, factors, target))
            count("finite_pole_coefficient_queries")

    # Reversal gives the regular series H_infinity=Z^(N-D)R(1/Z).
    formal_degree = (sum(1 << exponent for _, _, exponent in numerators)
                     - sum(1 << exponent for _, _, exponent in denominators))
    infinity = 0
    if formal_degree >= base - 1:
        reversed_numerators = [(b, a, exponent) for a, b, exponent in numerators]
        reversed_denominators = []
        reversed_scalar = scalar
        for constant, slope, exponent in denominators:
            inverse = field.inverse(slope)
            reversed_scalar = field.mul(reversed_scalar, inverse)
            reversed_denominators.append((field.mul(constant, inverse), exponent))
        for target in range(formal_degree - (base - 1), -1, -(base - 1)):
            infinity ^= field.mul(reversed_scalar, coefficient_carries(
                field, reversed_numerators, reversed_denominators, target))
            count("infinity_coefficient_queries")
    if infinity:
        count("nonzero_infinity_corrections")
    result = finite ^ infinity
    return result, result ^ rational_value(field, scalar, numerators, denominators, 0)
