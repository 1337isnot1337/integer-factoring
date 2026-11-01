"""Explicit unit-or-zero inverse acquisition for Teichmuller carry sums.

The Frobenius reindexing follows the credited field carry interface.
Norm coefficients use the lifted automorphism, not ordinary squaring.
The collected-denominator identity and inverse polynomial are credited in
research/precision-efficient-residue-indicators.md. This module acquires
their polynomial carry queries without expanding the denominator product.
"""
from teichmuller_cyclic_carries import parity_lift_coefficients, polynomial_sum


def unit_inverse_coefficients(precision):
    """Coefficients of the formal polynomial H_(h+1)(X)/X."""
    return parity_lift_coefficients(precision+1)[1:]


def conjugate_factor(ring, width, row_modulus, variables, factor, power):
    row, terms = factor
    t, columns = row_modulus, width//row_modulus
    transformed = []
    for coefficient, counts in terms:
        out = [0]*(variables*columns)
        for v in range(variables):
            for j in range(columns):
                exponent = (row+t*j+power) % width
                out[v*columns+exponent//t] += counts[v*columns+j]
        transformed.append((ring.frobenius(coefficient, power), tuple(out)))
    return ((row+power) % t, transformed)


def rational_sum(ring, width, row_modulus, variables, numerators, denominators, nonzero=False):
    """Sum products, replacing each nonunit denominator value by zero.

    Units contribute their actual ring inverses. All nonconstant denominator
    factors share at most precision+1 polynomial branches.
    """
    if width < 1 or ring.degree % width or row_modulus < 1 or width % row_modulus or variables < 1:
        raise ValueError('Incompatible ring, domain, or row dimensions')
    dimension = variables*(width//row_modulus)
    scalar, data = 1, []
    for factor in denominators:
        row, terms = factor
        if not 0 <= row < row_modulus or any(len(h) != dimension or min(h) < 0 for c, h in terms):
            raise ValueError('Bad inverse factor')
        if all(not c or not any(h) for c, h in terms):
            value = 0
            for c, h in terms:
                value = ring.add(value, c)
            if not ring.mod_two(value):
                return 0
            scalar = ring.mul(scalar, ring.unit_inverse(value))
            continue
        orbit = [conjugate_factor(ring, width, row_modulus, variables, factor, r)
                 for r in range(ring.degree)]
        data.append((orbit, orbit[1:]))
    choices = [(j, ring.integer(c)) for j, c in enumerate(unit_inverse_coefficients(ring.precision))
               if ring.integer(c)]
    answer = 0
    for exponent, c in choices if data else [(0, 1)]:
        coefficient, factors = ring.mul(scalar, c), list(numerators)
        for orbit, partial in data:
            factors.extend(partial)
            factors.extend(orbit*exponent)
        value = polynomial_sum(ring, width, row_modulus, variables, factors, nonzero)
        answer = ring.add(answer, ring.mul(coefficient, value))
    return answer
