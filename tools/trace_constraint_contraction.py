"""Affine binary trace constraints outside a synchronized product.

Constraints are (coefficients, bit), meaning sum_i Tr(a_i U_i) = bit.
The coefficients and variables lie in the width-bit subfield of field.
"""
from itertools import product

from joint_frobenius_carries import joint_trace


def independent_constraints(field, width, variables, constraints):
    """Return an equivalent independent system, or None if inconsistent."""
    basis = {}
    mask = field.size - 1
    for coefficients, bit in constraints:
        if len(coefficients) != variables or bit not in (0, 1):
            raise ValueError('Bad trace constraint')
        packed = 0
        for i, a in enumerate(coefficients):
            if not 0 <= a < field.size or field.power(a, 1 << width) != a:
                raise ValueError('Constraint coefficient outside subfield')
            packed |= a << (i * field.degree)
        while packed:
            pivot = packed.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = (packed, bit)
                break
            old, rhs = basis[pivot]
            packed ^= old
            bit ^= rhs
        if not packed and bit:
            return None
    return [
        (tuple((row >> (i * field.degree)) & mask for i in range(variables)), bit)
        for row, bit in basis.values()
    ]


def trace_packets(field, width, t, coefficients, bit):
    """Split the 0/1 indicator into t synchronized linear polynomials."""
    D = width // t
    dimension = len(coefficients) * D
    packets = []
    for row in range(t):
        terms = []
        if row == 0 and bit == 0:
            terms.append((1, (0,) * dimension))
        for i, a in enumerate(coefficients):
            if not a:
                continue
            for j in range(D):
                counts = [0] * dimension
                counts[i * D + j] = 1
                terms.append((field.power(a, 1 << (row + t * j)), tuple(counts)))
        if terms:
            packets.append((row, terms))
    return packets


def constrained_trace(field, width, t, variables, numerators, denominators,
                      constraints, nonzero=False):
    """Return the constrained sum and exact expansion statistics."""
    if width < 1 or t < 1 or field.degree % width or width % t or variables < 1:
        raise ValueError('Incompatible dimensions')
    independent = independent_constraints(field, width, variables, constraints)
    if independent is None:
        return 0, {'inconsistent': True, 'rank': None, 'summands': 0}
    packets = [trace_packets(field, width, t, a, bit) for a, bit in independent]
    answer = 0
    summands = 0
    for selected in product(*packets):
        answer ^= joint_trace(field, field.degree, width, t, variables,
                              list(numerators) + list(selected), denominators, nonzero)
        summands += 1
    return answer, {'inconsistent': False, 'rank': len(independent), 'summands': summands}
