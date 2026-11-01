"""Acquire the coefficient-extension system packing from the precision note.

The supplied basic irreducible ring representation and its larger degree
are charged inputs. Equations test residue zeros, not full ring zeros.
"""
from teichmuller_cyclic_carries import parity_lift_coefficients, polynomial_sum


def packed_system_root_count(ring, width, variables, equations, nonzero=False):
    """Count simultaneous residue roots in the Teichmuller subfield.

    Coefficients are represented in ring, and must reduce into F_(2^width).
    A degree-D ring can pack at most D/width equations without adding
    summation variables or division precision. Exponents are ordinary.
    """
    if width < 1 or ring.degree % width or variables < 1:
        raise ValueError('Incompatible ring or domain dimensions')
    if len(equations) > ring.degree//width:
        raise ValueError('Insufficient relative coefficient-extension degree')
    combined, tag = [], 1
    for equation in equations:
        for coefficient, exponents in equation:
            if len(exponents) != variables or min(exponents) < 0:
                raise ValueError('Bad ordinary monomial')
            residue = ring.mod_two(coefficient)
            if ring.mod_two(ring.frobenius(residue, width)) != residue:
                raise ValueError('Equation coefficient outside the summation subfield')
            combined.append((ring.mul(tag, coefficient), tuple(exponents)))
        tag = ring.mul(tag, ring.generator)
    answer = polynomial_sum(ring, width, width, variables, [], nonzero)
    if not equations:
        return answer
    norm = [(r % width, [(ring.frobenius(c, r), exponents) for c, exponents in combined])
            for r in range(ring.degree)]
    for exponent, coefficient in enumerate(parity_lift_coefficients(ring.precision)):
        if coefficient:
            moment = polynomial_sum(ring, width, width, variables, norm*exponent, nonzero)
            answer = ring.sub(answer, ring.mul(ring.integer(coefficient), moment))
    assert not any(ring.coefficients(answer)[1:])
    return answer
