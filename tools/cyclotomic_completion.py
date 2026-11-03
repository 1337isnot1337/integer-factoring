"""Recover bounded integral cyclotomic coefficients from one lifted field.

The original-index divisor coefficient must be supplied. This module supplies
Hensel acquisition and classical exact LLL/nearest-plane decoding, not a fast
coefficient engine. Costs are polynomial in the numerical root order and
precision. LLL reuses the repository's primitive-recipient implementation.
"""
from fractions import Fraction
from functools import cache
from odd_character_residues import cyclotomic_polynomial, root_of_unity_powers, reduce_cyclotomic
from cyclotomic_character_decoder import (binary_cyclotomic, modular_inverse,
    polynomial_product, polynomial_divmod, polynomial_remainder)
from check_factoring_primitive_recipients import lll_reduce, nearest
from check_carry_fourier_elimination import irreducible


def product(a, b):
    """Ordinary integer polynomial multiplication, before any reduction."""
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def bits(p, size):
    return tuple((p >> i) & 1 for i in range(size))


@cache
def hensel_factor(order, field, precision):
    """Acquire the chosen monic cyclotomic factor modulo 2**precision."""
    if not isinstance(order, int) or order < 3 or order % 2 == 0:
        raise ValueError('An odd root order at least three is required')
    if not isinstance(precision, int) or precision < 1:
        raise ValueError('Positive local precision is required')
    if not isinstance(field, int) or field < 4 or not irreducible(field):
        raise ValueError('A certified binary irreducible factor is required')
    phi = cyclotomic_polynomial(order)
    binary = binary_cyclotomic(order)
    companion, rem = polynomial_divmod(binary, field)
    if rem:
        raise ValueError('The field polynomial must divide the cyclotomic polynomial')
    degree, other = field.bit_length()-1, companion.bit_length()-1
    inverse = modular_inverse(companion, field)
    left, right = list(bits(field, degree+1)), list(bits(companion, other+1))
    for k in range(1, precision):
        value = product(left, right)
        assert all((a-b) % (1 << k) == 0 for a, b in zip(phi, value))
        error = sum((((a-b) >> k) & 1) << i for i, (a, b) in enumerate(zip(phi, value)))
        delta = polynomial_remainder(polynomial_product(error, inverse), field)
        quotient, rem = polynomial_divmod(error ^ polynomial_product(delta, companion), field)
        assert rem == 0 and quotient.bit_length() <= other
        for i in range(degree):
            left[i] += ((delta >> i) & 1) << k
        for i in range(other):
            right[i] += ((quotient >> i) & 1) << k
    assert all((a-b) % (1 << precision) == 0 for a, b in zip(phi, product(left, right)))
    return tuple(left)


def local_remainder(poly, factor, precision):
    return tuple(x % (1 << precision) for x in reduce_cyclotomic(poly, factor))


@cache
def completion_decoder(order, field, bound):
    """Acquire a decoder for sums of at most bound signed root powers."""
    if not isinstance(bound, int) or bound < 1:
        raise ValueError('A positive bound on the number of character terms is required')
    hensel_factor(order, field, 1)
    phi, powers = root_of_unity_powers(order)
    dimension, degree = len(phi)-1, field.bit_length()-1
    height = bound*max(abs(x) for row in powers for x in row)
    h = (height-1).bit_length()
    precision = dimension*(h+2*dimension+2)//degree+1
    factor = hensel_factor(order, field, precision)
    modulus = 1 << precision
    columns = [local_remainder((0,)*i+(1,), factor, precision) for i in range(dimension)]
    basis = []
    for i in range(dimension):
        row = [0]*dimension
        if i < degree:
            row[i] = modulus
        else:
            row[i] = 1
            for j, x in enumerate(columns[i]):
                row[j] = -x
        basis.append(row)
    reduced, stats = lll_reduce(basis)
    # Independently rebuild full Gram--Schmidt vectors for the decoder.
    orthogonal, norms = [], []
    for row in reduced:
        vector = [Fraction(x) for x in row]
        for previous, norm in zip(orthogonal, norms):
            scalar = sum(x*y for x, y in zip(row, previous))/norm
            vector = [x-scalar*y for x, y in zip(vector, previous)]
        orthogonal.append(vector)
        norms.append(sum(x*x for x in vector))
    assert all(4*height*height*dimension < norm for norm in norms)
    return factor, precision, height, reduced, orthogonal, norms, stats, basis


def reconstruct_coefficient(order, field, bound, observed):
    """Recover the promised bounded integral coefficient; no divisor search."""
    factor, precision, height, basis, orthogonal, norms, _, _ = completion_decoder(order, field, bound)
    degree, dimension = len(factor)-1, len(basis)
    if (len(observed) != degree or any(not isinstance(x, int) or not 0 <= x < 1 << precision
                                     for x in observed)):
        raise ValueError('Canonical local coordinates at the prescribed precision are required')
    initial = list(observed)+[0]*(dimension-degree)
    target = [-x for x in initial]
    for row, star, norm in reversed(list(zip(basis, orthogonal, norms))):
        coefficient = nearest(sum(x*y for x, y in zip(target, star))/norm)
        target = [x-coefficient*y for x, y in zip(target, row)]
    answer = tuple(-x for x in target)
    if max(abs(x) for x in answer) > height or local_remainder(answer, factor, precision) != tuple(observed):
        raise ValueError('The supplied residue has no bounded integral preimage')
    return answer
