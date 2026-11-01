"""A Gaussian-free inverse of the binary cyclotomic character map.

The input remains an unproved fast coefficient-engine interface. This module
only decodes supplied coefficients; it never obtains them from the index N.
Binary polynomials are integers whose bit i is the coefficient of T**i.
"""

from math import gcd, isqrt

from odd_character_residues import primitive_generator


def polynomial_product(left, right):
    result = 0
    while left:
        bit = left & -left
        result ^= right << (bit.bit_length()-1)
        left ^= bit
    return result


def polynomial_divmod(numerator, denominator):
    if not denominator:
        raise ZeroDivisionError("zero binary polynomial")
    quotient = 0
    while numerator.bit_length() >= denominator.bit_length():
        shift = numerator.bit_length()-denominator.bit_length()
        quotient ^= 1 << shift
        numerator ^= denominator << shift
    return quotient, numerator


def polynomial_remainder(numerator, denominator):
    return polynomial_divmod(numerator, denominator)[1]


def modular_product(left, right, modulus):
    return polynomial_remainder(polynomial_product(left, right), modulus)


def modular_power(value, exponent, modulus):
    result = 1
    while exponent:
        if exponent & 1:
            result = modular_product(result, value, modulus)
        value = modular_product(value, value, modulus)
        exponent >>= 1
    return result


def modular_inverse(value, modulus):
    previous, current = modulus, value
    previous_coefficient, coefficient = 0, 1
    while current:
        quotient, remainder = polynomial_divmod(previous, current)
        previous, current = current, remainder
        previous_coefficient, coefficient = (
            coefficient, previous_coefficient ^ polynomial_product(quotient, coefficient))
    if previous != 1:
        raise ValueError("binary polynomial is not a unit")
    return polynomial_remainder(previous_coefficient, modulus)


def modular_compose(polynomial, value, modulus):
    result = 0
    for index in range(polynomial.bit_length()-1, -1, -1):
        result = modular_product(result, value, modulus) ^ (polynomial >> index & 1)
    return result


def polynomial_bits(coefficients):
    return sum((value & 1) << index for index, value in enumerate(coefficients))


def binary_cyclotomic(order):
    """Construct Phi_order modulo two, using only divisors of order."""
    factors = {}
    for divisor in range(1, order+1):
        if order % divisor:
            continue
        value = (1 << divisor) ^ 1
        for smaller, factor in factors.items():
            if divisor % smaller == 0:
                value, remainder = polynomial_divmod(value, factor)
                if remainder:
                    raise ValueError("inexact cyclotomic division")
        factors[divisor] = value
    return factors[order]


def parameters(order):
    """order=2*e*u, with e a power of two and u odd; return Phi_u, Phi_order."""
    if not isinstance(order, int) or order < 2 or order % 2:
        raise ValueError("a positive even character order is required")
    e = (order & -order)//2
    u = order//(2*e)
    phi = binary_cyclotomic(u)
    modulus = sum(1 << (e*index) for index in range(phi.bit_length()) if phi >> index & 1)
    return e, u, phi, modulus


def unipotent_coefficients(value, e, u, phi, modulus):
    """Recover c_l from value=sum c_l T**(u*l), retaining every nilpotent jet."""
    complement, remainder = polynomial_divmod((1 << u) ^ 1, phi)
    if remainder:
        raise ValueError("invalid cyclotomic factor")
    inverse = modular_inverse(complement, modulus)
    precision, shifted = modulus, 0
    for level in range(e):
        constant = polynomial_remainder(value, phi)
        if constant not in (0, 1):
            raise ValueError("inconsistent character coefficient bits")
        shifted |= constant << level
        value, remainder = polynomial_divmod(value ^ constant, phi)
        if remainder:
            raise ValueError("inconsistent nilpotent coefficient")
        precision, remainder = polynomial_divmod(precision, phi)
        if remainder:
            raise ValueError("invalid nilpotent precision")
        value = modular_product(value, inverse, precision)
    if value:
        raise ValueError("inconsistent nilpotent coefficient bits")
    # If shifted encodes C(W), return C(V+1). Lucas' theorem makes this
    # the upper subset-zeta transform on the log2(e) index bits.
    step = 1
    while step < e:
        for index in range(e):
            if not index & step and shifted >> (index | step) & 1:
                shifted ^= 1 << index
        step <<= 1
    return shifted


def decode_character_polynomial(order, observations):
    """Invert h -> (h(T**j) mod Phi_order) for j dividing its odd part.

    h has degree <order/2. Observations map precisely those j to canonical
    binary polynomial masks. An inconsistent observation raises ValueError.
    No Gaussian elimination, field factorization or reduced quotient is used.
    """
    e, u, phi, modulus = parameters(order)
    width, dimension = order//2, modulus.bit_length()-1
    indices = tuple(index for index in range(1, u+1) if u % index == 0)
    if set(observations) != set(indices) or any(
            not isinstance(value, int) or value < 0 or value.bit_length() > dimension
            for value in observations.values()):
        raise ValueError("canonical observations for every divisor of the odd part are required")
    # T=omega*eta, omega**u=1, eta**e=1, eta**u=T**u.
    a = e*pow(e, -1, u) if u > 1 else 0
    b = (1-a) % width
    omega = modular_power(2, a, modulus)
    result = 0
    # Recompute evaluations for each s: this deliberately avoids an array of
    # u ring elements and keeps workspace linear in the supplied input bits.
    for s in range(u):
        fourier = 0
        for r in range(u):
            k = (a*r+b) % width
            j = gcd(k, u)
            # Lift k/j mod (e*u/j) to a unit modulo 2*e*u. Such a lift
            # exists by CRT, including r=0 and the unramified case e=1.
            step = width//j
            lift = k//j
            while gcd(lift, order) != 1:
                lift += step
            value = modular_compose(
                observations[j], modular_power(2, lift, modulus), modulus)
            weight = modular_power(omega, (-r*s) % u, modulus)
            fourier ^= modular_product(weight, value, modulus)
        fourier = modular_product(
            fourier, modular_power(2, (-b*s) % width, modulus), modulus)
        coefficients = unipotent_coefficients(fourier, e, u, phi, modulus)
        for level in range(e):
            result |= (coefficients >> level & 1) << (s+u*level)
    # Enforce consistency on arbitrary supplied observations, rather than
    # relying on them having been generated by a genuine divisor instance.
    for j, observed in observations.items():
        if modular_compose(result, modular_power(2, j, modulus), modulus) != observed:
            raise ValueError("inconsistent character coefficient bits")
    return result


def signed_residue_from_algebraic_character_bits(prime, coefficient):
    """Drop-in alternative to the Gaussian decoder on supplied M_j(N) mod 2."""
    if prime == 2:
        return (0, 0)
    if (not isinstance(prime, int) or prime < 3 or prime % 2 == 0
            or any(prime % divisor == 0 for divisor in range(3, isqrt(prime)+1, 2))):
        raise ValueError("an odd prime level is required")
    e, u, _, modulus = parameters(prime-1)
    generator = primitive_generator(prime)
    dimension = modulus.bit_length()-1
    observations = {}
    for index in range(1, u+1):
        if u % index:
            continue
        values = tuple(coefficient(generator, index))
        if len(values) != dimension or any(value not in (0, 1) for value in values):
            raise ValueError("canonical binary cyclotomic coordinate vectors are required")
        observations[index] = polynomial_bits(values)
    paired = decode_character_polynomial(prime-1, observations)
    output, residue = [0]*prime, 1
    for exponent in range(prime-1):
        output[residue] = paired >> (exponent % (e*u)) & 1
        residue = residue*generator % prime
    return tuple(output)
