"""Exact odd-character Fourier conversion for weight-one divisor coefficients.

The coefficient provider is an unproved fast-algorithm input. Cyclotomic
arithmetic and conversion to signed cyclic divisor residues are explicit.
"""

from functools import lru_cache


def exact_polynomial_quotient(numerator, denominator):
    output = list(numerator)
    quotient = [0]*(len(numerator)-len(denominator)+1)
    for offset in range(len(quotient)-1, -1, -1):
        coefficient = output[offset+len(denominator)-1]
        quotient[offset] = coefficient
        for i, value in enumerate(denominator):
            output[offset+i] -= coefficient*value
    assert not any(output)
    return tuple(quotient)


@lru_cache(None)
def cyclotomic_polynomial(order):
    polynomial = (-1,)+(0,)*(order-1)+(1,)
    for divisor in range(1, order):
        if order % divisor == 0:
            polynomial = exact_polynomial_quotient(polynomial, cyclotomic_polynomial(divisor))
    return polynomial


def reduce_cyclotomic(polynomial, modulus):
    output = list(polynomial)+[0]*max(0, len(modulus)-1-len(polynomial))
    for position in range(len(output)-1, len(modulus)-2, -1):
        coefficient = output[position]
        offset = position-len(modulus)+1
        for i, value in enumerate(modulus):
            output[offset+i] -= coefficient*value
    return tuple(output[:len(modulus)-1])


def multiply_cyclotomic(left, right, modulus):
    output = [0]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            output[i+j] += a*b
    return reduce_cyclotomic(output, modulus)


@lru_cache(None)
def root_of_unity_powers(order):
    modulus = cyclotomic_polynomial(order)
    powers = []
    current = (1,)+(0,)*(len(modulus)-2)
    for _ in range(order):
        powers.append(current)
        current = reduce_cyclotomic((0,)+current, modulus)
    assert current == powers[0]
    return modulus, tuple(powers)


def primitive_generator(prime):
    remaining, factors, candidate = prime-1, [], 2
    while candidate*candidate <= remaining:
        if remaining % candidate == 0:
            factors.append(candidate)
            while remaining % candidate == 0:
                remaining //= candidate
        candidate += 1
    if remaining > 1:
        factors.append(remaining)
    return next(value for value in range(2, prime)
                if all(pow(value, (prime-1)//factor, prime) != 1 for factor in factors))


def signed_residue_from_odd_characters(prime, coefficient, characteristic=2):
    """Convert c_j=sum_(d|N) chi_j(d) for odd j to c_a-c_-a.

    Characters are extended by zero on multiples of the prime. coefficient(g,j)
    returns an integer vector in Q[zeta_(q-1)]'s standard power basis. q=2 has the
    zero signed image. A primitive generator is found by small trial work.
    """
    if prime == 2:
        return (0, 0)
    if prime < 3 or characteristic < 2:
        raise ValueError("an odd prime level and characteristic>=2 are required")
    generator = primitive_generator(prime)
    order = prime-1
    modulus, powers = root_of_unity_powers(order)
    values = {index: tuple(coefficient(generator, index)) for index in range(1, order, 2)}
    if any(len(value) != len(modulus)-1 or any(not isinstance(x, int) for x in value)
           for value in values.values()):
        raise ValueError("exact integral cyclotomic coefficient vectors are required")
    output, residue = [0]*prime, 1
    for exponent in range(order):
        total = [0]*(len(modulus)-1)
        for index, value in values.items():
            summand = multiply_cyclotomic(value, powers[-index*exponent % order], modulus)
            total = [left+right for left, right in zip(total, summand)]
        if any(total[1:]) or (2*total[0]) % order:
            raise ValueError("the odd-character inverse transform is not an integer residue count")
        output[residue] = (2*total[0]//order) % characteristic
        residue = residue*generator % prime
    return tuple(output)
