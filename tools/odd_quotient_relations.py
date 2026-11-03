"""Odd-residual relation certificates; no smooth-pair search or class number."""
from math import isqrt, lcm

from arithmetic_two_group import power
from check_positive_relation_kernel import determinant


def word_value(generators, row, identity, multiply, inverse):
    value = identity
    for generator, exponent in zip(generators, row):
        base = generator if exponent >= 0 else inverse(generator)
        value = multiply(value, power(base, abs(exponent), identity, multiply))
    return value


def odd_relation_annihilator(generators, rows, odd_witnesses,
                             identity, multiply, inverse):
    """Validate row residues of certified odd order, then return an annihilator.

    If R is nonsingular and each row residue raised to its supplied positive
    odd witness is identity, lcm(witnesses)*abs(det(R)) kills the generated
    finite abelian group. No coordinates or factorization of its order enter.
    """
    size = len(generators)
    if (len(rows) != size or any(len(row) != size for row in rows)
            or len(odd_witnesses) != size
            or any(not isinstance(x, int) for row in rows for x in row)):
        raise ValueError('A square integer matrix and one odd witness per row are required')
    if any(not isinstance(x, int) or x <= 0 or x % 2 == 0 for x in odd_witnesses):
        raise ValueError('Positive odd row-annihilator witnesses are required')
    det = determinant(rows)
    if not det:
        raise ValueError('The relation matrix must be nonsingular')
    for row, witness in zip(rows, odd_witnesses):
        residue = word_value(generators, row, identity, multiply, inverse)
        if power(residue, witness, identity, multiply) != identity:
            raise ValueError('A supplied odd witness does not kill its row residue')
    return abs(det)*lcm(*odd_witnesses)


def uniform_odd_multiplier(order_bit_bound, smooth_bound):
    """Kill every smooth odd order <=2**order_bit_bound, by charged trial search.

    Complexity is polynomial in the NUMERICAL smooth_bound and the bit bound.
    This constructs no class number and no smooth-pair extraction procedure.
    """
    if (not isinstance(order_bit_bound, int) or order_bit_bound < 0
            or not isinstance(smooth_bound, int) or smooth_bound < 0):
        raise ValueError('Nonnegative integer bounds are required')
    limit, result = 1 << order_bit_bound, 1
    for prime in range(3, min(smooth_bound, limit)+1, 2):
        if any(prime % d == 0 for d in range(3, isqrt(prime)+1, 2)):
            continue
        prime_power = prime
        while prime_power <= limit//prime:
            prime_power *= prime
        result *= prime_power
    return result
