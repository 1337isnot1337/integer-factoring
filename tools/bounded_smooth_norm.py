"""Convert a bounded smooth represented norm into a smooth-form quotient.

The representation must be supplied. This is not an algorithm for finding
smooth values of a quadratic form; see research/bounded-smooth-norm.md.
"""
from math import gcd

from arithmetic_two_group import power
from binary_quadratic_forms import compose_forms, inverse_form, reduce_form


def smooth_factors(integer, bound):
    factors = []
    for candidate in range(2, bound+1):
        while integer % candidate == 0:
            integer //= candidate
            factors.append(candidate)
    return factors if integer == 1 else None


def bezout(left, right):
    a, b, x, y, u, v = left, right, 1, 0, 0, 1
    while b:
        quotient = a//b
        a, b, x, u, y, v = b, a-quotient*b, u, x-quotient*u, v, y-quotient*v
    return (x, y) if a == 1 else (-x, -y)


def balanced_divisor(factors, radius):
    """Find a two-bin norm split; exact divisor search is sometimes necessary."""
    loads = [1, 1]
    for prime in factors:
        side = 0 if loads[0] <= loads[1] else 1
        loads[side] *= prime
    if max(loads) <= radius:
        return loads[0]
    integer = loads[0]*loads[1]
    divisors = [1]
    for prime in sorted(set(factors)):
        powers = [prime**exponent for exponent in range(factors.count(prime)+1)]
        divisors = [value*multiplier for value in divisors for multiplier in powers if value*multiplier <= radius]
    return next((value for value in divisors if integer//value <= radius), None)


def norm_to_pair(form, x, y, smooth_bound, radius):
    a, b, c = form
    discriminant = b*b-4*a*c
    norm = a*x*x+b*x*y+c*y*y
    assert norm > 0 and smooth_bound >= 2
    assert 4*radius*radius < -discriminant
    assert norm <= radius*radius and gcd(norm, discriminant) == 1
    divisor = gcd(x, y)
    x, y = x//divisor, y//divisor
    primitive_norm = a*x*x+b*x*y+c*y*y
    factors = smooth_factors(primitive_norm, smooth_bound)
    assert factors is not None
    left_norm = balanced_divisor(factors, radius)
    if left_norm is None:
        return None
    v, negative_u = bezout(x, y)
    u = -negative_u
    assert x*v-y*u == 1
    middle = 2*a*x*u+b*(x*v+y*u)+2*c*y*v
    trailing = a*u*u+b*u*v+c*v*v
    assert middle*middle-4*primitive_norm*trailing == discriminant
    identity = reduce_form((1, discriminant % 2, ((discriminant % 2)-discriminant)//4))
    pair = [identity, identity]
    remainder = left_norm
    for prime in sorted(set(factors)):
        exponent = factors.count(prime)
        root = middle % (2*prime)
        assert (root*root-discriminant) % (4*prime) == 0
        value = reduce_form((prime, root, (root*root-discriminant)//(4*prime)))
        left_exponent = 0
        while remainder % prime == 0:
            remainder //= prime
            left_exponent += 1
        right_exponent = exponent-left_exponent
        pair[0] = compose_forms(pair[0], power(value, left_exponent, identity, compose_forms))
        pair[1] = compose_forms(pair[1], power(inverse_form(value), right_exponent, identity, compose_forms))
    assert remainder == 1
    assert [value[0] for value in pair] == [left_norm, primitive_norm//left_norm]
    assert compose_forms(pair[0], inverse_form(pair[1])) == reduce_form(form)
    return tuple(pair)
