"""Exact small-discriminant controls, not an optimized class-group package.

Composition searches for a coprime equivalent leading coefficient. This
finite search is intentionally explicit and has no polynomial-bit-cost claim.
"""
from math import gcd, isqrt


def reduce_form(form):
    a, b, c = form
    discriminant = b*b-4*a*c
    assert discriminant < 0 and a > 0 and gcd(gcd(a, b), c) == 1
    while True:
        shift = (a-b)//(2*a)
        c, b = c+b*shift+a*shift*shift, b+2*a*shift
        if a > c:
            a, b, c = c, -b, a
            continue
        if a == c and b < 0:
            b = -b
        assert b*b-4*a*c == discriminant and -a < b <= a <= c
        return a, b, c


def reduced_forms(discriminant):
    assert discriminant < 0 and discriminant % 4 in (0, 1)
    result = []
    for a in range(1, isqrt((-discriminant)//3)+1):
        for b in range(-a+1, a+1):
            numerator = b*b-discriminant
            if numerator % (4*a):
                continue
            c = numerator//(4*a)
            if a > c or (a == c and b < 0) or gcd(gcd(a, b), c) != 1:
                continue
            result.append((a, b, c))
    return result


def inverse_form(form):
    a, b, c = form
    return reduce_form((a, -b, c))


def compose_forms(left, right):
    a, b, c = left
    other_a, other_b, other_c = right
    discriminant = b*b-4*a*c
    assert other_b*other_b-4*other_a*other_c == discriminant
    # Q(1,k) can vanish identically modulo 2 only for the obstructed
    # even-leading, odd-middle, odd-trailing case. Swapping fixes it.
    if other_a % 2 == 0 and a % 2 == 0 and b % 2 == 1 and c % 2 == 1:
        a, b, c = c, -b, a
    shift = next(k for k in range(other_a) if gcd(a+b*k+c*k*k, other_a) == 1)
    a, b = a+b*shift+c*shift*shift, b+2*c*shift
    correction = ((other_b-b)//2)*pow(a, -1, other_a) % other_a
    middle = b+2*a*correction
    leading = a*other_a
    assert (middle*middle-discriminant) % (4*leading) == 0
    return reduce_form((leading, middle, (middle*middle-discriminant)//(4*leading)))


def ambiguous_factor(form, integer):
    a, b, c = form
    assert b == 0 or b == a or a == c
    candidates = (a, c) if b == 0 else ((a, 4*c-a) if b == a else (2*a-b, 2*a+b))
    return next((factor for value in candidates
                 for factor in [gcd(value, integer)] if 1 < factor < integer), None)


def smooth_prime_word(form, smooth_bound):
    """Return (prime form, signed exponent) factors, or None if not smooth."""
    a, b, c = form
    discriminant = b*b-4*a*c
    if gcd(a, discriminant) != 1:
        return None
    remaining, factors = a, []
    for prime in range(2, smooth_bound+1):
        if remaining % prime:
            continue
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        middle = next(value for value in range(prime+1)
                      if (value*value-discriminant) % (4*prime) == 0)
        sign = 1 if (b-middle) % (2*prime) == 0 else -1
        assert (b-sign*middle) % (2*prime) == 0
        prime_form = reduce_form((prime, middle, (middle*middle-discriminant)//(4*prime)))
        factors.append((prime_form, sign*exponent))
    return factors if remaining == 1 else None
