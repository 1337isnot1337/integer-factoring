"""Polynomial supplied-norm conversion with an inner/outer radius gap.

The factorization or a smoothness cutoff must be supplied.  The module does
not find represented smooth norms and does not decide exact-radius balance.
See research/smooth-norm-gap.md for the all-input guarantee and bit cost.
"""

from math import gcd, isqrt, prod

from bounded_smooth_norm import bezout, smooth_factors


def trimmed_lower_divisor(factors, inverse_slack=1, *, with_trace=False):
    """Approximate the largest divisor at most sqrt(n) from prime factors.

    For z=inverse_slack >= 1, the answer d obeys
        d <= d_opt <= (1 + 1/z) * d.
    Comparisons and trimming use integers only.  A nonprime item list instead
    gives the same guarantee over its subset products, not all divisors.
    """
    factors = tuple(factors)
    if not isinstance(inverse_slack, int) or inverse_slack < 1:
        raise ValueError("inverse_slack must be a positive integer")
    if any(not isinstance(value, int) or value < 2 for value in factors):
        raise ValueError("factors must be integers at least two")
    integer = prod(factors)
    target = isqrt(integer)
    if not factors:
        result = {"divisor": 1, "integer": 1, "target": 1,
                  "retained_sizes": [1], "trim_denominator": 0,
                  "candidate_count": 0}
        return result if with_trace else 1
    denominator = 2 * len(factors) * inverse_slack
    retained = [1]
    sizes = [1]
    candidate_count = 0
    for factor in factors:
        candidates = sorted(set(retained + [value * factor for value in retained
                                            if value * factor <= target]))
        candidate_count += len(candidates)
        trimmed = []
        for value in candidates:
            if not trimmed or denominator * value > (denominator + 1) * trimmed[-1]:
                trimmed.append(value)
        retained = trimmed
        sizes.append(len(retained))
    result = {"divisor": retained[-1], "integer": integer, "target": target,
              "retained_sizes": sizes, "trim_denominator": denominator,
              "candidate_count": candidate_count}
    return result if with_trace else result["divisor"]


def balanced_divisor_with_gap(factors, outer_radius, inverse_slack=1):
    """Return a valid outer-radius split, or None without an exact negative claim.

    Every split with both parts <= R_inner is accepted if
        (z + 1) * R_inner <= z * outer_radius.
    """
    if not isinstance(outer_radius, int) or outer_radius < 1:
        raise ValueError("outer_radius must be a positive integer")
    factors = tuple(factors)
    integer = prod(factors)
    divisor = trimmed_lower_divisor(factors, inverse_slack)
    if divisor <= outer_radius and integer // divisor <= outer_radius:
        return divisor
    return None


def pair_from_primitive_divisor(form, x, y, divisor, radius):
    """Construct the pair directly from one balanced primitive represented norm.

    This uses no class-group composition or search for modular square roots.
    The input vector is primitive, the norm is coprime to the discriminant,
    and divisor and norm/divisor are both at most radius < sqrt(|Delta|)/2.
    """
    a, b, c = form
    discriminant = b*b - 4*a*c
    norm = a*x*x + b*x*y + c*y*y
    if (discriminant >= 0 or a <= 0 or gcd(gcd(a, b), c) != 1
            or gcd(x, y) != 1 or norm <= 0
            or gcd(norm, discriminant) != 1 or divisor <= 0
            or norm % divisor or max(divisor, norm // divisor) > radius
            or 4*radius*radius >= -discriminant):
        raise ValueError("invalid primitive balanced-norm input")
    v, negative_u = bezout(x, y)
    u = -negative_u
    assert x*v - y*u == 1
    middle = 2*a*x*u + b*(x*v + y*u) + 2*c*y*v
    pair = []
    for leading, root in ((divisor, middle), (norm // divisor, -middle)):
        centered = root + 2*leading*((leading-root)//(2*leading))
        numerator = centered*centered - discriminant
        assert numerator % (4*leading) == 0
        trailing = numerator // (4*leading)
        assert -leading < centered <= leading < trailing
        assert gcd(gcd(leading, centered), trailing) == 1
        pair.append((leading, centered, trailing))
    return tuple(pair)


def norm_to_pair_with_gap(form, x, y, smooth_bound, radius, inverse_slack=1):
    """Convert a supplied smooth norm when the trimmed split fits the radius.

    None means that this approximation did not supply a split.  It is not a
    certificate that no exact-radius split or smooth quotient exists.
    """
    a, b, c = form
    discriminant = b*b - 4*a*c
    if discriminant >= 0 or a <= 0 or smooth_bound < 2 or radius < 1:
        raise ValueError("invalid form, cutoff, or radius")
    common = gcd(x, y)
    if common == 0:
        raise ValueError("the zero vector is not a norm witness")
    x, y = x // common, y // common
    norm = a*x*x + b*x*y + c*y*y
    if norm <= 0 or gcd(norm, discriminant) != 1:
        raise ValueError("the primitive norm must be positive and coprime to Delta")
    if 4*radius*radius >= -discriminant:
        raise ValueError("the outer radius must satisfy 4 R^2 < |Delta|")
    factors = smooth_factors(norm, smooth_bound)
    if factors is None:
        return None
    divisor = balanced_divisor_with_gap(factors, radius, inverse_slack)
    if divisor is None:
        return None
    return pair_from_primitive_divisor(form, x, y, divisor, radius)
