"""Enumerate the complete smaller-norm ellipse of a prime-power target.

The near-height running-time bound requires a rejected smooth leading norm.
It is not a theorem that a smaller smooth norm exists. See the research note.
"""
from math import gcd, isqrt

from factored_norm_search import lift_prime_root
from quadratic_form_arithmetic import form_discriminant, reduce_with_transform


def prime_power_target(discriminant, prime, exponent, root):
    """Return the reduced oriented p-power form and its primitive norm witness."""
    norm = prime ** exponent
    middle = lift_prime_root(discriminant, prime, exponent, root)
    form, transform = reduce_with_transform(
        (norm, middle, (middle * middle - discriminant) // (4 * norm)))
    # reduced(z)=original(U*z); U^-1*(1,0) is the known norm witness.
    vector = (transform[1][1], -transform[1][0])
    x, y = vector
    a, b, c = form
    assert gcd(x, y) == 1 and a*x*x + b*x*y + c*y*y == norm
    return form, vector


def smaller_norm_vectors(form, cap):
    """Yield every 0<Q(x,y)<cap once modulo simultaneous negation.

    Include imprimitive vectors: the caller can reject them explicitly.
    Integer square roots and interval endpoints make this an exact exhaustion.
    """
    a, b, _ = form
    discriminant = form_discriminant(form)
    if cap < 1:
        raise ValueError("The strict norm cap must be positive.")
    for y in range(isqrt(4 * a * (cap - 1) // -discriminant) + 1):
        radius = isqrt(4 * a * (cap - 1) + discriminant * y*y)
        low = -((radius + b*y) // (2*a))
        high = (radius - b*y) // (2*a)
        for x in range(low, high + 1):
            if y == 0 and x <= 0:
                continue
            norm = ((2*a*x + b*y)**2 - discriminant*y*y) // (4*a)
            assert 0 < norm < cap
            yield x, y, norm
