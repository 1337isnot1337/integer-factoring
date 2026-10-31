"""Find a primitive representation of one supplied factored candidate norm.

The exact search enumerates prime-ideal orientations by meet-in-the-middle.
It does not choose the candidate norm.  See research/factored-norm-search.md.
"""

from bisect import bisect_left
from math import gcd, prod

from bounded_smooth_norm import smooth_factors
from quadratic_form_arithmetic import (
    compose_primitive, form_discriminant, inverse_primitive, reduce_with_transform,
)


def lift_prime_root(discriminant, prime, exponent, root):
    """Lift a valid root modulo 4p to the oriented norm-p^e form root."""
    if exponent < 1 or prime < 2 or gcd(prime, discriminant) != 1:
        raise ValueError("a positive unramified prime power is required")
    if (root*root-discriminant) % (4*prime):
        raise ValueError("the supplied root is invalid modulo 4p")
    if prime == 2:
        value = root % 4
        for level in range(1, exponent):
            scale = 1 << (level+2)
            assert (discriminant-value*value) % scale == 0
            bit = ((discriminant-value*value)//scale) % 2
            value += bit*(1 << (level+1))
        modulus = 1 << (exponent+1)
        value %= modulus
    else:
        value, scale = root % prime, prime
        inverse_derivative = pow(2*value, -1, prime)
        for _ in range(1, exponent):
            assert (discriminant-value*value) % scale == 0
            digit = ((discriminant-value*value)//scale)*inverse_derivative % prime
            value += digit*scale
            scale *= prime
        # Parity supplies the additional factor four in the form congruence.
        if value % 2 != discriminant % 2:
            value += scale
        modulus = 2*scale
    norm = prime**exponent
    assert modulus == 2*norm and 0 <= value < modulus
    assert (value*value-discriminant) % (4*norm) == 0
    return value


def combine_roots(discriminant, prime_powers, lifted_roots, positive_mask):
    """CRT for a root modulo 2n satisfying the stronger modulus-4n test."""
    residues = []
    contains_two = False
    for index, ((prime, exponent), root) in enumerate(zip(prime_powers, lifted_roots)):
        modulus = prime**exponent
        if prime == 2:
            modulus *= 2
            contains_two = True
        residue = root if positive_mask & (1 << index) else -root
        residues.append((residue % modulus, modulus))
    if not contains_two:
        residues.append((discriminant % 2, 2))
    value, modulus = 0, 1
    for residue, next_modulus in residues:
        assert gcd(modulus, next_modulus) == 1
        digit = (residue-value)*pow(modulus, -1, next_modulus) % next_modulus
        value += digit*modulus
        modulus *= next_modulus
    integer = prod(prime**exponent for prime, exponent in prime_powers)
    assert modulus == 2*integer
    assert (value*value-discriminant) % (4*integer) == 0
    return value


def represent_with_prime_powers(form, prime_powers, prime_roots=None):
    """Decide primitive representation of a supplied coprime factored norm.

    The prime_powers input is assumed to be a complete prime factorization;
    represent_smooth_norm constructs it by trial division.  Optional roots
    avoid the finite root search.  A positive answer is checked directly.
    """
    discriminant = form_discriminant(form)
    target, target_transform = reduce_with_transform(form)
    prime_powers = tuple(sorted(prime_powers))
    if (any(prime < 2 or exponent < 1 for prime, exponent in prime_powers)
            or len({prime for prime, _ in prime_powers}) != len(prime_powers)):
        raise ValueError("distinct positive prime powers are required")
    integer = prod(prime**exponent for prime, exponent in prime_powers)
    if gcd(integer, discriminant) != 1:
        raise ValueError("the candidate norm must be coprime to the discriminant")
    metadata = {"norm": integer, "omega": len(prime_powers),
                "norm_bits": integer.bit_length(), "prime_powers": prime_powers}
    parity = discriminant % 2
    identity = (1, parity, (parity-discriminant)//4)
    root_forms = []
    lifted_roots = []
    for prime, exponent in prime_powers:
        if prime_roots is None:
            root = next((value for value in range(prime+1)
                         if (value*value-discriminant) % (4*prime) == 0), None)
        else:
            if prime not in prime_roots:
                raise ValueError("a supplied prime root is missing")
            root = prime_roots[prime]
        if root is None:
            return {"status": "no_primitive_representation", "reason": "inert_prime",
                    "inert_prime": prime, **metadata}
        lifted = lift_prime_root(discriminant, prime, exponent, root)
        prime_norm = prime**exponent
        root_form, _ = reduce_with_transform((prime_norm, lifted,
            (lifted*lifted-discriminant)//(4*prime_norm)))
        root_forms.append(root_form)
        lifted_roots.append(lifted)
    composition_calls = 0

    def multiply(left, right):
        nonlocal composition_calls
        composition_calls += 1
        return compose_primitive(left, right)

    def half_states(start, stop):
        states = [(identity, 0)]
        peak = 1
        for index in range(start, stop):
            positive = root_forms[index]
            negative = inverse_primitive(positive)
            candidates = []
            for value, mask in states:
                candidates.append((multiply(value, negative), mask))
                candidates.append((multiply(value, positive), mask | (1 << index)))
            candidates.sort()
            next_states = []
            for value, mask in candidates:
                if not next_states or next_states[-1][0] != value:
                    next_states.append((value, mask))
            states = next_states
            peak = max(peak, len(states))
        return states, peak

    middle = len(prime_powers)//2
    left, left_peak = half_states(0, middle)
    right, right_peak = half_states(middle, len(prime_powers))
    right_values = [value for value, _ in right]
    found = None
    for value, mask in left:
        needed = multiply(inverse_primitive(value), target)
        position = bisect_left(right_values, needed)
        if position < len(right) and right[position][0] == needed:
            found = mask | right[position][1]
            break
    metadata.update({"half_state_counts": (len(left), len(right)),
                     "peak_half_states": (left_peak, right_peak),
                     "composition_calls": composition_calls})
    if found is None:
        return {"status": "no_primitive_representation", "reason": "orientation_classes",
                **metadata}
    root = combine_roots(discriminant, prime_powers, lifted_roots, found)
    auxiliary = (integer, root, (root*root-discriminant)//(4*integer))
    reduced, transform = reduce_with_transform(auxiliary)
    assert reduced == target
    # reduced = auxiliary o transform.  The inverse image of (1,0)
    # under transform is a primitive representation in the reduced target.
    reduced_x, reduced_y = transform[1][1], -transform[1][0]
    x = target_transform[0][0]*reduced_x + target_transform[0][1]*reduced_y
    y = target_transform[1][0]*reduced_x + target_transform[1][1]*reduced_y
    a, b, c = form
    assert gcd(x, y) == 1 and a*x*x+b*x*y+c*y*y == integer
    return {"status": "represented", "primitive_vector": (x, y),
            "positive_orientation_mask": found, "combined_root": root,
            "auxiliary_form": auxiliary, **metadata}


def represent_smooth_norm(form, integer, smooth_bound):
    """Recognize smoothness, then exactly test this one candidate norm."""
    if integer < 1 or smooth_bound < 2:
        raise ValueError("a positive norm and a smoothness bound at least two are required")
    factors = smooth_factors(integer, smooth_bound)
    if factors is None:
        return {"status": "not_smooth", "norm": integer, "smooth_bound": smooth_bound}
    prime_powers = tuple((prime, factors.count(prime)) for prime in sorted(set(factors)))
    return represent_with_prime_powers(form, prime_powers)
