"""Fixed smooth-square certificates from supplied bounded class-group words.

This module does not find those words.  Its forward constructor accepts one,
and its reverse converter accepts a representation of the fixed square.
See research/universal-smooth-square.md for the exact conditional gate.
"""

from math import gcd, prod

from factored_norm_search import combine_roots, lift_prime_root
from quadratic_form_arithmetic import form_discriminant, reduce_with_transform
from smooth_norm_relation import norm_to_prime_word


def universal_square_parameters(primes, radius, index_bound):
    """Choose the cap from density data before any relation target is chosen.

    The caller certifies |G| <= index_bound * |S_radius|, with S_radius
    containing the identity and the prime forms at the supplied primes.
    The function checks the elementary numerical input conditions only.
    """
    primes = tuple(primes)
    if (radius < 2 or index_bound < 1 or len(set(primes)) != len(primes)
            or any(prime < 2 or prime > radius for prime in primes)):
        raise ValueError("distinct primes at most radius and a positive density bound are required")
    leading_word_bound = radius.bit_length()-1
    cap = 2*(2*index_bound-1)*leading_word_bound
    root_norm = prod(prime**cap for prime in primes)
    target_exponent = 2*(cap*len(primes)+1)
    return {"primes": primes, "radius": radius, "index_bound": index_bound,
            "leading_word_bound": leading_word_bound, "exponent_cap": cap,
            "root_norm": root_norm, "fixed_norm": root_norm*root_norm,
            "target_exponent": target_exponent,
            "primitive_word_bound": 2*cap*len(primes)}


def representation_from_square_word(form, primes, roots, word, exponent_cap):
    """Construct Q(x,y)=M^2 from h=product(g_p^(2*word[p])).

    The supplied word must have each absolute exponent at most exponent_cap.
    This is a certificate constructor, not a bounded-word search algorithm.
    The target form is checked by reduction and the final equality directly.
    """
    discriminant = form_discriminant(form)
    primes, word = tuple(primes), tuple(word)
    if (exponent_cap < 1 or len(primes) != len(word)
            or len(set(primes)) != len(primes)
            or any(abs(value) > exponent_cap for value in word)
            or any(prime < 2 or gcd(prime, discriminant) != 1 for prime in primes)):
        raise ValueError("a bounded word on distinct unramified primes is required")
    prime_powers, lifted_roots = [], []
    positive_mask = 0
    for prime, exponent in zip(primes, word):
        if prime not in roots or (roots[prime]**2-discriminant) % (4*prime):
            raise ValueError("a valid splitting root is required for each prime")
        if not exponent:
            continue
        position = len(prime_powers)
        prime_powers.append((prime, 2*abs(exponent)))
        lifted_roots.append(lift_prime_root(discriminant, prime,
                                           2*abs(exponent), roots[prime]))
        if exponent > 0:
            positive_mask |= 1 << position
    primitive_norm = prod(prime**exponent for prime, exponent in prime_powers)
    middle = combine_roots(discriminant, prime_powers, lifted_roots, positive_mask)
    auxiliary = (primitive_norm, middle,
                 (middle*middle-discriminant)//(4*primitive_norm))
    reduced_auxiliary, transform = reduce_with_transform(auxiliary)
    reduced_target, target_transform = reduce_with_transform(form)
    if reduced_auxiliary != reduced_target:
        raise ValueError("the supplied squared word does not equal the target class")
    reduced_x, reduced_y = transform[1][1], -transform[1][0]
    primitive_x = target_transform[0][0]*reduced_x + target_transform[0][1]*reduced_y
    primitive_y = target_transform[1][0]*reduced_x + target_transform[1][1]*reduced_y
    content = prod(prime**(exponent_cap-abs(exponent))
                   for prime, exponent in zip(primes, word))
    x, y = content*primitive_x, content*primitive_y
    root_norm = prod(prime**exponent_cap for prime in primes)
    a, b, c = form
    assert gcd(primitive_x, primitive_y) == 1
    assert a*primitive_x**2+b*primitive_x*primitive_y+c*primitive_y**2 == primitive_norm
    assert gcd(x, y) == content and a*x*x+b*x*y+c*y*y == root_norm*root_norm
    return {"vector": (x, y), "primitive_vector": (primitive_x, primitive_y),
            "content": content, "primitive_norm": primitive_norm,
            "root_norm": root_norm, "fixed_norm": root_norm*root_norm,
            "supplied_half_word": word, "positive_orientation_mask": positive_mask}


def word_from_square_representation(form, x, y, primes, exponent_cap):
    """Normalize any exact fixed-square witness and recover its bounded word."""
    primes = tuple(primes)
    if (exponent_cap < 1 or len(set(primes)) != len(primes)
            or any(prime < 2 for prime in primes)):
        raise ValueError("distinct primes and a positive cap are required")
    root_norm = prod(prime**exponent_cap for prime in primes)
    a, b, c = form
    if a*x*x+b*x*y+c*y*y != root_norm*root_norm:
        raise ValueError("the coordinates do not represent the prescribed fixed square")
    content = gcd(x, y)
    assert content > 0 and root_norm % content == 0
    result = norm_to_prime_word(form, x, y, max(primes, default=2))
    assert result is not None and result["norm"] == (root_norm//content)**2
    positions = {prime: index for index, prime in enumerate(primes)}
    vector = [0]*len(primes)
    for entry in result["word"]:
        assert entry["prime"] in positions
        assert entry["exponent"] % 2 == 0 and abs(entry["exponent"]) <= 2*exponent_cap
        vector[positions[entry["prime"]]] = entry["exponent"]
    assert sum(abs(value) for value in vector) <= 2*exponent_cap*len(primes)
    return {"word": vector, "half_word": [value//2 for value in vector],
            "word_length": sum(abs(value) for value in vector),
            "primitive_norm": result["norm"], "content": content,
            "primitive_vector": result["primitive_vector"]}


def dominant_square_relation(generator_index, primes, half_word, exponent_cap):
    """Build a dominant row after its fixed-square group identity is checked.

    The caller certifies product(g_p^(2*half_word[p]))=g_i^(2K), where
    K=exponent_cap*len(primes)+1.  The target is fixed after the cap is fixed.
    """
    size = len(primes)
    if (not 0 <= generator_index < size or len(half_word) != size
            or exponent_cap < 1 or any(abs(value) > exponent_cap for value in half_word)):
        raise ValueError("a valid index and a bounded half-word are required")
    half_diagonal = exponent_cap*size+1
    row = [-2*value for value in half_word]
    row[generator_index] += 2*half_diagonal
    slack = row[generator_index]-sum(abs(value) for index, value in enumerate(row)
                                     if index != generator_index)
    assert slack >= 2
    return row
