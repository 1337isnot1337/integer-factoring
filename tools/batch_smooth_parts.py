"""Exact batched smooth-part tests by product and remainder trees.

This is an application of D. J. Bernstein's 2004 smooth-parts algorithm;
see research/batch-smooth-filtering.md for attribution and bit costs.

No coprimality promise is needed for the candidate moduli.  The caller supplies
the product of the allowed primes and an exponent bounding every candidate's
prime valuations.  This module does not find or certify the allowed primes.
"""

from math import gcd


def smoothness_exponent(bit_bound):
    """A power-of-two exponent at least the positive bit bound."""
    if not isinstance(bit_bound, int) or bit_bound < 1:
        raise ValueError('The bit bound must be a positive integer.')
    return 1 << (bit_bound - 1).bit_length()


def _positive_integers(values):
    values = tuple(values)
    if any(not isinstance(value, int) or value < 1 for value in values):
        raise ValueError('Candidate moduli must be positive integers.')
    return values


def product_tree(moduli):
    """Return all balanced product levels, beginning with the leaves."""
    moduli = _positive_integers(moduli)
    if not moduli:
        return []
    tree = [list(moduli)]
    while len(tree[-1]) > 1:
        level = tree[-1]
        tree.append([level[i] * (level[i + 1] if i + 1 < len(level) else 1)
                     for i in range(0, len(level), 2)])
    return tree


def simultaneous_remainders(value, moduli):
    """Return value modulo every modulus, including repeats/shared factors."""
    if not isinstance(value, int) or value < 0:
        raise ValueError('The common dividend must be a nonnegative integer.')
    tree = product_tree(moduli)
    if not tree:
        return []
    remainders = [value % tree[-1][0]]
    for level in reversed(tree[:-1]):
        remainders = [remainders[i // 2] % modulus
                      for i, modulus in enumerate(level)]
    return remainders


def _smooth_part_from_residue(value, residue, exponent):
    """Use two bounded gcd saturations, then the worst-case power bound."""
    smooth = gcd(value, residue)
    remaining = value // smooth
    if smooth == 1 or remaining == 1:
        return smooth
    for _ in range(2):
        common = gcd(smooth, remaining)
        if common == 1:
            return smooth
        smooth *= common
        remaining //= common
        if remaining == 1:
            return smooth
    # Every supported prime now has valuation at least min(4,v_p(value)).
    # A quarter of the original safe exponent therefore suffices.
    return gcd(value, pow(smooth, (exponent + 3) // 4, value))


def smooth_parts(moduli, prime_product, *, exponent=None):
    """Return each modulus's largest divisor supported on prime_product.

    Every positive candidate a has v_p(a) <= bit_length(a).  Thus an exponent
    at least the largest candidate bit length is always safe.  A supplied
    exponent must satisfy this conservative bound; no heuristic is used.
    """
    moduli = _positive_integers(moduli)
    if not isinstance(prime_product, int) or prime_product < 1:
        raise ValueError('The prime product must be a positive integer.')
    bound = max((value.bit_length() for value in moduli), default=1)
    if exponent is None:
        exponent = smoothness_exponent(bound)
    if not isinstance(exponent, int) or exponent < bound:
        raise ValueError('The exponent must bound all candidate bit lengths.')
    residues = simultaneous_remainders(prime_product, moduli)
    return [_smooth_part_from_residue(value, residue, exponent)
            for value, residue in zip(moduli, residues)]


def smooth_flags(moduli, prime_product, *, exponent=None):
    """Return exact base-smoothness flags without factoring any candidate."""
    moduli = _positive_integers(moduli)
    if not isinstance(prime_product, int) or prime_product < 1:
        raise ValueError('The prime product must be a positive integer.')
    bound = max((value.bit_length() for value in moduli), default=1)
    if exponent is None:
        exponent = smoothness_exponent(bound)
    if not isinstance(exponent, int) or exponent < bound:
        raise ValueError('The exponent must bound all candidate bit lengths.')
    residues = simultaneous_remainders(prime_product, moduli)
    return [_smooth_part_from_residue(value, residue, exponent) == value
            for value, residue in zip(moduli, residues)]


def smooth_sequence_batches(initial, successor, count, prime_product, batch_size,
                            *, norm=lambda value: value, exponent=None):
    """Yield offset, states, and smooth parts for a fixed-product sequence.

    Blocks grow 1, 2, 4, ... up to batch_size, then retain that cap.  A caller
    stopping at its first accepted state discards at most batch_size-1 states.
    The supplied count must already include any remaining global query budget.
    """
    if not isinstance(count, int) or count < 0:
        raise ValueError('The sequence count must be a nonnegative integer.')
    if not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError('The batch size must be a positive integer.')
    if not isinstance(prime_product, int) or prime_product < 1:
        raise ValueError('The prime product must be a positive integer.')
    offset, block_size, state = 0, 1, initial
    while offset < count:
        take = min(block_size, count - offset)
        states = []
        for index in range(take):
            states.append(state)
            if index + 1 < take:
                state = successor(state)
        values = [norm(value) for value in states]
        yield offset, states, smooth_parts(values, prime_product, exponent=exponent)
        offset += take
        block_size = min(2 * block_size, batch_size)
        if offset < count:
            state = successor(state)
