#!/usr/bin/env python3
"""Scale-one F2 crossings on partitions into even-sized blocks."""

from collections import Counter
from functools import lru_cache
from itertools import combinations

from fermionic_piece_polynomial import polynomial_product


def partition_labels(vertices, block_size):
    if block_size < 2 or block_size % 2 or vertices % block_size:
        raise ValueError("The block size must be even and divide the vertex count.")

    def generate(remaining):
        if not remaining:
            yield ()
            return
        first = remaining & -remaining
        rest = [i for i in range(vertices) if (remaining ^ first) >> i & 1]
        for chosen in combinations(rest, block_size - 1):
            block = first | sum(1 << i for i in chosen)
            for tail in generate(remaining ^ block):
                yield tuple(sorted((block,) + tail))

    return tuple(generate((1 << vertices) - 1))


def validate_labels(labels, vertices, block_size):
    full = (1 << vertices) - 1
    for label in labels:
        used = 0
        for block in label:
            if block <= 0 or block.bit_count() != block_size or block & used or block & ~full:
                raise ValueError("Each label must be a uniform set partition.")
            used |= block
        if used != full:
            raise ValueError("Each label must cover every vertex.")


def kernel(first, second):
    """I plus the physical-vertex incidence Gram matrix, over F2."""
    return int(first == second) ^ ((first & second).bit_count() & 1)


def particle_matrix(operator_a, operator_b, output, input_, vertices, block_size):
    """Nonzero entries as pairs (z-block mask, t-block mask)."""
    validate_labels((operator_a, operator_b, output, input_), vertices, block_size)
    particles = (
        [("U", e) for e in output] + [("B", e) for e in operator_b]
        + [("A", e) for e in operator_a] + [("Y", e) for e in input_]
    )
    count = len(particles)
    entries = [[None] * count for _ in range(count)]
    for i, (left, e) in enumerate(particles):
        for j in range(i + 1, count):
            right, f = particles[j]
            pair = left, right
            if pair == ("U", "B") and e == f:
                entries[i][j] = (0, f)
            elif pair == ("U", "A") and e == f:
                entries[i][j] = (e, f)
            elif pair == ("U", "Y") and e == f:
                entries[i][j] = (f, 0)
            elif pair == ("B", "A") and kernel(e, f):
                entries[i][j] = (f, 0)
            elif pair == ("B", "Y") and kernel(e, f):
                entries[i][j] = (f, e)
            elif pair == ("A", "Y") and kernel(e, f):
                entries[i][j] = (0, e)
    return entries


def coefficient_terms(operator_a, operator_b, output, input_, vertices, block_size):
    """Collect products L^z(1+L)^t prod X_i^a_i(1+X_i)^b_i.

    This deliberately simple finite-check routine enumerates nonzero
    Pfaffian matchings; the full-polynomial routine below uses subset DP.
    """
    entries = particle_matrix(operator_a, operator_b, output, input_, vertices, block_size)
    count = len(entries)
    terms = Counter()

    def recurse(remaining, z_masks, t_masks):
        if not remaining:
            z_counts = tuple(sum(mask >> i & 1 for mask in z_masks) for i in range(vertices))
            t_counts = tuple(sum(mask >> i & 1 for mask in t_masks) for i in range(vertices))
            term = len(z_masks), len(t_masks), z_counts, t_counts
            terms[term] ^= 1
            return
        first_bit = remaining & -remaining
        first = first_bit.bit_length() - 1
        rest = remaining ^ first_bit
        candidates = rest
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            other = bit.bit_length() - 1
            entry = entries[first][other]
            if entry is not None:
                z, t = entry
                recurse(rest ^ bit, z_masks + ((z,) if z else ()),
                        t_masks + ((t,) if t else ()))

    recurse((1 << count) - 1, (), ())
    return tuple(term for term, parity in terms.items() if parity)


def odd_binomial(total, chosen):
    return 0 <= chosen <= total and not chosen & ~total


def coefficient(terms, l_degree, degrees):
    result = 0
    for z, t, z_counts, t_counts in terms:
        if not odd_binomial(t, l_degree - z):
            continue
        if all(odd_binomial(b, degree - a)
               for degree, a, b in zip(degrees, z_counts, t_counts)):
            result ^= 1
    return result


def univariate(terms, physical_weights, l_weight=1):
    result = 0
    for z, t, z_counts, t_counts in terms:
        polynomial = 1 << (z * l_weight + sum(
            weight * exponent for weight, exponent in zip(physical_weights, z_counts)))
        for _ in range(t):
            polynomial ^= polynomial << l_weight
        for weight, exponent in zip(physical_weights, t_counts):
            for _ in range(exponent):
                polynomial ^= polynomial << weight
        result ^= polynomial
    return result


def lifted_polynomial(operator_a, operator_b, output, input_, vertices, block_size):
    """Exact multigraded polynomial, with radix4q+1 followed by radix5."""
    entries = particle_matrix(operator_a, operator_b, output, input_, vertices, block_size)
    count = len(entries)
    radix = count + 1
    places = tuple(radix * 5**i for i in range(vertices))

    @lru_cache(None)
    def monomial(mask):
        return 1 << sum(places[i] for i in range(vertices) if mask >> i & 1)

    @lru_cache(None)
    def binomial(mask):
        result = 1
        for i in range(vertices):
            if mask >> i & 1:
                result ^= result << places[i]
        return result

    polynomials = [[0] * count for _ in range(count)]
    for i in range(count):
        for j in range(i + 1, count):
            if entries[i][j] is not None:
                z, t = entries[i][j]
                value = monomial(z) << int(bool(z))
                if t:
                    value = polynomial_product(value, polynomial_product(3, binomial(t)))
                polynomials[i][j] = value

    @lru_cache(None)
    def recurse(remaining):
        if not remaining:
            return 1
        first_bit = remaining & -remaining
        first = first_bit.bit_length() - 1
        rest = remaining ^ first_bit
        result = 0
        candidates = rest
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            other = bit.bit_length() - 1
            if polynomials[first][other]:
                result ^= polynomial_product(polynomials[first][other], recurse(rest ^ bit))
        return result

    return recurse((1 << count) - 1)
