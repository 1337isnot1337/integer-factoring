"""Exact all-one-weight profiles on explicit binary subsets of word labels."""

from functools import lru_cache
from itertools import product

from hamming_tensor_power import local_terms
from modular_pfaffian import polynomial_product


@lru_cache(maxsize=4096)
def polynomial_power(value, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = polynomial_product(result, value)
        value = polynomial_product(value, value)
        exponent >>= 1
    return result


@lru_cache(maxsize=4096)
def histogram_profile(histogram):
    """Four-bit pattern counts; fixed low and leading bits stay unchanged."""
    assert len(histogram) == 16 and all(n >= 0 for n in histogram)
    result = 0
    for local, scalar in local_terms():
        value = polynomial_product(scalar, local[0])
        for entry, count in zip(local, histogram):
            if count:
                if not entry:
                    value = 0
                    break
                value = polynomial_product(value, polynomial_power(entry, count))
        result ^= value
    return result


def word_histogram(bit_count, indices):
    assert len(indices) == 4
    assert all(0 <= word < 1 << bit_count for word in indices)
    histogram = [0] * 16
    for bit in range(bit_count):
        pattern = sum(((word >> bit) & 1) << (3 - j) for j, word in enumerate(indices))
        histogram[pattern] += 1
    return tuple(histogram)


def restricted_matrices(bit_count, words):
    assert bit_count >= 1 and len(set(words)) == len(words)
    dimension = len(words)
    result = {}
    for indices in product(range(dimension), repeat=4):
        polynomial = histogram_profile(
            word_histogram(bit_count, tuple(words[i] for i in indices))
        )
        x, y, u, v = indices
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        while polynomial:
            bit = polynomial & -polynomial
            degree = bit.bit_length() - 1
            result.setdefault(degree, [0] * dimension**2)[x * dimension + y] |= 1 << (
                u * dimension + v
            )
            polynomial ^= bit
    return result


def linear_code(generators):
    words = [0]
    for row in generators:
        words += [word ^ row for word in words]
    assert len(set(words)) == len(words)
    return tuple(words)


def affine_histogram(variables, indices):
    """Index bits encode b then a in the affine function a dot t+b."""
    assert variables >= 1 and len(indices) == 4
    assert all(0 <= i < 1 << (variables + 1) for i in indices)
    constant = sum((i & 1) << (3 - j) for j, i in enumerate(indices))
    span = {0}
    for bit in range(variables):
        vector = sum(((i >> (bit + 1)) & 1) << (3 - j) for j, i in enumerate(indices))
        span |= {value ^ vector for value in span}
    histogram = [0] * 16
    multiplicity = (1 << variables) // len(span)
    for value in span:
        histogram[value ^ constant] = multiplicity
    return tuple(histogram)


def affine_profile(variables, indices):
    return histogram_profile(affine_histogram(variables, indices))
