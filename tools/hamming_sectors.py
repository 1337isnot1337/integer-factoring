"""Crossing profiles on fixed-size subsets of antipodal mode pairs."""

from itertools import combinations
from math import comb

from hamming_pfaffian import coefficient_matrices, labels


def sector_labels(free_width, pair_count, limit=None):
    assert free_width >= 1 and 1 <= pair_count <= 1 << free_width
    pairs = labels(free_width)
    choices = combinations(range(len(pairs)), pair_count)
    result = []
    for selected in choices:
        result.append(tuple(sorted(part for a in selected for part in pairs[a])))
        if limit is not None and len(result) == limit:
            break
    assert limit is None or len(result) == limit
    return tuple(result)


def binary_sector_labels(width, pair_count):
    assert width >= 1 and pair_count >= 1
    free_width = 1
    while comb(1 << free_width, pair_count) < 1 << width:
        free_width += 1
    return free_width, sector_labels(free_width, pair_count, 1 << width)


def sector_matrices(free_width, pair_count, limit=None, weights=None):
    if weights is None:
        weights = (1,)*(free_width+2)
    assert len(weights) == free_width+2 and weights[0] == 1 and min(weights) >= 1
    domain = sector_labels(free_width, pair_count, limit)
    weighted_size = pair_count*(2+sum(weights[1:]))
    return coefficient_matrices(domain, 4*weighted_size, weights)


def complement_scalar(weights):
    """D=TQ+Q+T for a full 2^(len(weights)-1)-mode core."""
    from modular_pfaffian import polynomial_product
    assert weights and weights[0] == 1 and min(weights) >= 1
    t = 1 << (2+sum(weights[1:]))
    q = 5  # (1+z)^2, the low-bit contribution.
    for weight in weights[1:]:
        q = polynomial_product(q, 1 | (1 << weight))
    return polynomial_product(t, q) ^ q ^ t
