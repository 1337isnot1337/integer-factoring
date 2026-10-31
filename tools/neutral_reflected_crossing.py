"""Exact reflected neutral-fermion crossings in characteristic two.

Finite Laurent polynomials use (lowest exponent, packed F2 coefficients).
The kernel and fixed-sector endpoint are proved in
research/neutral-reflected-crossings.md.  Source comparison uses the
Appendix B W2 quotient, with parameter descent proved separately.
"""

from functools import lru_cache
from itertools import combinations
from math import gcd


def normalize(shift, bits):
    if not bits:
        return (0, 0)
    valuation = (bits & -bits).bit_length() - 1
    return (shift + valuation, bits >> valuation)


def laurent_add(left, right):
    if not left[1]:
        return right
    if not right[1]:
        return left
    shift = min(left[0], right[0])
    return normalize(shift, (left[1] << (left[0] - shift))
                     ^ (right[1] << (right[0] - shift)))


@lru_cache(None)
def polynomial_multiply(left, right):
    result = 0
    while right:
        bit = right & -right
        result ^= left << (bit.bit_length() - 1)
        right ^= bit
    return result


def laurent_multiply(left, right):
    return normalize(left[0] + right[0],
                     polynomial_multiply(left[1], right[1]))


@lru_cache(None)
def binomial_polynomial(exponent):
    """Packed (1+z)^exponent over F2, by Frobenius."""
    if exponent < 0:
        raise ValueError("negative binomial exponent")
    result, bit = 1, 0
    while exponent:
        if exponent & 1:
            result ^= result << (1 << bit)
        exponent >>= 1
        bit += 1
    return result


def sector_labels(period, length, weight):
    if period < 2 or period & (period - 1):
        raise ValueError("period must be a positive power of two, at least two")
    if not 0 <= length <= period // 2:
        raise ValueError("invalid occupation length")
    return tuple(parts for parts in combinations(range(1, period, 2), length)
                 if sum(parts) == weight)


def kernel_entry(period, first_group, first_index, second_group, second_index):
    """One entry of the symmetric zero-diagonal Wick matrix.

    Groups 0,1,2,3 are the reflected parameter, second parameter,
    input state, and output coefficient. All indices are positive odd
    integers strictly below period.
    """
    if first_group > second_group:
        return kernel_entry(period, second_group, second_index,
                            first_group, first_index)
    a, b = first_group, second_group
    k, t = first_index, second_index
    if a == 0 and b == 1:
        return (0, int(t & ~k == 0))
    if a == 0 and b == 2 and t & ~k == 0:
        return (k - period - t, binomial_polynomial(period - k))
    if a == 0 and b == 3 and (period - k - 1) & ~(t - 1) == 0:
        return (k - period + t, binomial_polynomial(period - k))
    if a == 1 and b == 2 and (k - 1) & t == 0:
        return (0, binomial_polynomial(k))
    if a == 1 and b == 3 and (k - 1) & ~(t - 1) == 0:
        return (0, binomial_polynomial(k))
    if a == 2 and b == 3 and k == t:
        return (0, 1)
    return (0, 0)


def crossing_entry(period, x, y, v, u):
    """The exact Laurent entry in row/column order (x,y)|(v,u)."""
    length, weight = len(x), sum(x)
    if any(len(label) != length or sum(label) != weight
           for label in (x, y, v, u)):
        raise ValueError("all labels must have the same length and weight")
    if any(tuple(sorted(set(label))) != tuple(label)
           or any(k <= 0 or k >= period or not k % 2 for k in label)
           for label in (x, y, v, u)):
        raise ValueError("labels must be increasing distinct odd indices below period")
    fields = (tuple((0, k) for k in x) + tuple((1, k) for k in v)
              + tuple((2, k) for k in y) + tuple((3, k) for k in u))
    size = len(fields)
    entries = {(i, j): kernel_entry(period, *fields[i], *fields[j])
               for i in range(size) for j in range(i + 1, size)}

    @lru_cache(None)
    def pfaffian(mask):
        if not mask:
            return (0, 1)
        first = (mask & -mask).bit_length() - 1
        rest = mask ^ (1 << first)
        remaining, result = rest, (0, 0)
        while remaining:
            second = (remaining & -remaining).bit_length() - 1
            remaining ^= 1 << second
            entry = entries[first, second]
            if entry[1]:
                result = laurent_add(result, laurent_multiply(
                    entry, pfaffian(rest ^ (1 << second))))
        return result

    return pfaffian((1 << size) - 1)


def tensor_coefficients(period, labels):
    """Return literal Laurent coefficient matrices as packed F2 rows."""
    labels = tuple(labels)
    if not labels:
        raise ValueError("empty label set")
    dimension = len(labels)
    matrices = {}
    for xi, x in enumerate(labels):
        for yi, y in enumerate(labels):
            for vi, v in enumerate(labels):
                for ui, u in enumerate(labels):
                    shift, bits = crossing_entry(period, x, y, v, u)
                    while bits:
                        bit = bits & -bits
                        exponent = shift + bit.bit_length() - 1
                        bits ^= bit
                        if exponent not in matrices:
                            matrices[exponent] = [0] * dimension ** 2
                        matrices[exponent][xi * dimension + yi] ^= 1 << (vi * dimension + ui)
    return dict(sorted(matrices.items()))


def dense_polynomial_seed(period, labels, compress=True):
    """Return a strict polynomial seed for the native Smith oracles.

    Common monomial translation is always exact. Optional exponent
    compression is checked from every nonzero coefficient; it is not
    extrapolated from sampled entries or an asserted all-size pattern.
    Dense order is x + d*y + d**2*u + d**3*v.
    """
    labels = tuple(labels)
    dimension = len(labels)
    matrices = tensor_coefficients(period, labels)
    least = min(matrices)
    step = 0
    if compress:
        for exponent in matrices:
            step = gcd(step, exponent - least)
    step = max(step, 1)
    values = [0] * dimension ** 4
    for exponent, rows in matrices.items():
        assert (exponent - least) % step == 0
        bit = 1 << ((exponent - least) // step)
        for xi in range(dimension):
            for yi in range(dimension):
                row = rows[xi * dimension + yi]
                for vi in range(dimension):
                    for ui in range(dimension):
                        if row >> (vi * dimension + ui) & 1:
                            values[xi + dimension * yi + dimension ** 2 * ui
                                   + dimension ** 3 * vi] ^= bit
    return values, {"common_laurent_shift": least, "verified_exponent_step": step,
                    "original_nonzero_exponents": list(matrices),
                    "polynomial_degree": max(value.bit_length() - 1 for value in values)}
