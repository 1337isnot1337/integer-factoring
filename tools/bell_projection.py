"""Exact signed intertwiners of the folded crossing, in characteristic three."""

from collections import Counter
from fractions import Fraction
from math import log2
import random

from code_folded_bell import span
from folded_bell import polynomials
from search_bell_profiles import optimum, pressure


def trim(a):
    a = list(a)
    while a and not a[-1]:
        a.pop()
    return a


def divide(a, b):
    a, b = trim(a), trim(b)
    result = [0] * max(0, len(a) - len(b) + 1)
    while a and len(a) >= len(b):
        j = len(a) - len(b)
        scale = a[-1] * pow(b[-1], -1, 3) % 3
        result[j] = scale
        for i, value in enumerate(b):
            a[i + j] = (a[i + j] - scale * value) % 3
        a = trim(a)
    return trim(result), a


def gcd(a, b):
    while b:
        a, b = b, divide(a, b)[1]
    return [x * pow(a[-1], -1, 3) % 3 for x in a]


def quadratic_transform(values, vectors):
    """Quadratic convolution, via Walsh transform and the exact polar map."""
    q = len(vectors)
    bits = q.bit_length() - 1
    phases = [(-1) ** (v.bit_count() // 2 % 2) for v in vectors]
    d = 1 << (bits // 2)
    gaussian = sum(phases)
    assert abs(gaussian) == d
    result = [x * s % 3 for x, s in zip(values, phases)]
    step = 1
    while step < q:
        for block in range(0, q, 2 * step):
            for j in range(block, block + step):
                a, b = result[j], result[j + step]
                result[j], result[j + step] = (a + b) % 3, (a - b) % 3
        step *= 2
    basis = [vectors[1 << j] for j in range(bits)]
    output = []
    scale = gaussian // d * pow(d, -1, 3)
    for v, phase in zip(vectors, phases):
        index = sum(((v & b).bit_count() % 2) << j for j, b in enumerate(basis))
        output.append(scale * phase * result[index] % 3)
    return tuple(output)


def split(n, discarded_bits, seed, mixes):
    remaining = [(1 << i) | (1 << (n - 1)) for i in range(n - 1)]
    pairs = []
    dot = lambda a, b: (a & b).bit_count() % 2
    while remaining:
        a = remaining.pop(0)
        index = next(i for i, b in enumerate(remaining) if dot(a, b))
        b = remaining.pop(index)
        pairs.extend((a, b))
        remaining = [v ^ (a if dot(v, b) else 0) ^ (b if dot(v, a) else 0)
                     for v in remaining]
    rng = random.Random(seed)
    transvections = []
    for _ in range(mixes):
        a = rng.randrange(1, 1 << n)
        while a.bit_count() % 2 != 0:
            a = rng.randrange(1, 1 << n)
        pairs = [v ^ a if dot(v, a) else v for v in pairs]
        transvections.append(a)
    return (span(tuple(pairs[discarded_bits:])), span(tuple(pairs[:discarded_bits])),
            transvections)


def project(physical_bits, boundary_bits, seed=20261031, mixes=8, boundary="cubic",
            supplied_weights=None):
    assert boundary_bits % 2 == 0 and boundary_bits >= 2
    n = 2 * physical_bits + boundary_bits + 1
    k = (n - 1) // 2
    v, d, transvections = split(n, boundary_bits, seed, mixes)
    if boundary == "supplied":
        weights = tuple(supplied_weights)
    elif boundary == "cubic":
        assert boundary_bits % 6 == 0
        f = []
        for x in range(len(d)):
            phase = 0
            for start in range(0, boundary_bits, 6):
                b = x >> start
                phase ^= ((b >> 0 & 1) & (b >> 1 & 1) & (b >> 2 & 1))
                phase ^= ((b >> 2 & 1) & (b >> 3 & 1) & (b >> 4 & 1))
            f.append(2 if phase else 1)
    elif boundary == "delta":
        f = [1] + [0] * (len(d) - 1)
    else:
        raise ValueError("unknown boundary")
    if boundary != "supplied":
        tf = quadratic_transform(f, d)
        normalization = 2 * sum(f) % 3
        assert normalization
        weights = tuple((a + b) * pow(normalization, -1, 3) % 3 for a, b in zip(f, tf))
    assert sum(weights) % 3 == 1 and quadratic_transform(weights, d) == weights
    local = polynomials(k)
    source = []
    for x in v:
        polynomial = [0] * (k + 1)
        for coefficient, y in zip(weights, d):
            if coefficient:
                w = (x ^ y).bit_count()
                ell = min(w, n - w)
                for i, value in enumerate(local[k - ell]):
                    polynomial[ell + i] = (polynomial[ell + i] + coefficient * value) % 3
        source.append(trim(polynomial))
    common = source[0]
    for row in source[1:]:
        common = gcd(common, row)
    initial_order = next(i for i, value in enumerate(common) if value)
    degree = k - initial_order - (len(common) - 1)
    common = [x * pow(sum(common) % 3, -1, 3) % 3 for x in common]
    reduced = []
    for row in source:
        polynomial, remainder = divide(row, common)
        assert not remainder and sum(polynomial) % 3 == 1
        assert len(polynomial) <= degree + 1
        reduced.append(tuple(polynomial + [0] * (degree + 1 - len(polynomial))))
    coefficients = tuple(tuple(row[i] for row in reduced) for i in range(degree + 1))
    assert tuple(quadratic_transform(row, v) for row in coefficients) == tuple(reversed(coefficients))
    histogram = Counter(next(i for i, value in enumerate(row) if value) for row in reduced)
    low = tuple(histogram[i] for i in range(degree + 1))
    _, approximate_x = optimum(low, tuple(reversed(low)))
    x = Fraction(approximate_x).limit_denominator(256)
    bound = pressure(low, tuple(reversed(low)), x)
    return {"physical_bits": physical_bits, "boundary_bits": boundary_bits,
            "seed": seed, "mixes": mixes, "boundary": boundary,
            "transvections": transvections, "ambient_degree": k,
            "common_divisor": common, "removed_initial_order": initial_order,
            "degree": degree, "valuation_histogram": low,
            "pressure": [bound.numerator, bound.denominator],
            "scale": [x.numerator, x.denominator],
            "rho_pressure": (log2(bound.numerator) - log2(bound.denominator)) / (4 * physical_bits),
            "strict_target": bound ** 577 < 2 ** (4 * physical_bits * 408),
            "coefficients_by_degree": coefficients}
