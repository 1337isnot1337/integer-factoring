"""An explicit folded-cube Bell crossing over every odd characteristic."""

from fractions import Fraction
from math import comb

def polynomials(k: int) -> tuple[tuple[int, ...], ...]:
    f = [(1,), (2, -1)]
    for m in range(2, k + 1):
        row = [0] * (m + 1)
        for i, value in enumerate(f[-1]):
            row[i] += 2 * value
        for i, value in enumerate(f[-2]):
            row[i + 2] -= value
        f.append(tuple(row))
    return tuple(f[:k + 1])


def even_weight(v: int) -> int:
    w = v.bit_count()
    return w + w % 2


def quadratic(v: int) -> int:
    return even_weight(v) // 2 % 2


def polar(v: int, w: int) -> int:
    return quadratic(v ^ w) ^ quadratic(v) ^ quadratic(w)


def coordinates(k: int) -> tuple[tuple[int, ...], int]:
    """Isometry to hyperbolic coordinates, allowing a column translation."""
    remaining = [1 << i for i in range(2 * k)]
    pairs = []
    while remaining:
        e = remaining.pop(0)
        index = next(i for i, f in enumerate(remaining) if polar(e, f))
        f = remaining.pop(index)
        pairs.append((e, f))
        remaining = [r ^ (e if polar(r, f) else 0) ^ (f if polar(r, e) else 0)
                     for r in remaining]
    ga = sum(quadratic(f) << i for i, (e, f) in enumerate(pairs))
    gb = sum(quadratic(e) << i for i, (e, f) in enumerate(pairs))
    d = 1 << k
    image = []
    for a in range(d):
        for b in range(d):
            v = 0
            for i, (e, f) in enumerate(pairs):
                if a >> i & 1:
                    v ^= e
                if b >> i & 1:
                    v ^= f
            image.append(v)
    return tuple(image), ga * d + gb


def coefficient_vectors(k: int, prime: int) -> tuple[tuple[int, ...], ...]:
    f = polynomials(k)
    image, shift = coordinates(k)
    n, q = 2 * k + 1, 1 << (2 * k)
    result = [[0] * q for _ in range(k + 1)]
    for j in range(q):
        w = even_weight(image[j ^ shift])
        level = min(w, n - w)
        for i, value in enumerate(f[k - level]):
            result[level + i][j] = value % prime
    return tuple(map(tuple, result))


def radial_transform(k: int, output_level: int, input_level: int) -> Fraction:
    n = 2 * k + 1
    a = input_level if input_level % 2 == 0 else n - input_level
    b = output_level if output_level % 2 == 0 else n - output_level
    character_sum = sum((-1) ** h * comb(b, h) * comb(n - b, a - h)
                        for h in range(max(0, a - n + b), min(a, b) + 1))
    return Fraction((-1) ** (k * (k + 1) // 2 + (a + b) // 2) * character_sum,
                    1 << k)


def radial_identity(k: int) -> int:
    f = polynomials(k)
    for level in range(k + 1):
        output = [Fraction(0)] * (k + 1)
        for ell in range(k + 1):
            scalar = radial_transform(k, level, ell)
            for j, value in enumerate(f[k - ell]):
                output[j + ell] += scalar * value
        assert output == list(reversed(f[k - level])) + [0] * level
    return k + 1


def rank_profiles(k: int, prime: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    f = polynomials(k)
    row = [0] * (k + 1)
    for level in range(k + 1):
        multiplicity = comb(2 * k + 1, level)
        for i, value in enumerate(f[k - level]):
            if value % prime:
                row[i + level] += multiplicity
    return tuple(row), tuple(reversed(row))
