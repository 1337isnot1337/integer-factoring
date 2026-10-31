"""Exact Pauli/Bell-diagonal crossings in odd prime characteristic.

The physical crossing is E(x,y,u,v)=[x=v][y=u]. Both normalized
cuts are diagonal, with their eigenvalue vectors related by transform().
"""

from collections import Counter
from functools import lru_cache
from itertools import product
from math import comb, factorial


def sign(value: int) -> int:
    return -1 if value.bit_count() % 2 else 1


def direct_transform(values: tuple[int, ...], bits: int, prime: int) -> tuple[int, ...]:
    d = 1 << bits
    if prime == 2 or len(values) != d * d:
        raise ValueError("requires odd characteristic and d squared entries")
    inverse = pow(d, -1, prime)
    return tuple(inverse * sum(sign((a ^ c) & (b ^ e)) * values[a * d + b]
                               for a in range(d) for b in range(d)) % prime
                 for c in range(d) for e in range(d))


def transform(values: tuple[int, ...], bits: int, prime: int) -> tuple[int, ...]:
    """The same signed transform in O(d^2 log d) field operations."""
    d = 1 << bits
    if prime == 2 or len(values) != d * d:
        raise ValueError("requires odd characteristic and d squared entries")
    work = [sign(a & b) * values[a * d + b] % prime
            for a in range(d) for b in range(d)]
    step = 1
    while step < len(work):
        for block in range(0, len(work), 2 * step):
            for j in range(block, block + step):
                a, b = work[j], work[j + step]
                work[j], work[j + step] = (a + b) % prime, (a - b) % prime
        step *= 2
    inverse = pow(d, -1, prime)
    return tuple(sign(c & e) * inverse * work[e * d + c] % prime
                 for c in range(d) for e in range(d))


def entry(values: tuple[int, ...], bits: int, prime: int,
          x: int, y: int, u: int, v: int) -> int:
    d = 1 << bits
    if x ^ y != u ^ v:
        return 0
    a = x ^ y
    return pow(d, -1, prime) * sum(
        sign(b & (x ^ v)) * values[a * d + b] for b in range(d)) % prime


def physical_matrix(values: tuple[int, ...], bits: int, prime: int,
                    column: bool = False) -> list[list[int]]:
    d = 1 << bits
    if column:
        return [[entry(values, bits, prime, x, y, u, v)
                 for v in range(d) for y in range(d)]
                for x in range(d) for u in range(d)]
    return [[entry(values, bits, prime, x, y, u, v)
             for v in range(d) for u in range(d)]
            for x in range(d) for y in range(d)]


def bell_basis(bits: int, prime: int) -> list[list[int]]:
    d = 1 << bits
    return [[sign(b & x) % prime if y == x ^ a else 0
             for a in range(d) for b in range(d)]
            for x in range(d) for y in range(d)]


def matmul(left: list[list[int]], right: list[list[int]],
           prime: int) -> list[list[int]]:
    columns = list(zip(*right))
    return [[sum(a * b for a, b in zip(row, col)) % prime
             for col in columns] for row in left]


def rank_mod(matrix: list[list[int]], prime: int) -> int:
    a = [row[:] for row in matrix]
    rank = 0
    for j in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(rank, len(a)) if a[i][j] % prime), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = pow(a[rank][j] % prime, -1, prime)
        a[rank] = [x * inverse % prime for x in a[rank]]
        for i in range(rank + 1, len(a)):
            scale = a[i][j] % prime
            if scale:
                a[i] = [(x - scale * y) % prime for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


@lru_cache(None)
def linear_profile_ternary(counts: tuple[int, int, int], length: int) -> tuple[int, ...]:
    """Ranks of coefficients of (1+(z-1)v)^tensor length, from v counts.

    Nonzero scalar factors disappear from the rank. The three local
    eigenpolynomials are 1, z, and -(1+z) over F3. Integer word
    multiplicities are never reduced modulo three.
    """
    a, b, c = counts
    result = [0] * (length + 1)
    for nz in range(length + 1):
        for nb in range(length - nz + 1):
            na = length - nz - nb
            multiplicity = (comb(length, nz) * comb(length - nz, nb)
                            * a ** na * b ** nb * c ** nz)
            if not multiplicity:
                continue
            for j in range(nz + 1):
                if comb(nz, j) % 3:
                    result[nb + j] += multiplicity
    return tuple(result)


def pencil_profiles(values: tuple[int, ...], bits: int, length: int
                    ) -> tuple[tuple[int, ...], tuple[int, ...]]:
    other = transform(values, bits, 3)
    row_counts, col_counts = Counter(values), Counter(other)
    return (linear_profile_ternary(tuple(row_counts[i] for i in range(3)), length),
            linear_profile_ternary(tuple(col_counts[i] for i in range(3)), length))


def direct_pencil_profiles(values: tuple[int, ...], bits: int, length: int
                          ) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Independent literal scalar-eigenvalue products for small alphabets."""
    result = []
    for local in (values, transform(values, bits, 3)):
        ranks = [0] * (length + 1)
        for word in product(local, repeat=length):
            polynomial = [1]
            for value in word:
                new = [0] * (len(polynomial) + 1)
                for j, coefficient in enumerate(polynomial):
                    new[j] = (new[j] + (1 - value) * coefficient) % 3
                    new[j + 1] = (new[j + 1] + value * coefficient) % 3
                polynomial = new
            for j, coefficient in enumerate(polynomial):
                ranks[j] += bool(coefficient)
        result.append(tuple(ranks))
    return tuple(result)


def common_determinant_floor(symbols: int, repetitions: int) -> int:
    t = 2 * symbols * repetitions
    return factorial(t) // (factorial((symbols + 1) * repetitions)
                             * factorial(repetitions) ** (symbols - 1))


def entropy_exceeds(symbols: int, numerator: int, denominator: int) -> bool:
    """Exact gamma_N > numerator/denominator, no floating logarithms."""
    return ((2 * symbols) ** (2 * symbols * denominator)
            > (symbols + 1) ** ((symbols + 1) * denominator)
            * symbols ** (2 * symbols * numerator))
