"""LCD binary-code quotients of the folded Bell crossing.

The exact duality, maximal common division and characteristic-zero
pressure restriction are proved in research/code-quotient-crossings.md.
"""

from collections import Counter
from fractions import Fraction
from math import comb, log2

from search_bell_profiles import optimum, pressure


def span(basis: tuple[int, ...]) -> tuple[int, ...]:
    result = [0]
    for row in basis:
        result += [x ^ row for x in result]
    return tuple(result)


def binary_basis(rows: tuple[int, ...]) -> dict[int, int]:
    basis = {}
    for row in rows:
        while row:
            pivot = (row & -row).bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return basis


def nullspace(rows: tuple[int, ...], n: int) -> tuple[int, ...]:
    basis = binary_basis(rows)
    result = []
    for j in range(n):
        if j in basis:
            continue
        value = 1 << j
        for pivot in sorted(basis, reverse=True):
            if (basis[pivot] & value).bit_count() % 2:
                value |= 1 << pivot
        result.append(value)
    return tuple(result)


def multiply(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return tuple(result)


def code_path(n: int, generators: tuple[int, ...], prime: int) -> dict:
    assert n % 2 == 1 and prime != 2
    c = span(generators)
    assert len(set(c)) == len(c) and (1 << n) - 1 in c
    v_basis = nullspace(generators, n)
    v = span(v_basis)
    assert set(c).intersection(v) == {0}
    assert all(x.bit_count() % 2 == 0 for x in v)
    m, h = len(v_basis), max(x.bit_count() for x in v) // 2
    assert m >= 2 and m % 2 == 0 and h >= m // 2
    # Fourier diagonal form after maximal common scalar division.
    spectra = []
    for y in v:
        j = y.bit_count() // 2
        a = tuple(comb(h - j, i) for i in range(h - j + 1))
        b = tuple((-1) ** i * comb(j, i) for i in range(j + 1))
        spectra.append(multiply(a, b))
    coefficients = []
    for degree in range(h + 1):
        coefficients.append(tuple(pow(2, -h, prime) * sum(
            (-1) ** ((x & y).bit_count() % 2) * spectrum[degree]
            for y, spectrum in zip(v, spectra)) % prime for x in v))
    coefficients = tuple(coefficients)
    assert all(sum(col) % prime == 1 for col in zip(*coefficients))
    gaussian = sum((-1) ** (x.bit_count() // 2 % 2) for x in v)
    d = 1 << (m // 2)
    assert abs(gaussian) == d
    dual = tuple(tuple(gaussian // d * pow(d, -1, prime) * sum(
        (-1) ** ((x ^ y).bit_count() // 2 % 2) * values[j]
        for j, y in enumerate(v)) % prime for x in v) for values in coefficients)
    assert dual == tuple(reversed(coefficients))
    orders = tuple(next(i for i in range(h + 1) if coefficients[i][j]) for j in range(len(v)))
    distances = tuple(min((x ^ y).bit_count() for y in c) for x in v)
    assert all(a >= b for a, b in zip(orders, distances))
    histogram = Counter(orders)
    low = tuple(histogram[i] for i in range(h + 1))
    high = tuple(reversed(low))
    _, approximate_x = optimum(low, high)
    x = Fraction(approximate_x).limit_denominator(256)
    bound = pressure(low, high, x)
    return {"ambient_length": n, "code_generators": generators,
            "dual_basis": v_basis, "dual_dimension": m, "physical_bits": m // 2,
            "prime": prime, "degree_after_maximal_common_division": h,
            "common_divisor_degree": (n - 1) // 2 - h,
            "valuation_histogram": low, "distance_histogram":
                tuple(Counter(distances)[i] for i in range(h + 1)),
            "raised_valuation_cosets": sum(a > b for a, b in zip(orders, distances)),
            "scale": [x.numerator, x.denominator],
            "uniform_pressure": [bound.numerator, bound.denominator],
            "rho_pressure_approx": (log2(bound.numerator) - log2(bound.denominator)) / (2 * m),
            "coefficients_by_degree": coefficients}
