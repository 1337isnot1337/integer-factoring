#!/usr/bin/env python3
"""Exact finite modular Casimir crossings and all-scale pressure exclusions.

The tested spins have n=p-1, p in {3,5,7,11}. Rational projectors are
constructed before reduction modulo p. Each actual coefficient is checked
for a denominator coprime to p. No modular spectral projector is assumed.
"""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, factorial, lcm
from pathlib import Path

from check_all_field_pascal_crossings import smith_orders, trim, weighted_orders


def matrix_product(left, right):
    size = len(left)
    return [[sum((left[i][k] * right[k][j] for k in range(size)), Fraction(0))
             for j in range(size)] for i in range(size)]


def polynomial_product(left, right):
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def rational_casimir_crossing(n):
    """Return the explicit rational polynomial R(u), indexed by (x,y,u,v)."""
    values, denominator = {}, 1
    for occupation in range(2 * n + 1):
        labels = [(i, occupation - i)
                  for i in range(max(0, occupation - n), min(n, occupation) + 1)]
        size = len(labels)
        positions = {label: i for i, label in enumerate(labels)}
        casimir = [[Fraction(0)] * size for _ in range(size)]
        for column, (i, j) in enumerate(labels):
            casimir[column][column] = Fraction(
                (n - occupation) * (n - occupation + 1)
                + i * (n - i + 1) + j * (n - j + 1))
            if i and j < n:
                casimir[positions[i - 1, j + 1]][column] += i * (n - j)
            if j and i < n:
                casimir[positions[i + 1, j - 1]][column] += j * (n - i)
        entries = [[[Fraction(0)] * (n + 1) for _ in range(size)] for _ in range(size)]
        for spin in range(abs(n - occupation), n + 1):
            projector = [[Fraction(i == j) for j in range(size)] for i in range(size)]
            divisor = 1
            for other in range(abs(n - occupation), n + 1):
                if other == spin:
                    continue
                factor = [[casimir[i][j] - other * (other + 1) * (i == j)
                           for j in range(size)] for i in range(size)]
                projector = matrix_product(projector, factor)
                divisor *= spin * (spin + 1) - other * (other + 1)
            eigenvalue = [Fraction(1, factorial(n))]
            for k in range(1, n + 1):
                eigenvalue = polynomial_product(
                    eigenvalue, [Fraction(k if k <= spin else -k), Fraction(1)])
            for i in range(size):
                for j in range(size):
                    coefficient = projector[i][j] / divisor
                    for degree, value in enumerate(eigenvalue):
                        entries[i][j][degree] += coefficient * value
        for i, (x, y) in enumerate(labels):
            for j, (u, v) in enumerate(labels):
                value = tuple(entries[i][j])
                assert value[0] == (x == v and y == u)
                assert value[-1] == Fraction(x == u and y == v, factorial(n))
                values[x, y, u, v] = value
                for coefficient in value:
                    denominator = lcm(denominator, coefficient.denominator)
    return values, denominator


def modular_crossing(rational, prime, spectral_scale=1):
    """T(z)=R(c(z-1)), with the literal swap endpoint T(1)."""
    n = prime - 1
    result = {}
    for key, coefficients in rational.items():
        modular = [value.numerator * pow(value.denominator, -1, prime) % prime
                   for value in coefficients]
        shifted = [sum(modular[i] * pow(spectral_scale, i, prime)
                       * comb(i, degree) * (-1)**(i - degree)
                       for i in range(degree, n + 1))
                   for degree in range(n + 1)]
        result[key] = trim(shifted, prime)
    return result


def zero_blocks(local, prime):
    n = prime - 1
    result = []
    for occupation in range(2 * n + 1):
        labels = list(range(max(0, occupation - n), min(n, occupation) + 1))
        matrix = [[local[x, occupation - x, u, occupation - u]
                   for u in labels] for x in labels]
        orders, precision = smith_orders(matrix, prime)
        result.append((occupation, orders, precision))
    return result


def infinity_profile(local, prime, gauge):
    n = prime - 1
    result = []
    for difference in range(-n, n + 1):
        left = list(range(max(0, -difference), min(n, n - difference) + 1))
        right = list(range(max(0, difference), min(n, n + difference) + 1))
        matrix = [[local[u + difference, v - difference, u, v]
                   for v in right] for u in left]
        orders, _ = weighted_orders(
            matrix, [gauge * u for u in left], [gauge * v for v in right],
            prime, infinity=True)
        result.extend(orders)
    return tuple(sorted(result))


def pressure_polynomial(zero, infinity, energies):
    """x is s^(-1/3); return x^c(B(s)-Z(s)^3) with nonzero constant."""
    coefficients = Counter()
    for first, first_count in Counter(zero).items():
        for second, second_count in Counter(infinity).items():
            coefficients[3 * (first + second)] += first_count * second_count
    for a in energies:
        for b in energies:
            for c in energies:
                coefficients[a + b + c] -= 1
    coefficients = {degree: value for degree, value in coefficients.items() if value}
    minimum, maximum = min(coefficients), max(coefficients)
    return tuple(coefficients.get(degree, 0) for degree in range(minimum, maximum + 1)), minimum


def bernstein_positive(coefficients, maximum_degree=8192):
    """All signs are checked as exact integers; elevation is only a search."""
    support = [(j, value) for j, value in enumerate(coefficients) if value]
    degree = len(coefficients) - 1
    elevated = max(1, degree)
    while elevated <= maximum_degree:
        # b_k * binom(N,k) = sum_j c_j binom(N-j,k-j).
        scaled = [sum(value * comb(elevated - j, k - j)
                      for j, value in support if j <= k)
                  for k in range(elevated + 1)]
        if min(scaled) > 0:
            return {
                "polynomial_degree": degree, "bernstein_degree": elevated,
                "strictly_positive_coefficients": len(scaled),
                "scaled_coefficient_sha256": sha256(
                    json.dumps(scaled, separators=(",", ":")).encode()).hexdigest(),
            }
        elevated *= 2
    raise ArithmeticError("the proposed pressure exclusion has no certificate at the tested elevations")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--primes", nargs="+", type=int, default=(3, 5, 7, 11))
    args = parser.parse_args()
    records, endpoints, pressure_controls, full_cut_controls = [], 0, 0, 0
    for prime in args.primes:
        assert prime in (3, 5, 7, 11)
        rational, denominator = rational_casimir_crossing(prime - 1)
        assert denominator == factorial(prime - 1) and denominator % prime
        reference = modular_crossing(rational, prime)
        infinity = {gauge: infinity_profile(reference, prime, gauge)
                    for gauge in range(-2, 3)}
        shifts = []
        for spectral_scale in range(1, prime):
            local = modular_crossing(rational, prime, spectral_scale)
            for x in range(prime):
                for y in range(prime):
                    for u in range(prime):
                        for v in range(prime):
                            assert sum(local.get((x, y, u, v), ())) % prime == (x == v and y == u)
                            endpoints += 1
            blocks = zero_blocks(local, prime)
            assert infinity_profile(local, prime, 0) == infinity[0]
            if prime in (3, 5):
                # Independently flatten the literal four-index tensor. These
                # controls do not use either weight-block indexing formula.
                row_cut = [[()] * (prime * prime) for _ in range(prime * prime)]
                column_cut = [[()] * (prime * prime) for _ in range(prime * prime)]
                for (x, y, u, v), value in local.items():
                    row_cut[x * prime + y][u * prime + v] = value
                    column_cut[x * prime + u][y * prime + v] = value
                for gauge in range(-2, 3):
                    actual_zero, _ = weighted_orders(
                        row_cut, [gauge * (x + y) for x in range(prime) for y in range(prime)],
                        [0] * (prime * prime), prime)
                    expected_zero = tuple(sorted(order + gauge * occupation
                                                for occupation, orders, _ in blocks for order in orders))
                    assert actual_zero == expected_zero
                    actual_infinity, _ = weighted_orders(
                        column_cut, [gauge * u for x in range(prime) for u in range(prime)],
                        [gauge * v for y in range(prime) for v in range(prime)], prime, infinity=True)
                    assert actual_infinity == infinity[gauge]
                    full_cut_controls += 2
            gauges = []
            for zero_gauge in range(-2, 3):
                zero = tuple(sorted(order + zero_gauge * occupation
                                    for occupation, orders, _ in blocks for order in orders))
                for infinity_gauge in range(-2, 3):
                    energies = [2 * (zero_gauge + infinity_gauge) * j for j in range(prime)]
                    polynomial, minimum = pressure_polynomial(
                        zero, infinity[infinity_gauge], energies)
                    certificate = bernstein_positive(polynomial)
                    certificate.update(zero_gauge=zero_gauge, infinity_gauge=infinity_gauge,
                                       laurent_minimum=minimum)
                    gauges.append(certificate)
                    pressure_controls += 1
            shifts.append({
                "spectral_scale": spectral_scale,
                "zero_blocks": [{"occupation": k, "orders": orders, "precision": precision}
                                for k, orders, precision in blocks],
                "all_scale_gauge_certificates": gauges,
            })
        records.append({"prime": prime, "spin_parameter": prime - 1,
                        "coefficient_denominator_lcm": denominator,
                        "infinity_profiles": infinity, "shifts": shifts})
    output = {
        "scope": "Finite exact crossing family; natural integer gauges in [-2,2]; every real scale s>1",
        "certified_pressure_lower_bound": [3, 4],
        "factoring_lower_bound_claimed": False,
        "literal_endpoint_entries": endpoints,
        "independent_full_cut_controls": full_cut_controls,
        "all_scale_pressure_certificates": pressure_controls,
        "records": records,
    }
    encoded = json.dumps(output, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
