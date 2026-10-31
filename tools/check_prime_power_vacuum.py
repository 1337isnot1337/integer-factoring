#!/usr/bin/env python3
"""Exact Artin--Hasse coefficients and finite rational-localizer diagnostics.

The all-size nonrationality statement uses the cited transcendence theorem.
These computations verify series transcription and finite linear complexity.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


def prime_powers(p: int, limit: int) -> list[int]:
    result, value = [], 1
    while value <= limit:
        result.append(value)
        value *= p
    return result


def artin_hasse_recurrence(p: int, degree: int) -> list[Fraction]:
    """n*a_n = sum_{p^j <= n} a_{n-p^j}, from the logarithmic derivative."""
    powers = prime_powers(p, degree)
    coefficients = [Fraction(1)]
    for n in range(1, degree + 1):
        coefficients.append(sum((coefficients[n - q] for q in powers if q <= n),
                                Fraction()) / n)
    return coefficients


def mobius(n: int) -> int:
    parity, divisor = 0, 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            n //= divisor
            parity ^= 1
            if n % divisor == 0:
                return 0
        divisor += 1
    return -1 if parity ^ (n > 1) else 1


def artin_hasse_euler_product(p: int, degree: int) -> list[Fraction]:
    """Independent Euler product product_(p not | m) (1-x^m)^(-mu(m)/m)."""
    output = [Fraction(1)] + [Fraction() for _ in range(degree)]
    for m in range(1, degree + 1):
        if m % p == 0 or (mu := mobius(m)) == 0:
            continue
        block = [Fraction(1)]
        for k in range(1, degree // m + 1):
            block.append(block[-1] * (Fraction(k - 1) + Fraction(mu, m)) / k)
        previous = output
        output = [sum((previous[n - k * m] * block[k]
                       for k in range(min(n // m, len(block) - 1) + 1)), Fraction())
                  for n in range(degree + 1)]
    return output


def reduce_coefficients(coefficients: list[Fraction], p: int) -> list[int]:
    assert all(value.denominator % p for value in coefficients), "nonintegral coefficient"
    return [(value.numerator % p) * pow(value.denominator, -1, p) % p
            for value in coefficients]


def berlekamp_massey(sequence: list[int], p: int) -> tuple[int, list[int]]:
    current, previous = [1], [1]
    length, shift, saved_discrepancy = 0, 1, 1
    for n in range(len(sequence)):
        discrepancy = sum(current[j] * sequence[n - j]
                          for j in range(len(current)) if j <= n) % p
        if discrepancy == 0:
            shift += 1
            continue
        old_current = current[:]
        ratio = discrepancy * pow(saved_discrepancy, -1, p) % p
        current += [0] * max(0, len(previous) + shift - len(current))
        for j, value in enumerate(previous):
            current[j + shift] = (current[j + shift] - ratio * value) % p
        if 2 * length <= n:
            length = n + 1 - length
            previous, saved_discrepancy, shift = old_current, discrepancy, 1
        else:
            shift += 1
    current += [0] * max(0, length + 1 - len(current))
    assert all(sum(current[j] * sequence[n - j] for j in range(length + 1)) % p == 0
               for n in range(length, len(sequence)))
    return length, current[:length + 1]


def linear_system_consistent(rows: list[list[int]], variables: int, p: int) -> bool:
    """Independent row elimination of an augmented recurrence system."""
    pivots: dict[int, list[int]] = {}
    for original in rows:
        row = [value % p for value in original]
        for column in range(variables):
            if row[column] == 0:
                continue
            if column not in pivots:
                inverse = pow(row[column], -1, p)
                pivots[column] = [value * inverse % p for value in row]
                break
            factor = row[column]
            row = [(value - factor * other) % p
                   for value, other in zip(row, pivots[column])]
        else:
            if row[variables]:
                return False
    return True


def recurrence_feasible(sequence: list[int], length: int, p: int) -> bool:
    return linear_system_consistent(
        [[sequence[n - j] for j in range(1, length + 1)] + [-sequence[n]]
         for n in range(length, len(sequence))], length, p)


def verify(maximum: int = 256) -> dict:
    if maximum < 64:
        raise ValueError("the diagnostic requires at least 64 coefficients")
    fields = (2, 3, 5, 7)
    records = []
    checked_coefficients, feasibility_checks = 0, 0
    lengths = sorted({8, 16, 32, 64, maximum})
    for p in fields:
        coefficients = artin_hasse_recurrence(p, max(maximum - 1, 64))
        assert coefficients[:65] == artin_hasse_euler_product(p, 64)
        checked_coefficients += 65
        reduced = reduce_coefficients(coefficients, p)
        diagnostics = []
        for count in lengths:
            sequence = reduced[:count]
            order, denominator = berlekamp_massey(sequence, p)
            assert recurrence_feasible(sequence, order, p)
            assert order > 0 and not recurrence_feasible(sequence, order - 1, p)
            feasibility_checks += 2
            diagnostics.append({"coefficient_count": count,
                                "minimum_recurrence_length": order,
                                "one_denominator": denominator})
        rational = [pow(1 if p == 2 else 2, j, p) for j in range(maximum)]
        assert berlekamp_massey(rational, p)[0] == 1
        assert recurrence_feasible(rational, 1, p)
        assert not recurrence_feasible(rational, 0, p)
        records.append({"characteristic": p, "first_32_coefficients": reduced[:32],
                        "diagnostics": diagnostics,
                        "geometric_series_control_length": 1})
    return {
        "status": "exact finite-series controls; no new factoring exponent",
        "scope": "A fixed finite polynomial cannot localize the specified Artin-Hasse vacuum factor. Finite linear complexities are diagnostics only, with no claimed asymptotic rate.",
        "external_all_size_input": {
            "author": "Joe Kramer-Miller", "theorem": "Theorem 1.1",
            "source": "https://arxiv.org/abs/2404.06968v1",
            "statement": "The reduction of the Artin-Hasse exponential is transcendental over F_p(x) for every prime p."},
        "independent_euler_product_coefficients_checked": checked_coefficients,
        "independent_recurrence_feasibility_checks": feasibility_checks,
        "fields": records,
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--maximum", type=int, default=256)
    args = parser.parse_args()
    result = verify(args.maximum)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
