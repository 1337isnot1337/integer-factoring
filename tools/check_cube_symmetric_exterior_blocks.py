#!/usr/bin/env python3
"""Exact exterior blocks and natural Cayley profile exclusions through twelve axes.

The retained pressure certificates use only integer and rational arithmetic.
Specialization gives an optimistic lower bound on the generic profile method.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from math import factorial
from pathlib import Path
import random

from cube_symmetric_exterior_blocks import (
    canonical_vectors, exterior_product, linear_form, multiplication_matrix,
    natural_cayley_profile, series_smith_orders, symmetric_blocks, symmetric_profile,
)
from check_multiaxis_gaussian_profile_floor import (
    BITS, independent_enclosure_controls, softplus_negative,
)


def binary_rank(rows: list[int]) -> int:
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return len(basis)


def canonical_form_controls() -> dict:
    identities = 0
    for r in range(1, 13):
        u, pairs = canonical_vectors(r)
        vectors = [value for pair in pairs for value in pair]
        if r & 1:
            vectors.append(u)
        assert binary_rank(vectors) == r
        even = [1]
        for first, second in pairs:
            pair = exterior_product(linear_form(first), linear_form(second))
            even.append(0)
            for j in range(len(even) - 1, 0, -1):
                even[j] ^= exterior_product(even[j - 1], pair)
        for j, polynomial in enumerate(even):
            expected = sum(1 << mask for mask in range(1 << r) if mask.bit_count() == 2 * j)
            assert polynomial == expected
            expected_odd = sum(1 << mask for mask in range(1 << r) if mask.bit_count() == 2 * j + 1)
            assert exterior_product(linear_form(u), polynomial) == expected_odd
            identities += 2
    return {"generator_dimensions": list(range(1, 13)), "complete_form_identities": identities}


def toeplitz_rank(matrix: list[list[int]], precision: int) -> int:
    n = len(matrix)
    rows = [0] * (n * precision)
    for i in range(n):
        for j in range(n):
            polynomial = matrix[i][j]
            while polynomial:
                bit = polynomial & -polynomial
                degree = bit.bit_length() - 1
                for power in range(degree, precision):
                    rows[i * precision + power] ^= 1 << (j * precision + power - degree)
                polynomial ^= bit
    return binary_rank(rows)


def profile_controls() -> dict:
    rng = random.Random(931037)
    full_controls = toeplitz_controls = 0
    for r in range(1, 8):
        degree = 2 + r % 5
        for _ in range(4):
            coefficients = [1 << degree] + [rng.randrange(1 << (degree + 1)) for _ in range(r)]
            full = multiplication_matrix([coefficients[mask.bit_count()] for mask in range(1 << r)])
            literal = series_smith_orders(full, 2 * degree + 1)
            blocked = symmetric_profile(coefficients, 2 * degree + 1)
            assert Counter(literal) == blocked
            assert sum(literal) == (1 << r) * degree
            assert Counter(literal) == Counter(2 * degree - order for order in literal)
            full_controls += 1
            if r <= 4:
                for precision in (1, degree, degree + 1, 2 * degree + 1):
                    expected = sum(max(precision - order, 0) for order in literal)
                    assert toeplitz_rank(full, precision) == expected
                    toeplitz_controls += 1
    return {"full_matrix_profile_comparisons": full_controls,
            "independent_binary_toeplitz_ranks": toeplitz_controls}


def pressure_certificate(d: int, degree: int, profile: Counter, rho: Fraction) -> dict:
    m = 1 << (d - 1)
    assert min(profile) * d == degree

    @lru_cache(None)
    def lhs_lower(x):
        return rho.denominator * sum(count * softplus_negative(x * Fraction(order, degree))[0]
                                     for order, count in profile.items())

    @lru_cache(None)
    def rhs_upper(x):
        return m * rho.numerator * softplus_negative(x / (d * rho))[1]

    limit = 64
    tail_exponent = limit * (1 - rho) / (d * rho)
    assert sum(tail_exponent ** j / factorial(j) for j in range(65)) > 2 * m * rho
    stack, leaves, gaps = [(0, 0)], [], []
    while stack:
        depth, index = stack.pop()
        left = Fraction(limit * index, 1 << depth)
        right = Fraction(limit * (index + 1), 1 << depth)
        gap = lhs_lower(right) - rhs_upper(left)
        if gap > 0:
            leaves.append((depth, index))
            gaps.append(gap)
        else:
            assert depth < 24, (d, depth, index)
            stack.extend(((depth + 1, 2 * index + 1), (depth + 1, 2 * index)))
    maximum_depth = max(depth for depth, _ in leaves)
    previous = 0
    for depth, index in leaves:
        assert index << (maximum_depth - depth) == previous
        previous = (index + 1) << (maximum_depth - depth)
    assert previous == 1 << maximum_depth
    threshold_gap = d * rho.numerator ** 2 - rho.denominator ** 2
    assert threshold_gap > 0
    packed = json.dumps(leaves, separators=(",", ":")).encode()
    return {"axes": d, "degree": degree,
            "specialized_profile": [[order, count] for order, count in sorted(profile.items())],
            "strict_pressure_floor": [rho.numerator, rho.denominator],
            "threshold_integer_gap": threshold_gap, "compact_interval": [0, limit],
            "cover": leaves, "cover_sha256": hashlib.sha256(packed).hexdigest(),
            "cover_intervals": len(leaves), "maximum_depth": maximum_depth,
            "minimum_scaled_interval_gap": min(gaps),
            "tail_exponent": [tail_exponent.numerator, tail_exponent.denominator],
            "tail_positive_taylor_terms": 65}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    controls = {"canonical_forms": canonical_form_controls(), "smith_profiles": profile_controls(),
                "transcendental_enclosures": independent_enclosure_controls()}
    floors = [(2, Fraction(3, 4)), (3, Fraction(3, 5)), (4, Fraction(11, 20)),
              (5, Fraction(1, 2)), (6, Fraction(19, 40)), (7, Fraction(9, 20)),
              (8, Fraction(43, 100)), (9, Fraction(2, 5)), (10, Fraction(2, 5)),
              (11, Fraction(19, 50)), (12, Fraction(3, 8))]
    certificates = []
    for d, rho in floors:
        degree, profile = natural_cayley_profile(d)
        certificate = pressure_certificate(d, degree, profile, rho)
        coefficients = [1 << (degree // (j + 1)) for j in range(d)]
        certificate["full_matrix_order"] = 1 << (d - 1)
        certificate["block_orders_and_multiplicities"] = [
            [len(local), multiplicity] for multiplicity, local in symmetric_blocks(coefficients)]
        certificates.append(certificate)
    result = {"status": "passed", "scope": "exact symmetric exterior blocks; natural Cayley common-scale profile exclusion through twelve axes",
              "fixed_point_bits": BITS, "controls": controls, "certificates": certificates,
              "factoring_score": "1/2 unchanged",
              "generic_profile_scope": "specialization supplies only an optimistic lower pressure bound"}
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(json.dumps({"status": "passed", "dimensions": [d for d, _ in floors],
                      "cover_intervals": [c["cover_intervals"] for c in certificates],
                      "controls": controls, "factoring_score": "1/2 unchanged"}, sort_keys=True))


if __name__ == "__main__":
    main()
