#!/usr/bin/env python3
"""Exact radial identities, physical Bell controls and complete rank profiles."""

import argparse
from collections import Counter
import json
from math import comb, log2

from bell_crossing import transform
from folded_bell import (coefficient_vectors, coordinates, even_weight,
                         polynomials, radial_identity, rank_profiles)


def check(maximum: int = 24) -> dict:
    rational_identities = sum(radial_identity(k) for k in range(1, maximum + 1))
    full_coefficients = 0
    for k in range(1, 5):
        image, _ = coordinates(k)
        n = 2 * k + 1
        expected = [comb(n, ell) for ell in range(k + 1)]
        histogram = Counter(min(even_weight(v), n - even_weight(v)) for v in image)
        assert [histogram[i] for i in range(k + 1)] == expected
        for prime in (3, 5, 7):
            source = coefficient_vectors(k, prime)
            dual = tuple(transform(row, k, prime) for row in source)
            f = polynomials(k)
            for j, v in enumerate(image):
                ell = min(even_weight(v), n - even_weight(v))
                target = tuple(x % prime for x in reversed(f[k - ell])) + (0,) * ell
                assert tuple(row[j] for row in dual) == target
                full_coefficients += k + 1
            assert all(sum(col) % prime == 1 for col in zip(*source))
            actual_ranks = tuple(sum(x != 0 for x in row) for row in source)
            actual_dual = tuple(sum(x != 0 for x in row) for row in dual)
            assert (actual_ranks, actual_dual) == rank_profiles(k, prime)
    samples = []
    for k in (2, 4, 8, 16, 32, 64, 128, 256):
        row, col = rank_profiles(k, 3)
        maximum_product = max(a * b for a, b in zip(row, col))
        middle = k // 2
        floor = comb(2 * k + 1, middle)
        assert maximum_product >= floor ** 2
        # Exact fixed-scale uniform upper bound at x=3.
        a = sum(comb(2 * k + 1, ell) * 3 ** (k - ell) for ell in range(k + 1))
        b = sum(comb(2 * k + 1, ell) * 3 ** (k - ell) for ell in range(k + 1))
        assert maximum_product * 3 ** k <= a * b
        samples.append({"bits": k, "coefficient_prime": 3,
                        "maximum_rank_product": maximum_product,
                        "central_rank_floor": floor,
                        "finite_rho": log2(maximum_product) / (4 * k)})
    assert 3 ** 5 < 2 ** 8  # H_2(1/4)>4/5>1/sqrt(2).
    assert 2 * 4 ** 2 > 5 ** 2
    return {"status": "passed", "rational_radial_identities": rational_identities,
            "maximum_radial_bits": maximum,
            "full_Bell_coefficient_checks": full_coefficients,
            "physical_dimensions": [2, 4, 8, 16], "coefficient_primes": [3, 5, 7],
            "rank_samples": samples,
            "optimal_leading_rho_in_displayed_coefficient_span":
                "H_2(1/4)=2-(3/4)log2(3)=0.811278124459...",
            "factoring_score": "1/2; unchanged"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum", type=int, default=24)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = check(args.maximum)
    if args.output:
        with open(args.output, "w") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
