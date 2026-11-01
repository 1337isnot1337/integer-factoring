#!/usr/bin/env python3
"""Independent exact checks of the lattice terminal and its complete search."""

from __future__ import annotations

import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path
import random

from circuit_box_ratio_lattice import (
    lattice_factor_rectangle,
    ratio_lattice_factors,
    reduced_affine_lattice,
    smallest_factor_by_ratio_lattice,
)


def _smallest_factor(number):
    if number % 2 == 0:
        return 2
    for divisor in range(3, isqrt(number) + 1, 2):
        if number % divisor == 0:
            return divisor
    return number


def run_checks(large=True):
    rng = random.Random(2026110151921)
    counts = {
        "weighted_lattice_reductions": 0,
        "literal_ratio_factor_lists": 0,
        "literal_complete_rectangles": 0,
        "smallest_prime_factor_controls": 0,
        "known_factorization_large_class_controls": 0,
    }
    for _ in range(1000):
        modulus = 1 << rng.randrange(1, 101)
        slope = rng.randrange(modulus) | 1
        width = rng.randrange(1, 1 << rng.randrange(1, 101))
        height = rng.randrange(1, 1 << rng.randrange(1, 101))
        result = reduced_affine_lattice(modulus, slope, width, height)
        first, second = result["first"], result["second"]
        assert first[0] * second[1] - first[1] * second[0] == modulus
        assert (first[0] + slope * first[1]) % modulus == 0
        assert (second[0] + slope * second[1]) % modulus == 0
        norm_first, norm_second = result["first_norm"], result["second_norm"]
        inner = result["inner_product"]
        determinant_square = (width * height * modulus) ** 2
        assert norm_first * norm_second - inner * inner == determinant_square
        assert 2 * abs(inner) <= norm_first
        assert 4 * norm_second >= 3 * norm_first
        assert norm_first * norm_first <= 2 * determinant_square
        swaps = result["swaps"]
        assert norm_first * 4**swaps <= (height * modulus) ** 2 * 3**swaps
        assert result["reductions"] == swaps + 1
        counts["weighted_lattice_reductions"] += 1

    maximum_terminal_lines = 0
    for _ in range(1200):
        first = tuple(sorted(rng.randrange(1, 101) for _ in range(2)))
        second = tuple(sorted(rng.randrange(1, 101) for _ in range(2)))
        number = rng.randrange(1, 10001) | 1
        precision = rng.randrange(1, 8)
        ratio = rng.randrange(1, 1 << precision, 2)
        exact = sorted(
            (a, b)
            for a in range(*first)
            for b in range(*second)
            if a * b == number and (a - ratio * b) % (1 << precision) == 0
        )
        actual = ratio_lattice_factors(number, first, second, precision, ratio)
        assert actual["pairs"] == exact
        counts["literal_ratio_factor_lists"] += 1
        whole = lattice_factor_rectangle(number, first, second)
        expected = sorted((a, b) for a in range(*first) for b in range(*second) if a * b == number)
        assert whole["count"] == len(expected)
        assert whole["least_row"] == (expected[0][0] if expected else None)
        assert whole["maximum_lines_per_root"] <= 5
        maximum_terminal_lines = max(maximum_terminal_lines, whole["maximum_lines_per_root"])
        counts["literal_complete_rectangles"] += 1

    for number in range(2, 1200):
        assert smallest_factor_by_ratio_lattice(number)["factor"] == _smallest_factor(number)
        counts["smallest_prime_factor_controls"] += 1

    known_factorizations = []
    for precision in (64, 128, 256, 512):
        exponent_a, exponent_b = 3 * precision // 4, precision // 2
        powers_a = [3**i for i in range(exponent_a + 1)]
        powers_b = [5**i for i in range(exponent_b + 1)]
        a, b = powers_a[-1], powers_b[-1]
        number = a * b
        first = (a - a // 5, a + a // 5 + 1)
        second = (b - b // 5, b + b // 5 + 1)
        modulus = 1 << precision
        assert (first[1] - 1) * (second[1] - 1) <= modulus**3
        ratio = a * pow(b, -1, modulus) % modulus
        expected = []
        for i in range(exponent_a + 1):
            for j in range(exponent_b + 1):
                row = powers_a[i] * powers_b[j]
                if not first[0] <= row < first[1]:
                    continue
                column = powers_a[exponent_a - i] * powers_b[exponent_b - j]
                if second[0] <= column < second[1] and (row - ratio * column) % modulus == 0:
                    expected.append((row, column))
        result = ratio_lattice_factors(number, first, second, precision, ratio)
        assert result["pairs"] == sorted(expected)
        assert result["maximum_lines_per_root"] <= 5
        counts["known_factorization_large_class_controls"] += 1
        known_factorizations.append({
            "ratio_precision": precision,
            "input_bits": number.bit_length(),
            "factorization": {"3": exponent_a, "5": exponent_b},
            "first_interval_length_bits": (first[1] - first[0]).bit_length(),
            "second_interval_length_bits": (second[1] - second[0]).bit_length(),
            "independent_divisor_candidates": (exponent_a + 1) * (exponent_b + 1),
            "count": result["count"],
            "line_queries": result["line_queries"],
            "maximum_lines_per_root": result["maximum_lines_per_root"],
            "lattice_swaps": result["lattice_swaps"],
        })

    empty_inputs = [
        (15, 16411, 262147),
        (18, 131101, 2097169),
        (21, 1048583, 16777259),
        (24, 8388617, 134217757),
    ]
    if not large:
        empty_inputs = empty_inputs[:2]
    empty_rectangles = []
    for exponent, prime_a, prime_b in empty_inputs:
        assert _smallest_factor(prime_a) == prime_a
        assert _smallest_factor(prime_b) == prime_b
        scale = 1 << exponent
        number = prime_a * prime_b
        first, second = (scale, 3 * scale), (3 * scale, 9 * scale)
        assert all(not first[0] <= divisor < first[1] for divisor in (1, prime_a, prime_b, number))
        result = lattice_factor_rectangle(number, first, second)
        assert result["count"] == 0 and result["least_row"] is None
        assert result["line_queries"] <= 20 * result["ratio_queries"]
        empty_rectangles.append({
            "number": number,
            "input_bits": number.bit_length(),
            "prime_factors": [prime_a, prime_b],
            "first_interval": first,
            "second_interval": second,
            **result,
        })

    complete_factorings = []
    for a, b in ((1031, 1091), (8209, 8293), (65537, 65579)):
        assert _smallest_factor(a) == a and _smallest_factor(b) == b
        number = a * b
        result = smallest_factor_by_ratio_lattice(number)
        assert result["factor"] == a
        complete_factorings.append({"number": number, "input_bits": number.bit_length(), "prime_factors": [a, b], **result})

    root = Path(__file__).resolve().parent.parent
    paths = [
        "tools/circuit_box_kloosterman_blocks.py",
        "tools/circuit_box_ratio_lattice.py",
        "tools/check_circuit_box_ratio_lattice.py",
    ]
    return {
        "status": "passed",
        "seed": 2026110151921,
        "large_controls_included": large,
        "counts": counts,
        "maximum_terminal_lines_per_root_in_small_controls": maximum_terminal_lines,
        "known_factorization_controls": known_factorizations,
        "empty_semiprime_rectangles": empty_rectangles,
        "complete_factorings": complete_factorings,
        "proved_search_cost": "2**(n/3)*poly(n) time, polynomial bit space, for positive rectangles with maximum product <=C*N for fixed C",
        "global_comparison": "Classical cube-root scale; weaker than established general factoring bounds and the repository score 1/2 theorem.",
        "source_sha256": {path: hashlib.sha256((root / path).read_bytes()).hexdigest() for path in paths},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--small-only", action="store_true")
    args = parser.parse_args()
    result = run_checks(not args.small_only)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        **result["counts"],
        "empty_rectangles": [
            {key: entry[key] for key in ("input_bits", "ratio_queries", "line_queries", "count")}
            for entry in result["empty_semiprime_rectangles"]
        ],
        "complete_factorings": result["complete_factorings"],
    }))


if __name__ == "__main__":
    main()
