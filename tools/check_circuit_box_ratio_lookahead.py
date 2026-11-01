#!/usr/bin/env python3
"""Exact controls for positive ratio lookahead and its complete factoring search."""

from __future__ import annotations

import argparse
import hashlib
import json
from math import isqrt
from pathlib import Path
import random

from circuit_box_ratio_lookahead import (
    adaptive_factor_rectangle,
    ratio_lookahead_count,
    smallest_factor_by_ratio_lookahead,
)


def _smallest_factor(number: int) -> int:
    if number % 2 == 0:
        return 2
    for candidate in range(3, isqrt(number) + 1, 2):
        if number % candidate == 0:
            return candidate
    return number


def _next_prime(start: int) -> int:
    candidate = max(3, start | 1)
    while _smallest_factor(candidate) != candidate:
        candidate += 2
    return candidate


def _tau_sieve(stop: int) -> list[int]:
    result = [0] * (stop + 1)
    for divisor in range(1, stop + 1):
        for multiple in range(divisor, stop + 1, divisor):
            result[multiple] += 1
    return result


def run_checks(large: bool = True) -> dict:
    rng = random.Random(2026110148801)
    counts = {
        "literal_ratio_lookahead_counts": 0,
        "positive_child_nesting_checks": 0,
        "literal_adaptive_factor_rectangles": 0,
        "exact_level_profiles": 0,
        "positive_product_divisor_bounds": 0,
        "smallest_prime_factor_controls": 0,
        "huge_interval_affine_counts": 0,
    }
    for _ in range(800):
        precision = rng.randrange(1, 10)
        modulus = 1 << precision
        ratio = rng.randrange(1, modulus, 2)
        number = rng.randrange(-1000, 1000) | 1
        first = tuple(sorted(rng.randrange(-20, 41) for _ in range(2)))
        second = tuple(sorted(rng.randrange(-20, 41) for _ in range(2)))
        expected = sum(
            (a * b - number) % (modulus * modulus) == 0
            and (a - ratio * b) % modulus == 0
            for a in range(*first)
            for b in range(*second)
        )
        parent = ratio_lookahead_count(number, first, second, precision, ratio)
        assert parent["count"] == expected
        assert parent["root_count"] <= 4 and parent["affine_pair_queries"] <= 4
        counts["literal_ratio_lookahead_counts"] += 1
        children = sum(
            ratio_lookahead_count(number, first, second, precision + 1, child)["count"]
            for child in (ratio, ratio + modulus)
        )
        next_expected = sum(
            (a * b - number) % (4 * modulus * modulus) == 0
            and (a - ratio * b) % modulus == 0
            for a in range(*first)
            for b in range(*second)
        )
        assert children == next_expected <= parent["count"]
        counts["positive_child_nesting_checks"] += 1

    divisors = _tau_sieve(6400)
    for _ in range(400):
        first = tuple(sorted(rng.randrange(1, 81) for _ in range(2)))
        second = tuple(sorted(rng.randrange(1, 81) for _ in range(2)))
        number = rng.randrange(1, 6400) | 1
        expected = sum(a * b == number for a in range(*first) for b in range(*second))
        result = adaptive_factor_rectangle(number, first, second)
        assert result["count"] == expected
        counts["literal_adaptive_factor_rectangles"] += 1
        if result["terminal_ratio_precision"] == 0:
            continue
        assert result["maximum_stack_size"] <= result["terminal_ratio_precision"]
        previous_live = 0
        for level in result["levels"]:
            precision = level["ratio_precision"]
            modulus = 1 << precision
            pairs = [
                (a, b)
                for a in range(*first)
                for b in range(*second)
                if (a * b - number) % (modulus * modulus) == 0
            ]
            ratios = {a * pow(b, -1, modulus) % modulus for a, b in pairs}
            assert level["congruence_pairs"] == len(pairs)
            assert level["live"] == len(ratios)
            assert level["visited"] == (1 if precision == 1 else 2 * previous_live)
            previous_live = level["live"]
            counts["exact_level_profiles"] += 1
            minimum, maximum = result["minimum_product"], result["maximum_product"]
            first_product = minimum + (number - minimum) % (modulus * modulus)
            positive_products = range(first_product, maximum + 1, modulus * modulus)
            assert all(product > 0 for product in positive_products)
            divisor_bound = sum(divisors[product] for product in positive_products)
            assert level["live"] <= len(pairs) <= divisor_bound
            counts["positive_product_divisor_bounds"] += 1

    for number in range(2, 1200):
        actual = smallest_factor_by_ratio_lookahead(number)
        assert actual["factor"] == _smallest_factor(number)
        counts["smallest_prime_factor_controls"] += 1

    huge = []
    for precision in (20, 80, 200, 512):
        modulus = 1 << precision
        copies_a = (1 << 120) + 3
        copies_b = (1 << 101) + 5
        first = (-17, -17 + copies_a * modulus * modulus)
        second = (23, 23 + copies_b * modulus * modulus)
        number = 5
        ratio = 5
        expected = 4 * modulus * copies_a * copies_b
        result = ratio_lookahead_count(number, first, second, precision, ratio)
        assert result["count"] == expected
        counts["huge_interval_affine_counts"] += 1
        huge.append({
            "ratio_precision": precision,
            "row_length_bits": (first[1] - first[0]).bit_length(),
            "column_length_bits": (second[1] - second[0]).bit_length(),
            "count": result["count"],
            "affine_pair_queries": result["affine_pair_queries"],
            "control": "four low roots; R affine solutions among R^2 high-bit pairs per full product period",
        })

    empty_rectangles = []
    scales = (10, 12, 15, 18, 21, 24) if large else (10, 12, 15)
    for scale_exponent in scales:
        scale = 1 << scale_exponent
        prime_small = _next_prime(scale // 2 + 3)
        prime_large = _next_prime(8 * scale + 3)
        assert scale // 2 < prime_small < scale
        assert 8 * scale < prime_large < 16 * scale
        number = prime_small * prime_large
        first = (scale, 3 * scale)
        second = (3 * scale, 9 * scale)
        result = adaptive_factor_rectangle(number, first, second)
        assert result["count"] == 0
        assert result["minimum_product"] < number < result["maximum_product"]
        assert result["maximum_product"] < 7 * number
        assert all(not first[0] <= divisor < first[1] for divisor in (1, prime_small, prime_large, number))
        empty_rectangles.append({
            "scale_exponent": scale_exponent,
            "number": number,
            "input_bits": number.bit_length(),
            "prime_factors": [prime_small, prime_large],
            "primality_certificate": "complete deterministic trial division through each integer square root",
            "first_interval": first,
            "second_interval": second,
            "fixed_terminal_ratio_classes": 1 << (result["terminal_ratio_precision"] - 1),
            "earlier_corner_tree_states_certified_through_depth_s_minus_3": scale // 8 - 1,
            **result,
        })

    complete_factorings = []
    for start_a, start_b in ((1 << 10, (1 << 10) + 64), (1 << 13, (1 << 13) + 100), (1 << 16, (1 << 16) + 32)):
        prime_a = _next_prime(start_a)
        prime_b = _next_prime(start_b)
        number = prime_a * prime_b
        result = smallest_factor_by_ratio_lookahead(number)
        assert result["factor"] == prime_a
        complete_factorings.append({
            "number": number,
            "input_bits": number.bit_length(),
            "prime_factors": [prime_a, prime_b],
            **result,
        })

    root = Path(__file__).resolve().parent.parent
    sources = [
        "tools/circuit_box_kloosterman_blocks.py",
        "tools/circuit_box_ratio_lookahead.py",
        "tools/check_circuit_box_ratio_lookahead.py",
    ]
    return {
        "status": "passed",
        "seed": 2026110148801,
        "large_controls_included": large,
        "counts": counts,
        "huge_interval_controls": huge,
        "empty_semiprime_rectangles": empty_rectangles,
        "complete_factorings": complete_factorings,
        "proved_search_cost": "2**(n/3 + O(n/log(n+2))) time, polynomial bit space, on positive rectangles with products at most C*N for fixed C",
        "global_comparison": "Weaker than the repository's exp(O(n/sqrt(log n))) complete factoring theorem; score remains 1/2.",
        "source_sha256": {path: hashlib.sha256((root / path).read_bytes()).hexdigest() for path in sources},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--small-only", action="store_true")
    args = parser.parse_args()
    result = run_checks(not args.small_only)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(json.dumps({
        "status": result["status"],
        **result["counts"],
        "empty_rectangles": [
            {key: item[key] for key in ("input_bits", "visited_nodes", "affine_pair_queries", "count")}
            for item in result["empty_semiprime_rectangles"]
        ],
        "complete_factorings": result["complete_factorings"],
    }))


if __name__ == "__main__":
    main()
