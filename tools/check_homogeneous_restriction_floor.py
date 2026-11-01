#!/usr/bin/env python3
"""Certify an all-scale 5/7 pressure floor from exact saved profiles."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import accumulate, product
import json
from pathlib import Path

from reflection_homogeneous_restrictions import middle_fiber_retraction


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "evidence/reflection-homogeneous-restrictions.json"


def multiply(first, second):
    output = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        if a:
            for j, b in enumerate(second):
                if b:
                    output[i + j] += a * b
    return output


def power(polynomial, exponent):
    result = [1]
    while exponent:
        if exponent & 1:
            result = multiply(result, polynomial)
        exponent >>= 1
        if exponent:
            polynomial = multiply(polynomial, polynomial)
    return result


def counting_polynomial(values):
    counts = Counter(values)
    assert values and min(values) >= 0
    return [counts[j] for j in range(max(values) + 1)]


def certify(record):
    zero = counting_polynomial(record["zero_orders"])
    infinity = counting_polynomial(record["infinity_orders"])
    heights = counting_polynomial(record["letter_heights"])
    left = power(multiply(zero, infinity), 7)
    right = power(heights, 20)
    coefficients = [0] * max(5 * (len(left) - 1) + 1,
                             7 * (len(right) - 1) + 1)
    for j, value in enumerate(left):
        coefficients[5 * j] += value
    for j, value in enumerate(right):
        coefficients[7 * j] -= value
    while coefficients[-1] == 0:
        coefficients.pop()
    first = list(accumulate(coefficients))
    second = list(accumulate(first))
    assert min(second) >= 0
    assert first[-1] == 15**28 - 15**20 > 0
    assert second[-1] >= 0
    # Independent evaluations check the polynomial substitutions and powers.
    for x in (1, 2, 3):
        def evaluate(polynomial, value):
            result = 0
            for coefficient in reversed(polynomial):
                result = result * value + coefficient
            return result
        assert evaluate(coefficients, x) == (
            (evaluate(zero, x**5) * evaluate(infinity, x**5))**7
            - evaluate(heights, x**7)**20)
    encoded = json.dumps(coefficients, separators=(",", ":")).encode()
    return {
        "representative_pair": [record["kernel_vector"], record["image_functional"]],
        "polynomial_degree": len(coefficients) - 1,
        "minimum_second_cumulative_coefficient": min(second),
        "terminal_first_sum": first[-1],
        "terminal_second_sum": second[-1],
        "polynomial_sha256": sha256(encoded).hexdigest(),
        "strictly_positive_for_zero_less_than_x_at_most_one": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    source_bytes = SOURCE.read_bytes()
    data = json.loads(source_bytes)
    assert data["status"] == "pass" and data["complete_fiber_enumeration"]
    for key, filename in (
            ("engine_sha256", "reflection_homogeneous_restrictions.py"),
            ("checker_sha256", "check_reflection_homogeneous_restrictions.py")):
        assert data[key] == sha256((ROOT / "tools" / filename).read_bytes()).hexdigest()
    expected_pairs = {(v, f) for v, f in product(range(1, 16), repeat=2)
                      if (v & f).bit_count() % 2}
    seen_pairs, profiles, certificates = set(), set(), []
    for record in data["records"]:
        pair = record["kernel_vector"], record["image_functional"]
        assert pair in expected_pairs and pair not in seen_pairs
        seen_pairs.add(pair)
        embedding, retraction, heights = middle_fiber_retraction(*pair)
        assert list(embedding) == record["embedding"]
        assert list(retraction) == record["retraction"]
        assert list(heights) == record["letter_heights"]
        assert [4 * h for h in heights] == record["type_energies"]
        for key in ("zero_orders", "infinity_orders"):
            assert len(record[key]) == 225
            assert record[key] == sorted(record[key]) and record[key][0] == 0
        profile = tuple(record["zero_orders"]), tuple(record["infinity_orders"])
        if profile not in profiles:
            profiles.add(profile)
            certificates.append(certify(record))
    assert seen_pairs == expected_pairs and len(profiles) == 24
    result = {
        "status": "pass",
        "scope": "All 120 stored homogeneous middle-fiber restrictions, with inherited type energies",
        "rho_floor": "5/7",
        "strictly_above_one_over_sqrt_two": 2 * 5**2 > 7**2,
        "retractions_checked": len(seen_pairs),
        "distinct_profile_certificates": len(certificates),
        "maximum_polynomial_degree": max(c["polynomial_degree"] for c in certificates),
        "certificates": certificates,
        "source_evidence_sha256": sha256(source_bytes).hexdigest(),
        "profile_engine_sha256": data["engine_sha256"],
        "profile_checker_sha256": data["checker_sha256"],
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "profile_matrices_recomputed_by_this_replay": False,
        "floating_arithmetic_used_for_certification": False,
        "factoring_exponent_improved": False,
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
        print(json.dumps({key: result[key] for key in (
            "status", "rho_floor", "retractions_checked",
            "distinct_profile_certificates", "maximum_polynomial_degree")}))
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
