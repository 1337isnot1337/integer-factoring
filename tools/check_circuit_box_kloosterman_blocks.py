#!/usr/bin/env python3
"""Independent exact checks for shared dyadic Kloosterman frequency blocks."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random

from circuit_box_kloosterman_blocks import (
    affine_interval_pairs,
    count_inverse_rectangle,
    odd_square_roots,
    stationary_nodes,
    weighted_frequency_block,
)


def _phase(vector: list[int], exponent: int, coefficient: int = 1) -> None:
    """Add coefficient*zeta**exponent in Z[zeta], zeta**(q/2)=-1."""
    half = len(vector)
    exponent %= 2 * half
    if exponent < half:
        vector[exponent] += coefficient
    else:
        vector[exponent - half] -= coefficient


def _literal_kloosterman(first: int, second: int, exponent: int) -> list[int]:
    modulus = 1 << exponent
    vector = [0] * (modulus // 2)
    for node in range(1, modulus, 2):
        _phase(vector, first * node + second * pow(node, -1, modulus))
    return vector


def _literal_block(number, first, second, exponent, slope, precision) -> int:
    """Complete x and c scan, after the independent elementary Ramanujan sum."""
    modulus = 1 << exponent
    answer = 0
    for coefficient in range(slope, modulus, 1 << precision):
        full = [0] * modulus
        half = [0] * (modulus // 2)
        for a in range(*first):
            for b in range(*second):
                full[(a + coefficient * b) % modulus] += 1
                half[(a + coefficient * b) % (modulus // 2)] += 1
        for node in range(1, modulus, 2):
            phase = node + coefficient * number * pow(node, -1, modulus)
            answer += modulus * full[phase % modulus]
            answer -= (modulus // 2) * half[phase % (modulus // 2)]
    return answer


def _raw_cyclotomic_block(number, first, second, exponent, slope, precision):
    modulus = 1 << exponent
    vector = [0] * (modulus // 2)
    for frequency in range(1, modulus, 2):
        for coefficient in range(slope, modulus, 1 << precision):
            for node in range(1, modulus, 2):
                phase = node + coefficient * number * pow(node, -1, modulus)
                for a in range(*first):
                    for b in range(*second):
                        _phase(vector, frequency * (phase - a - coefficient * b))
    return vector


def run_checks() -> dict:
    rng = random.Random(2026110141729)
    counts = {
        "affine_interval_counts": 0,
        "complete_odd_square_root_lists": 0,
        "exact_cyclotomic_stationary_phase_identities": 0,
        "opposite_parity_frequency_vanishing": 0,
        "literal_shared_frequency_blocks": 0,
        "raw_cyclotomic_weighted_blocks": 0,
        "ratio_filtered_conductor_differences": 0,
        "complete_inverse_rectangle_counts": 0,
        "no_wrap_factor_rectangle_counts": 0,
    }
    for _ in range(1000):
        modulus = 1 << rng.randrange(8)
        slope = rng.randrange(-100, 101) | 1
        target = rng.randrange(-200, 201)
        first = tuple(sorted(rng.randrange(-20, 31) for _ in range(2)))
        second = tuple(sorted(rng.randrange(-20, 31) for _ in range(2)))
        expected = sum(
            (a + slope * b - target) % modulus == 0
            for a in range(*first)
            for b in range(*second)
        )
        assert affine_interval_pairs(first, second, slope, target, modulus) == expected
        counts["affine_interval_counts"] += 1

    for exponent in range(1, 10):
        modulus = 1 << exponent
        for value in range(1, min(modulus, 34), 2):
            expected = [root for root in range(1, modulus, 2) if (root * root - value) % modulus == 0]
            assert sorted(odd_square_roots(value, exponent)) == expected
            counts["complete_odd_square_root_lists"] += 1

    maximum_nodes = 0
    for exponent in range(1, 11):
        modulus = 1 << exponent
        for _ in range(32):
            number = rng.randrange(-10000, 10000) | 1
            slope = rng.randrange(modulus) | 1
            frequency = rng.randrange(modulus) | 1
            nodes = stationary_nodes(number, slope, exponent)
            maximum_nodes = max(maximum_nodes, len(nodes))
            compact = [0] * (modulus // 2)
            for node in nodes:
                phase = node + slope * number * pow(node, -1, modulus)
                _phase(compact, frequency * phase, 1 << (exponent // 2))
            assert compact == _literal_kloosterman(frequency, frequency * slope * number, exponent)
            counts["exact_cyclotomic_stationary_phase_identities"] += 1
            if exponent >= 2:
                even_frequency = 2 * rng.randrange(modulus // 2)
                assert not any(_literal_kloosterman(frequency, even_frequency, exponent))
                counts["opposite_parity_frequency_vanishing"] += 1

    for _ in range(240):
        exponent = rng.randrange(1, 10)
        modulus = 1 << exponent
        precision = rng.randrange(max(1, exponent // 2), exponent + 1)
        slope = rng.randrange(1, 1 << precision, 2)
        number = rng.randrange(-1000, 1000) | 1
        left_a = rng.randrange(-modulus, modulus)
        left_b = rng.randrange(-modulus, modulus)
        first = (left_a, left_a + rng.randrange(1, 13))
        second = (left_b, left_b + rng.randrange(1, 13))
        block = weighted_frequency_block(number, first, second, exponent, slope, precision)
        assert block["value"] == _literal_block(number, first, second, exponent, slope, precision)
        assert block["affine_pair_queries"] <= 16
        counts["literal_shared_frequency_blocks"] += 1

    for exponent in range(1, 5):
        for _ in range(12):
            precision = max(1, exponent // 2)
            slope = rng.randrange(1, 1 << precision, 2)
            number = rng.randrange(-50, 50) | 1
            first = (rng.randrange(-5, 5), rng.randrange(5, 8))
            second = (rng.randrange(-5, 5), rng.randrange(5, 8))
            block = weighted_frequency_block(number, first, second, exponent, slope, precision)
            exact = _raw_cyclotomic_block(number, first, second, exponent, slope, precision)
            assert exact == [block["value"]] + [0] * ((1 << (exponent - 1)) - 1)
            counts["raw_cyclotomic_weighted_blocks"] += 1

    for _ in range(600):
        exponent = rng.randrange(2, 11)
        modulus = 1 << exponent
        precision = exponent // 2
        slope = rng.randrange(1, 1 << precision, 2)
        number = rng.randrange(-1000, 1000) | 1
        first = tuple(sorted(rng.randrange(-10, 40) for _ in range(2)))
        second = tuple(sorted(rng.randrange(-10, 40) for _ in range(2)))
        full = half = 0
        for a in range(*first):
            for b in range(*second):
                if (a - slope * b) % (1 << precision) == 0:
                    full += (a * b - number) % modulus == 0
                    half += (a * b - number) % (modulus // 2) == 0
        block = weighted_frequency_block(number, first, second, exponent, slope, precision)
        assert block["value"] == (modulus * modulus // 2) * (2 * full - half)
        counts["ratio_filtered_conductor_differences"] += 1

    for _ in range(400):
        exponent = rng.randrange(1, 10)
        modulus = 1 << exponent
        number = rng.randrange(-1000, 1000) | 1
        first = tuple(sorted(rng.randrange(-20, modulus + 21) for _ in range(2)))
        second = tuple(sorted(rng.randrange(-20, modulus + 21) for _ in range(2)))
        expected = sum(
            (a * b - number) % modulus == 0
            for a in range(*first)
            for b in range(*second)
        )
        assert count_inverse_rectangle(number, first, second, exponent)["count"] == expected
        counts["complete_inverse_rectangle_counts"] += 1

    for number in range(3, 300, 2):
        first = (2, int(number**0.5) + 1)
        second = (2, number // 2 + 1)
        height = max(abs(number - a * b) for a in (first[0], first[1] - 1) for b in (second[0], second[1] - 1))
        exponent = max(2, height.bit_length() + 1)
        expected = sum(number % a == 0 for a in range(*first))
        actual = count_inverse_rectangle(number, first, second, exponent)["count"]
        assert actual == expected
        counts["no_wrap_factor_rectangle_counts"] += 1

    large = []
    for exponent in (127, 256, 513, 1024):
        precision = exponent // 2
        modulus = 1 << exponent
        low_modulus = 1 << precision
        a = (1 << (precision - 2)) + 17
        b = (1 << (precision - 3)) + 29
        slope = a * pow(b, -1, low_modulus) % low_modulus
        first = (a, a + 1)
        second = (b, b + 1)
        records = []
        for delta, expected_sign in ((0, 1), (modulus // 2, -1), (modulus // 4, 0)):
            block = weighted_frequency_block(a * b + delta, first, second, exponent, slope, precision)
            assert block["value"] == expected_sign * modulus * modulus // 2
            records.append({
                "target_difference": delta,
                "value": block["value"],
                "stationary_node_count": len(block["stationary_nodes"]),
                "active_nodes": block["active_nodes"],
                "affine_pair_queries": block["affine_pair_queries"],
            })
        # The adjacent wrong ratio class must vanish even at a true pair.
        wrong = weighted_frequency_block(a * b, first, second, exponent, (slope + 2) % low_modulus, precision)
        assert wrong["value"] == 0
        large.append({
            "conductor_exponent": exponent,
            "slope_precision": precision,
            "slope_residue": slope,
            "row": a,
            "column": b,
            "frequency_pair_count": (modulus // 2) * (1 << (exponent - precision)),
            "records": records,
            "wrong_ratio_value": wrong["value"],
            "independent_control": "singleton ratio class and exact product congruences at m and m-1",
        })

    source_root = Path(__file__).resolve().parent.parent
    paths = [
        "tools/circuit_box_kloosterman_blocks.py",
        "tools/check_circuit_box_kloosterman_blocks.py",
    ]
    return {
        "status": "passed",
        "seed": 2026110141729,
        "counts": counts,
        "maximum_stationary_nodes": maximum_nodes,
        "large_shared_blocks": large,
        "scope": "Polynomial cost per complete shared frequency block; full rectangle count still costs 2**(m/2)*poly(m).",
        "source_sha256": {
            path: hashlib.sha256((source_root / path).read_bytes()).hexdigest()
            for path in paths
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(json.dumps({"status": result["status"], **result["counts"], "large_shared_blocks": len(result["large_shared_blocks"])}))


if __name__ == "__main__":
    main()
