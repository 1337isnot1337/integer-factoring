#!/usr/bin/env python3
"""Finite controls for the primitive-recipient rank envelope.

Checks use exact integers, including orbit determinants and conservative
rational bounds for the Euclidean packing expressions. This is not a
factoring algorithm or a certificate of a larger factoring exponent.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path

from check_factoring_network_structure import Fanout
from check_factoring_primitive_recipients import determinant, inverse_mod, negacyclic_orbit


def integer_root(value: int, degree: int) -> int:
    """The exact floor of the positive integer degree-th root."""
    if value < 0 or degree < 1:
        raise ValueError("nonnegative value and positive degree required")
    low, high = 0, 1 << ((value.bit_length() + degree - 1) // degree)
    while low + 1 < high:
        middle = (low + high) // 2
        if middle ** degree <= value:
            low = middle
        else:
            high = middle
    answer = high if high ** degree <= value else low
    assert answer ** degree <= value < (answer + 1) ** degree
    return answer


def recipient_sets(h: int) -> list[tuple[int, ...]]:
    odds = tuple(range(1, 2 * h, 2))
    if h <= 4:
        return [indices for size in range(h + 1)
                for indices in combinations(odds, size)]
    return [(), (1,), (3,), (1, 3), (1, 7), (1, 9),
            (1, 3, 7), (1, 3, 5, 7), (1, 5, 9, 13), odds]


def image_checks(fan: Fanout, totals: Counter) -> list[dict]:
    sets = recipient_sets(fan.h)
    # Compute the original source messages once, then project onto each set.
    odds = tuple(range(1, fan.a, 2))
    rows = [tuple(fan.message_at(digits, j) for j in odds)
            for digits in product(range(fan.base), repeat=fan.h)]
    records = []
    sqrt_upper = isqrt(fan.h)
    sqrt_upper += sqrt_upper * sqrt_upper < fan.h
    envelope_factor = (1 + sqrt_upper) ** fan.h
    for indices in sets:
        positions = tuple((j - 1) // 2 for j in indices)
        fibers = Counter(tuple(row[i] for i in positions) for row in rows)
        rank, largest = len(fibers), max(fibers.values())
        cap = min(fan.base ** fan.h, fan.modulus ** len(indices))
        assert rank <= cap <= rank * envelope_factor

        # Replacing delta by floor(delta) and sqrt(h) by ceil(sqrt(h))
        # weakens the packing estimate and makes it exactly integral.
        separation = integer_root(fan.modulus ** len(indices), fan.h)
        assert largest * separation ** fan.h <= (
            separation + (fan.base - 1) * sqrt_upper) ** fan.h
        records.append({
            "h": fan.h,
            "digit_bits": fan.digit_bits,
            "source_parameters": fan.source_parameters,
            "recipients": list(indices),
            "image_size": rank,
            "maximum_fiber": largest,
            "cardinality_cap": str(cap),
            "integer_envelope_factor": envelope_factor,
            "integer_separation_lower_bound": str(separation),
        })
        totals["exhaustive_images"] += 1
        totals["enumerated_input_projections"] += len(rows)
    return records


def kernel_basis(fan: Fanout, indices: tuple[int, ...]) -> list[list[int]]:
    count, modulus = len(indices), fan.modulus
    evaluation = [[pow(fan.root, i * j, modulus) for i in range(fan.h)]
                  for j in indices]
    pivot_inverse = inverse_mod([row[:count] for row in evaluation], modulus)
    rows = []
    for i in range(fan.h):
        vector = [0] * fan.h
        if i < count:
            vector[i] = modulus
        else:
            vector[i] = 1
            for j in range(count):
                vector[j] = -sum(pivot_inverse[j][a] * evaluation[a][i]
                                 for a in range(count)) % modulus
        rows.append(vector)
    assert abs(determinant(rows)) == modulus ** count
    return rows


def verify_orbit(fan: Fanout, indices: tuple[int, ...], vector: list[int],
                 totals: Counter) -> tuple[int, int]:
    assert any(vector)
    for j in indices:
        assert sum(v * pow(fan.root, i * j, fan.modulus)
                   for i, v in enumerate(vector)) % fan.modulus == 0
    orbit = negacyclic_orbit(vector)
    for column in zip(*orbit):
        for j in indices:
            assert sum(v * pow(fan.root, i * j, fan.modulus)
                       for i, v in enumerate(column)) % fan.modulus == 0
    det = abs(determinant(orbit))
    index = fan.modulus ** len(indices)
    norm_squared = sum(v * v for v in vector)
    assert det and det % index == 0
    assert index ** 2 <= det ** 2 <= norm_squared ** fan.h
    totals["orbit_determinants"] += 1
    return det, norm_squared


def orbit_checks(totals: Counter) -> list[dict]:
    records = []
    for h, bits in ((2, 1), (4, 1), (4, 2), (8, 2)):
        fan = Fanout(h, bits)
        for indices in recipient_sets(h):
            if not indices:
                continue
            basis = kernel_basis(fan, indices)
            vectors = basis + [
                [sum(((-1) ** i) * (i + 1) * basis[i][j] for i in range(h))
                 for j in range(h)]
            ]
            maximum_index = 0
            for vector in vectors:
                det, _ = verify_orbit(fan, indices, vector, totals)
                maximum_index = max(maximum_index, det // (fan.modulus ** len(indices)))
            records.append({
                "h": h,
                "digit_bits": bits,
                "recipients": list(indices),
                "kernel_index": str(fan.modulus ** len(indices)),
                "vectors_checked": len(vectors),
                "largest_orbit_sublattice_index": str(maximum_index),
            })
    return records


def threshold_collision(totals: Counter) -> dict:
    fan = Fanout(8, 128)
    small = 1 << 64
    indices = (1, 7)
    left = (0, 0, small, 0, 0, 0, 0, 0)
    right = (small, fan.base - 1, 0, 0, 0, 0, 0, 0)
    assert left != right
    assert all(0 <= value < fan.base for value in left + right)
    messages = tuple(fan.message_at(left, j) for j in indices)
    assert messages == tuple(fan.message_at(right, j) for j in indices)
    vector = [x - y for x, y in zip(left, right)]
    det, norm_squared = verify_orbit(fan, indices, vector, totals)
    assert det == fan.modulus ** 2
    assert norm_squared == fan.base ** 2 + 1
    totals["actual_source_boundary_collisions"] += 1
    return {
        "h": fan.h,
        "digit_bits": fan.digit_bits,
        "source_parameters": True,
        "recipients": list(indices),
        "first_input": [str(x) for x in left],
        "second_input": [str(x) for x in right],
        "common_messages": [str(x) for x in messages],
        "kernel_vector": [str(x) for x in vector],
        "orbit_determinant_equals_kernel_index": True,
        "norm_squared": str(norm_squared),
    }


def compositions(total: int, length: int):
    if length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, length - 1):
                yield (first, *rest)


def grouping_checks(totals: Counter) -> list[dict]:
    records = []
    for h in (2, 4, 8):
        for count in (2, 3, 4):
            q_vectors = list(compositions(2 * h, count))
            lowest_slack = None
            for odds in compositions(h, count):
                capped = [min(4 * x, h) for x in odds]
                total_capped = sum(capped)
                weights = [total_capped - x for x in capped]
                for source_counts in q_vectors:
                    local = h * max(source_counts)
                    crossing = sum(q * w for q, w in zip(source_counts, weights))
                    slack = local + crossing - 2 * h * h
                    assert slack >= 0
                    lowest_slack = slack if lowest_slack is None else min(lowest_slack, slack)
                    totals["group_count_assignments"] += 1
            records.append({"h": h, "groups": count, "least_integer_slack": lowest_slack})
    for h in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024):
        width = 2 * h * h
        log_upper = h.bit_length()  # ceil(log2(h+1)) for a power of two.
        for count in range(2, h + 1):
            assert 2 * (count - 1) * log_upper <= width
            totals["source_error_budget_controls"] += 1
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    totals = Counter()
    images = []
    for h, bits in ((2, 1), (4, 1), (4, 2), (8, 2), (2, 8)):
        images.extend(image_checks(Fanout(h, bits), totals))
    report = {
        "theorem": "primitive-recipient rank envelope",
        "scope": "local original fan-out and declared independent layer interfaces; no new factoring exponent",
        "images": images,
        "orbit_controls": orbit_checks(totals),
        "source_threshold_collision": threshold_collision(totals),
        "grouping_controls": grouping_checks(totals),
        "summary": dict(sorted(totals.items())),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
