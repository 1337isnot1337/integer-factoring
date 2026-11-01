#!/usr/bin/env python3
"""Exhaust the homogeneous middle-fiber retractions and certify one gain."""

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_block_smith_crossings import toeplitz_rank
from reflection_homogeneous_restrictions import (
    WORD_HEIGHTS, middle_fiber_retraction, restricted_endpoint_profiles,
)


ROOT = Path(__file__).resolve().parents[1]


def order_polynomial(orders):
    counts = Counter(orders)
    return [counts[degree] for degree in range(max(counts) + 1)]


def multiply(first, second):
    output = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            output[i + j] += a * b
    return output


def pressure_difference(first, second):
    polys = []
    for record in (first, second):
        polys.append(multiply(order_polynomial(record["zero_orders"]),
                              order_polynomial(record["infinity_orders"])))
    difference = [(polys[0][i] if i < len(polys[0]) else 0)
                  - (polys[1][i] if i < len(polys[1]) else 0)
                  for i in range(max(map(len, polys)))]
    while difference[-1] == 0:
        difference.pop()
    original = list(difference)
    zero_order = 0
    while difference[0] == 0:
        difference.pop(0)
        zero_order += 1
    unit_order = 0
    while sum(difference) == 0:
        divided, running = [], 0
        for coefficient in difference[:-1]:
            running += coefficient
            divided.append(running)
        difference = divided
        unit_order += 1
    assert zero_order == 9 and unit_order == 2 and difference[0] == 1
    degree = len(difference) - 1
    negative_numerator = sum(-coefficient * 2**(degree - j)
                             for j, coefficient in enumerate(difference) if coefficient < 0)
    denominator = 2**degree
    assert negative_numerator < denominator
    # Reconstruct the original pressure difference independently.
    restored = [0] * zero_order + list(difference)
    for _ in range(unit_order):
        restored = multiply(restored, [1, -1])
    assert restored == original
    return {
        "coordinate_minus_nonorthogonal": original,
        "zero_order": zero_order,
        "one_minus_r_order": unit_order,
        "remaining_polynomial": difference,
        "interval": "0 <= r <= 1/2",
        "uniform_positive_quotient_lower_bound": str(Fraction(
            denominator - negative_numerator, denominator)),
        "pressure_strictly_improves_for_every_finite_s_at_least_two": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    pairs = [(vector, functional) for vector, functional in product(range(1, 16), repeat=2)
             if (vector & functional).bit_count() % 2 == 1]
    assert len(pairs) == 120
    if args.small_only:
        pairs = [(1, 1), (1, 15), (15, 1), (15, 2)]
    records, controls = [], 0
    for vector, functional in pairs:
        embedding, retraction, heights = middle_fiber_retraction(vector, functional)
        keep_blocks = (vector, functional) in ((1, 1), (1, 15))
        profile = restricted_endpoint_profiles(
            embedding, retraction, heights, keep_blocks=keep_blocks)
        if keep_blocks:
            for endpoint in profile.pop("blocks"):
                for block in endpoint:
                    for precision in range(1, 5):
                        assert toeplitz_rank(block["matrix"], precision) == sum(
                            max(0, precision - order) for order in block["orders"])
                        controls += 1
        assert Counter(profile["type_energies"]) == Counter(
            [4 * height for height in WORD_HEIGHTS if height != 4] + [16] * 3)
        records.append({"kernel_vector": vector, "image_functional": functional,
                        "embedding": embedding, "retraction": retraction, **profile})
    baseline = next(record for record in records
                    if (record["kernel_vector"], record["image_functional"]) == (1, 1))
    improved = next(record for record in records
                    if (record["kernel_vector"], record["image_functional"]) == (1, 15))
    comparison = pressure_difference(baseline, improved)
    distinct = len({(record["zero_orders"], record["infinity_orders"]) for record in records})
    if not args.small_only:
        assert distinct == 24
    rejected = 0
    for vector, functional in ((0, 1), (1, 0), (1, 2), (16, 1)):
        try:
            middle_fiber_retraction(vector, functional)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("An invalid fiber projector was accepted.")
    embedding = [1 << index for index in range(16)]
    retraction = list(embedding)
    embedding[0] ^= 1 << 1
    retraction[1] ^= 1 << 0
    try:
        restricted_endpoint_profiles(embedding, retraction, WORD_HEIGHTS)
    except ValueError as error:
        assert "height" in str(error)
        rejected += 1
    else:
        raise AssertionError("A nonhomogeneous physical retraction was accepted.")
    result = {
        "status": "pass",
        "scope": "All120 homogeneous rank-three retractions of the two-site reflection middle-height fiber; a strict gain over one coordinate deletion, no gain over the unrestricted seed claimed",
        "complete_fiber_enumeration": not args.small_only,
        "retractions_checked": len(records),
        "distinct_profile_pairs": distinct,
        "independent_toeplitz_rank_controls": controls,
        "malformed_or_nonhomogeneous_rejections": rejected,
        "pointwise_comparison": comparison,
        "records": records,
        "engine_sha256": sha256((ROOT / "tools/reflection_homogeneous_restrictions.py").read_bytes()).hexdigest(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "floating_arithmetic_used_for_certification": False,
        "factoring_exponent_improved": False,
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
        print(json.dumps({key: result[key] for key in (
            "status", "retractions_checked", "distinct_profile_pairs",
            "independent_toeplitz_rank_controls", "malformed_or_nonhomogeneous_rejections")}))
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
