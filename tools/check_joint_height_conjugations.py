#!/usr/bin/env python3
"""Exhaust all binary joint-height conjugations and certify no pressure gain."""

import argparse
from collections import Counter
from hashlib import sha256
from itertools import accumulate
import json
from pathlib import Path

from check_block_smith_crossings import toeplitz_rank
from joint_height_conjugations import (
    binary_product, conjugate, exact_profile, joint_height_centralizers,
    leading_infinity_operator, restored_endpoints, source_matrix,
)


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    source = source_matrix()
    zero, infinity = restored_endpoints(source)
    baseline = exact_profile(zero)
    assert exact_profile(infinity) == baseline
    leading = leading_infinity_operator(source)
    counts, records, matrices = Counter(), {}, set()
    rank_histogram = Counter()
    unchanged = leading_commutants = controls = 0
    for index, flags, matrix, inverse in joint_height_centralizers():
        assert matrix not in matrices
        matrices.add(matrix)
        changed = conjugate(source, matrix, inverse)
        changed_zero, changed_infinity = restored_endpoints(changed)
        assert changed_zero == conjugate(zero, matrix, inverse)
        # All 1536 zero profiles are resolved independently, beyond the unit argument.
        assert exact_profile(changed_zero) == baseline
        profile = exact_profile(changed_infinity)
        assert all(a <= b for a, b in zip(accumulate(profile), accumulate(baseline)))
        zero_rank = profile.count(0)
        rank_histogram[zero_rank] += 1
        same_tensor = changed == source
        commutes_leading = binary_product(matrix, leading) == binary_product(leading, matrix)
        assert (zero_rank == 1) == same_tensor == commutes_leading
        unchanged += same_tensor
        leading_commutants += commutes_leading
        counts[profile] += 1
        if profile not in records:
            for precision in range(1, 4):
                assert toeplitz_rank(changed_infinity, precision) == sum(
                    max(0, precision - order) for order in profile)
                controls += 1
            records[profile] = {
                "center_index": index, "swap_flags": flags,
                "joint_matrix_rows": matrix,
                "infinity_orders": profile,
                "prefix_sum_slacks": [b - a for a, b in zip(
                    accumulate(profile), accumulate(baseline))],
            }
    assert len(matrices) == 1536 and len(records) == 22
    assert unchanged == leading_commutants == 32
    for profile, record in records.items():
        record["multiplicity"] = counts[profile]
    result = {
        "status": "pass",
        "scope": "All invertible binary joint pair maps preserving total height and commuting with swap, for the four-letter reflection seed",
        "joint_conjugations": len(matrices),
        "distinct_infinity_profiles": len(records),
        "zero_profile_unchanged": baseline,
        "independent_toeplitz_controls": controls,
        "unchanged_tensors": unchanged,
        "leading_operator_commutants": leading_commutants,
        "infinity_zero_rank_histogram": dict(sorted(rank_histogram.items())),
        "all_scale_pressure_at_least_original_seed": True,
        "profiles": list(records.values()),
        "engine_sha256": sha256((ROOT / "tools/joint_height_conjugations.py").read_bytes()).hexdigest(),
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "seed_checker_sha256": sha256((ROOT / "tools/check_triangular_auxiliary_row_fiber.py").read_bytes()).hexdigest(),
        "floating_arithmetic_used_for_certification": False,
        "factoring_exponent_improved": False,
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
        print(json.dumps({key: result[key] for key in (
            "status", "joint_conjugations", "distinct_infinity_profiles",
            "independent_toeplitz_controls", "unchanged_tensors")}))
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
