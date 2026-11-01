#!/usr/bin/env python3
"""Replay exact shell ranks, projected entries, and the asymptotic witness."""

import argparse
import hashlib
import json
from pathlib import Path

from weighted_energy_pencil import (
    asymptotic_type_witness, calculate, calculate_tilting,
    central_detector_bound, layer_count, layer_labels,
    primitive_entry, primitive_polynomial, projected_entry,
    rank_two, reduced_entry_polynomial, tilting_blocks,
)


def check_projected_entries(R):
    checked = 0
    for on, off, a, b, shift, _, _ in tilting_blocks(R):
        offset, groups_on, groups_off = 0, [], []
        for k, c in zip(on, off):
            groups_on.append(tuple(range(offset, offset + k)))
            groups_off.append(tuple(range(offset + k, offset + k + c)))
            offset += k + c
        left = layer_labels(tuple(groups_on), a)
        right = layer_labels(tuple(groups_off), b)
        masks_on = tuple(sum(1 << i for i in group) for group in groups_on)
        masks_off = tuple(sum(1 << i for i in group) for group in groups_off)
        dg = {tuple(((x ^ y) & mask).bit_count() for mask in masks_on)
              for k, x in enumerate(left) for y in left[k:]}
        do = {tuple(((x ^ y) & mask).bit_count() for mask in masks_off)
              for k, x in enumerate(right) for y in right[k:]}
        for x in dg:
            for y in do:
                direct = primitive_polynomial(
                    reduced_entry_polynomial(on, off, x, y, shift), R)
                assert direct == projected_entry(a, b, shift, x, y)
                checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true", help="include R8 finite detector")
    parser.add_argument("--output")
    args = parser.parse_args()

    binary_partitions = [0] * 129
    binary_partitions[0] = 1
    for weight in (1, 2, 4, 8, 16, 32, 64, 128):
        for n in range(weight, len(binary_partitions)):
            binary_partitions[n] += binary_partitions[n - weight]
    assert all(layer_count(tuple(range(1, 9)), n) == binary_partitions[n]
               for n in range(129))

    basis_checks = 0
    for b in range(65):
        polynomials = [primitive_entry(0, t + e, b - t)
                       for t in range(b + 1) for e in (0, 1)]
        assert rank_two(polynomials) == 2 * b + 2
        assert max(p.bit_length() for p in polynomials) == 2 * b + 2
        basis_checks += 1

    profiles, entry_checks = [], 0
    for R in range(1, 6):
        direct = calculate(R, primitive=True)
        fast = calculate_tilting(R, primitive=True)
        assert direct["coefficient_ranks"] == fast["coefficient_ranks"]
        entry_checks += check_projected_entries(R)
        profiles.append(fast)
    bounds = [central_detector_bound(R) for R in range(1, 9 if args.full else 8)]
    assert 577 ** 2 > 2 * 408 ** 2
    for bound in bounds:
        assert bound["strict_seed_obstruction_certified"] == (bound["R"] >= 3)

    witnesses = []
    for R in (26, 32, 64, 128, 256, 512):
        witness = asymptotic_type_witness(R)
        witness.pop("multiplicity_lower_hex")
        witness.pop("types")
        witnesses.append(witness)

    result = {
        "status": "PASS",
        "field": "every characteristic-two field",
        "scope": "weighted-energy reciprocal coefficient/symmetry span only",
        "proved_asymptotic_optimal_rho": 1,
        "new_factoring_exponent": None,
        "binary_partition_counts_checked": 129,
        "polynomial_basis_instances_checked": basis_checks,
        "projected_entry_identities_checked": entry_checks,
        "all_coefficient_rank_agreements_through_R": 5,
        "profiles": profiles,
        "finite_central_functional_bounds": bounds,
        "large_energy_structural_witnesses": witnesses,
        "source_sha256": {
            name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ("weighted_energy_pencil.py", "check_weighted_energy.py")
        },
    }
    if args.output:
        Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ("profiles", "finite_central_functional_bounds",
                       "large_energy_structural_witnesses")}, indent=2))


if __name__ == "__main__":
    main()
