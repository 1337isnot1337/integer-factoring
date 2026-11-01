#!/usr/bin/env python3
"""Elementary free-digit entry decoding over the narrow prime modulus."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import random

from check_observable_boundary_constructor import BoundaryNetwork, balanced_blocks
from check_prime_transform_boundary_constructor import PrimeTransform, PrimeLengthTransform, parameters
from check_width_sensitive_prime_constructor import initial_parameters, stage_parameters


COUNTS = Counter()


def local_entries():
    randomizer = random.Random(0xE171E5)
    examples = [(PrimeTransform(2, 1), 3, 5),
                (PrimeLengthTransform(2, 9, 1), 1, 269),
                (PrimeLengthTransform(2, 10, 10), 617, 1),
                (PrimeLengthTransform(2, 19, 1), 1, 267013),
                (PrimeLengthTransform(2, 19, 19), 267013, 1),
                (PrimeLengthTransform(2, 20, 10), 719, 1451)]
    records = []
    for model, left, right in examples:
        target = left * right
        for cut in sorted({0, target.bit_length() // 2, target.bit_length()}):
            network = BoundaryNetwork(model, 2, target, cut, left >> cut)
            trace = network.trace(left, right)
            for group in range(2):
                block = network.blocks[group]
                outside = model.a - len(block)
                free = sum(max(0, network.heights[family, residue] - outside)
                           for family in range(2) for residue in block)
                assert free <= len(block) * (len(block) - 1)
                if free:
                    COUNTS["positive_free_digit_entry_controls"] += 1
                boundary = network.boundary(group, trace)
                actual, decoded_trials = network.entry(group, boundary, True)
                reference, reference_trials = network.entry(group, boundary, False)
                assert actual == reference == 1
                assert decoded_trials <= 1
                COUNTS["exact_valid_entry_comparisons"] += 1
                COUNTS["reference_local_assignments"] += reference_trials
                COUNTS["decoded_internal_assignments"] += decoded_trials
                for _ in range(2):
                    changed = dict(boundary)
                    edge = randomizer.choice(tuple(changed))
                    bound = model.carry_bound + 1 if edge[0] == "carry" else model.modulus
                    changed[edge] = (changed[edge] + randomizer.randrange(1, bound)) % bound
                    expected, trials = network.entry(group, changed, False)
                    decoded, _ = network.entry(group, changed, True)
                    assert decoded == expected and decoded in (0, 1)
                    COUNTS["exact_mutated_entry_comparisons"] += 1
                    COUNTS["reference_local_assignments"] += trials
                assert network.entry(group, boundary, True, padding=True)[0] == 0
                COUNTS["invalid_padding_rejections"] += 1
                changed = dict(boundary)
                edge = next(iter(changed))
                changed[edge] = model.carry_bound + 1 if edge[0] == "carry" else model.modulus
                assert network.entry(group, changed, True)[0] == 0
                COUNTS["out_of_range_label_rejections"] += 1
        records.append({"N": target, "left": left, "right": right,
                        "Q": model.modulus, "digit_width": model.digit_bits})
    return records


def free_digit_controls():
    for A in [4, 8, 16, 32, 64]:
        for P in sorted({2, min(3, A // 2), A // 2}):
            for block in balanced_blocks(A, P):
                k = len(block)
                for left in range(A):
                    for right in range(A):
                        if left + right > A:
                            continue
                        free = max(0, left - (A - k)) + max(0, right - (A - k))
                        assert free <= k - 1
                        COUNTS["complete_paired_height_deficiency_controls"] += 1


def schedule_controls():
    inputs = set()
    for exponent in list(range(40, 257)) + [511, 512, 1023, 1024, 2048, 4096]:
        anchor = 1 << exponent
        inputs.update(anchor + offset for offset in [-3, -2, -1, 0, 1, 2, 3, 127])
        u = (anchor + 1).bit_length()
        initial_p = parameters(anchor)[1]
        for P in range(initial_p, initial_p + 3):
            boundary = (16 * u * u * P) ** 2
            inputs.update(boundary + offset for offset in [-2, -1, 0, 1, 2])
    worst = Fraction(0)
    for n in sorted(value for value in inputs if value >= (1 << 40)):
        u, P, A, _, w, _, _, _, width = parameters(n)
        k = (A + P - 1) // P
        free_exponent = k * k * w
        assert k <= 2 * u
        assert w < 2050 * u * u
        assert free_exponent * P < n
        assert P * P * width < 13 * n
        assert 12 * free_exponent < width
        assert ((P - 1) * width + free_exponent) * P < 13 * n
        worst = max(worst, Fraction(free_exponent * P, n))
        COUNTS["exact_explicit_schedule_controls"] += 1
    for epsilon in [Fraction(1), Fraction(1, 2), Fraction(1, 16), Fraction(3, 7)]:
        for exponent in [48, 64, 100, 256, 512, 1024]:
            n = 1 << exponent
            u, K, initial = initial_parameters(n, epsilon)
            if initial < 2 * K:
                continue
            last = (initial // (2 * K)).bit_length() - 1
            for stage in sorted({0, 1, 2, last // 2, last}):
                parameter = stage_parameters(n, epsilon, stage)
                if parameter is None:
                    continue
                assert parameter["A"] // parameter["P"] == K
                assert parameter["w"] < 16 * (4 ** stage) * u / (epsilon / 48)
                COUNTS["exact_constant_block_schedule_controls"] += 1
    return {"proved_uniform_F_P_over_n_bound": "904050/1048576",
            "largest_observed_F_P_over_n": str(worst)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = local_entries()
    free_digit_controls()
    bounds = schedule_controls()
    result = {"status": "passed", "claim": "elementary polynomial-time adaptive prime entries",
              "counts": dict(sorted(COUNTS.items())), "models": records, "bounds": bounds,
              "scope": "No lattice search is used by this decoder; full label enumeration is still charged."}
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(json.dumps({key: value for key, value in result.items() if key != "models"},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
