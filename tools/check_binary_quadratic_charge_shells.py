#!/usr/bin/env python3
"""Exact checks for quadratic charges and paired homogeneous vacuum labels.

All-size proofs and precise limitations are in binary-quadratic-charge-shells.md.
No factoring exponent is inferred from the finite data.
"""

from collections import defaultdict
from itertools import combinations
import argparse
import hashlib
import json
from math import comb
from pathlib import Path

from check_modular_crossing import determinant_operator_two, modular_operator, regular_partitions


def coefficient(a: int, b: int, c: int) -> int:
    return int(a & ~c == 0 and b & c == 0) ^ int(b & ~c == 0 and a & c == 0)


def direct_coefficient(a: int, b: int, c: int) -> int:
    return (comb(c, a) * comb(c + b, b) + comb(c, b) * comb(c + a, a)) % 2


def syndrome(parts: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    values = [0] * (modulus // 4)
    for part in parts:
        residue = part % modulus
        representative = min(residue, modulus - residue)
        values[(representative - 1) // 2] ^= 1
    return tuple(values)


def balanced_pairs(parts: tuple[int, ...], modulus: int) -> tuple[tuple[int, int], ...]:
    assert len(set(parts)) == len(parts)
    assert not any(syndrome(parts, modulus))
    groups = defaultdict(lambda: [[], []])
    for part in parts:
        residue = part % modulus
        representative = min(residue, modulus - residue)
        groups[representative][int(residue != representative)].append(part)
    result = []
    for first, second in groups.values():
        first.sort()
        second.sort()
        assert len(first) % 2 == len(second) % 2
        if len(first) % 2:
            result.append((first.pop(), second.pop()))
        for values in (first, second):
            result.extend(zip(values[::2], values[1::2]))
    assert len(result) * 2 == len(parts)
    assert all(sum(pair) >= modulus for pair in result)
    return tuple(result)


def paired_labels(modulus: int, chosen: int) -> tuple[tuple[int, ...], ...]:
    pairs = [(k, modulus - k) for k in range(1, modulus // 2, 2)]
    return tuple(tuple(sorted(part for i in selection for part in pairs[i]))
                 for selection in combinations(range(len(pairs)), chosen))


def verify(max_weight: int = 64) -> dict:
    if not 16 <= max_weight <= 80:
        raise ValueError("use a maximum state weight from 16 through 80")
    coefficient_cases = 0
    for a in range(32):
        for b in range(32):
            for c in range(65):
                actual = coefficient(a, b, c)
                assert actual == direct_coefficient(a, b, c)
                if a & b:
                    assert not actual
                elif a != b:
                    mask = a | b
                    assert actual == int(c & mask in (a, b))
                else:
                    assert not actual
                coefficient_cases += 1

    zero_mode_entries = 0
    for weight in range(17):
        states = regular_partitions(2, weight)
        for a in range(8):
            for b in range(a + 1, 8):
                parameters = (2 * a + 1, 2 * b + 1)
                for input_ in states:
                    expected = sum(coefficient(a, b, (part - 1) // 2) for part in input_) % 2
                    for output in states:
                        assert modular_operator(2, parameters, output, input_) == (
                            expected if output == input_ else 0)
                        zero_mode_entries += 1

    residue_function_cases = 0
    for modulus in (4, 8, 16, 32, 64):
        full_mask = modulus // 2 - 1
        # These complementary full-mask pairs are exactly the syndrome coordinates.
        for representative in range(1, modulus // 2, 2):
            a = (representative - 1) // 2
            b = full_mask ^ a
            for c in range(modulus // 2):
                residue = 2 * c + 1
                assert coefficient(a, b, c) == int(residue in (representative, modulus - representative))
                residue_function_cases += 1
        # Every lower-weight charge is invariant under this full complement.
        for a in range(full_mask + 1):
            for b in range(a + 1, full_mask + 1):
                if 2 * (a + b + 1) <= modulus:
                    for c in range(full_mask + 1):
                        assert coefficient(a, b, c) == coefficient(a, b, full_mask ^ c)
                        residue_function_cases += 1

    shell_records = []
    paired_difference_checks = 0
    for modulus in (4, 8, 16, 32, 64):
        for weight in range(max_weight + 1):
            states = regular_partitions(2, weight)
            fibers = defaultdict(list)
            for state in states:
                fibers[syndrome(state, modulus)].append(state)
            for signature, fiber in fibers.items():
                reference = set(fiber[0])
                for state in fiber:
                    difference = tuple(sorted(reference ^ set(state)))
                    pairs = balanced_pairs(difference, modulus)
                    assert sum(difference) <= 2 * weight
                    assert len(difference) <= 4 * weight // modulus
                    assert len(pairs) <= 2 * weight // modulus
                    paired_difference_checks += 1
                if not any(signature):
                    for state in fiber:
                        balanced_pairs(state, modulus)
                        assert len(state) <= 2 * weight // modulus
            modes = (weight + 1) // 2
            bound = sum(comb(modes, j) for j in range(min(modes, 4 * weight // modulus) + 1))
            largest = max(map(len, fibers.values()), default=0)
            assert largest <= bound
            if weight in (16, 32, 48, 64) and weight <= max_weight:
                shell_records.append({
                    "M": modulus, "weight": weight, "state_dimension": len(states),
                    "nonempty_character_spaces": len(fibers),
                    "largest_character_dimension": largest,
                    "zero_character_dimension": len(fibers.get((0,) * (modulus // 4), ())),
                    "proved_Hamming_ball_bound": bound,
                })

    gram_records = []
    for modulus, chosen in ((8, 1), (16, 1), (16, 2), (32, 2), (32, 4)):
        labels = paired_labels(modulus, chosen)
        assert len(labels) == comb(modulus // 4, chosen)
        assert all(sum(label) == modulus * chosen and len(label) == 2 * chosen for label in labels)
        assert all(not any(syndrome(label, modulus)) for label in labels)
        for i, first in enumerate(labels):
            for j, second in enumerate(labels):
                assert determinant_operator_two(first, (), second) == int(i == j)
        gram_records.append({"M": modulus, "R": chosen, "h": modulus * chosen,
                             "dimension": len(labels), "exact_vacuum_entries": len(labels) ** 2})

    witness = (1, 5, 9, 17)
    assert sum(witness) == 32 and syndrome(witness, 8) == (1, 1)
    assert all(not any(syndrome(state, 8)) for state in regular_partitions(2, 16))
    # All quadratic charges on a paired shell can vanish while a quartic charge moves states.
    assert syndrome((1, 7), 8) == syndrome((3, 5), 8) == (0, 0)
    assert modular_operator(2, (1, 3, 5, 7), (3, 5), (1, 7)) == 1
    # The entropy conclusion is specific to all low charges, not one parity test.
    assert len({syndrome(state, 8) for state in regular_partitions(2, 32)}) == 2
    return {
        "scope": "exact checks for an all-size quadratic-charge theorem and paired homogeneous vacuum duality",
        "binary_Lucas_coefficient_cases": coefficient_cases,
        "independent_normal_operator_zero_mode_entries": zero_mode_entries,
        "residue_function_cases": residue_function_cases,
        "paired_symmetric_difference_checks": paired_difference_checks,
        "joint_charge_shells": shell_records,
        "paired_vacuum_families": gram_records,
        "weight_32_nonzero_weight_8_charge_witness": list(witness),
        "quartic_counterexample": {"operator": [1, 3, 5, 7], "input": [1, 7], "output": [3, 5], "coefficient": 1},
        "new_factoring_exponent_claimed": False,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-weight", type=int, default=64)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(args.max_weight)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")
