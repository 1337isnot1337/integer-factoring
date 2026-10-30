#!/usr/bin/env python3
"""Exact dimensions and factors for residue-folded high-mode channels."""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
import argparse
import hashlib
import json
from math import comb
from pathlib import Path
import random

from check_modular_crossing import determinant_operator_two, rank_two, regular_partitions
from check_partition_channels import operator


def component_dimension(doubles, half_singles, winding):
    """Rank of all folded vectors with one fixed residue occupancy pattern."""
    i, p, q = doubles, half_singles, winding
    if min(i, p, q) < 0 or q < i:
        return 0
    if i:
        return comb(2 * p, min(p, q - i))
    if not q:
        return 1
    if not p:
        return 0
    return comb(2 * p, min(p, q)) - int(q < (p & -p))


@lru_cache(None)
def folded_dimensions(modulus, length, maximum_energy):
    """All exact ambient folded dimensions up to the stated energy."""
    if modulus < 2 or modulus & (modulus - 1) or min(length, maximum_energy) < 0:
        raise ValueError("a power-of-two modulus and nonnegative bounds are required")
    # A residue can be unused, doubled, or single.  The single count is even
    # in a valid component.  Its total particle count is at most 2*length.
    counts = {(0, 0, 0): 1}
    for residue in range(1, modulus, 2):
        following = defaultdict(int)
        for (doubles, singles, energy), count in counts.items():
            following[doubles, singles, energy] += count
            if 2 * doubles + singles + 2 <= 2 * length and energy + 2 * residue <= maximum_energy:
                following[doubles + 1, singles, energy + 2 * residue] += count
            if 2 * doubles + singles + 1 <= 2 * length and energy + residue <= maximum_energy:
                following[doubles, singles + 1, energy + residue] += count
        counts = following
    result = [0] * (maximum_energy + 1)
    for (doubles, singles, energy), count in counts.items():
        if singles & 1:
            continue
        for winding in range((maximum_energy - energy) // modulus + 1):
            result[energy + winding * modulus] += count * component_dimension(
                doubles, singles // 2, winding)
    return tuple(result)


def component_rows(doubles, half_singles, winding):
    """Independent incidence-matrix realization of one residue component."""
    i, p, q = doubles, half_singles, winding
    if q < i:
        return []
    universe = tuple(range(2 * p))
    fixed = [sum(1 << x for x in subset) for subset in combinations(universe, p)]
    if i:
        sizes = range(min(p, q - i) + 1)
    elif not q:
        sizes = (0,)
    else:
        sizes = range(1, min(p, q) + 1)
    positive = [sum(1 << x for x in subset)
                for size in sizes for subset in combinations(universe, size)]
    return [sum(int(not (a & t)) << k for k, t in enumerate(positive)) for a in fixed]


def compositions(total, length):
    if not length:
        if not total:
            yield ()
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first, *tail)


def folded_states(modulus, unchanged, residues, winding):
    """Explicit finite expansion used only for small independent checks."""
    unchanged = tuple(unchanged)
    for heights in compositions(winding, len(residues)):
        created = tuple(r + modulus * height for r, height in zip(residues, heights))
        state = tuple(sorted(unchanged + created))
        if len(set(state)) == len(state):
            yield state


def folded_generators(modulus, length, energy):
    low = tuple(range(1, modulus, 2))
    for size in range(min(length, len(low)) + 1):
        for unchanged, residues in product(combinations(low, size), repeat=2):
            remaining = energy - sum(unchanged) - sum(residues)
            if remaining >= 0 and remaining % modulus == 0:
                yield unchanged, residues, remaining // modulus


def explicit_folded_rows(modulus, length, energy):
    states = regular_partitions(2, energy)
    indices = {state: i for i, state in enumerate(states)}
    result = []
    for unchanged, residues, winding in folded_generators(modulus, length, energy):
        value = 0
        for state in folded_states(modulus, unchanged, residues, winding):
            value ^= 1 << indices[state]
        result.append(value)
    return result, states


def positive_composition_parity(total, parts):
    if not parts:
        return int(total == 0)
    if total < parts:
        return 0
    # Lucas parity of binomial(total-1, parts-1).
    return int((parts - 1) & ~(total - 1) == 0)


def folded_output_coefficient(modulus, unchanged, residues, winding, parameter, output):
    """Evaluate the outgoing factor without expanding high-mode assignments."""
    fixed = set(unchanged)
    doubled = set(residues) & fixed
    variable_singles = tuple(sorted(set(residues) - fixed))
    answer = 0
    for mask in range(1 << len(variable_singles)):
        positive = doubled | {r for k, r in enumerate(variable_singles) if mask >> k & 1}
        count = len(positive)
        if not positive_composition_parity(winding, count):
            continue
        extra = winding - count
        distinguished = min(positive) if positive else None
        created = tuple(r + modulus * (1 + (extra if r == distinguished else 0))
                        if r in positive else r for r in residues)
        state = tuple(sorted(tuple(unchanged) + created))
        assert len(set(state)) == len(state)
        answer ^= determinant_operator_two(parameter, output, state)
    return answer


def verify():
    component_checks = 0
    for i, p, q in product(range(4), range(7), range(10)):
        assert rank_two(component_rows(i, p, q)) == component_dimension(i, p, q), (i, p, q)
        component_checks += 1
    explicit_checks = 0
    for modulus, length, maximum in ((4, 1, 24), (8, 2, 40), (8, 3, 40), (16, 3, 32)):
        dimensions = folded_dimensions(modulus, length, maximum)
        energies = range(maximum + 1) if modulus <= 8 else (0, 8, 16, 24, 32)
        for energy in energies:
            rows, _ = explicit_folded_rows(modulus, length, energy)
            assert rank_two(rows) == dimensions[energy], (modulus, length, energy)
            explicit_checks += 1
    # Every incoming operator vector is checked against the whole folded span.
    image_checks = 0
    modulus, length, maximum = 8, 2, 32
    low_labels = tuple(combinations(range(1, modulus, 2), length))
    for energy in range(0, maximum + 1, 2):
        rows, states = explicit_folded_rows(modulus, length, energy)
        expected = rank_two(rows)
        vectors = [sum(determinant_operator_two(parameter, state, input_) << k
                       for k, state in enumerate(states))
                   for parameter, input_ in product(low_labels, repeat=2)]
        assert rank_two(rows + vectors) == expected
        image_checks += len(vectors)
    # Compare the finite positive-composition formula with direct expansion.
    rng = random.Random(202610307)
    factor_checks = original_checks = 0
    candidates = [(energy, item) for energy in range(0, 49, 2)
                  for item in folded_generators(8, 2, energy)]
    for _ in range(300):
        _, (unchanged, residues, winding) = rng.choice(candidates)
        parameter, output = rng.choice(low_labels), rng.choice(low_labels)
        expected = 0
        original = 0
        for state in folded_states(8, unchanged, residues, winding):
            expected ^= determinant_operator_two(parameter, output, state)
            if factor_checks < 40:
                original ^= operator(parameter, output, state) & 1
        actual = folded_output_coefficient(8, unchanged, residues, winding, parameter, output)
        assert actual == expected
        if factor_checks < 40:
            assert original == expected
            original_checks += 1
        factor_checks += 1
    stability_checks = 0
    for modulus, length in ((4, 1), (8, 2), (8, 3), (16, 4)):
        threshold = 3 * length * modulus - 2 * length * length
        dimensions = folded_dimensions(modulus, length, threshold + 3 * modulus)
        for energy in range(threshold, threshold + 2 * modulus + 1):
            assert dimensions[energy] == dimensions[energy + modulus]
            stability_checks += 1
        assert sum(dimensions[threshold:threshold + modulus]) == sum(
            comb(modulus // 2, size) ** 2 for size in range(1, length + 1))
    return {
        "component_incidence_ranks": component_checks,
        "explicit_high_mode_span_ranks": explicit_checks,
        "complete_incoming_operator_vectors": image_checks,
        "outgoing_factor_expansion_comparisons": factor_checks,
        "original_integer_operator_factor_comparisons": original_checks,
        "eventual_periodicity_checks": stability_checks,
        "stable_residue_total_checks": 4,
    }


def partition_counts(maximum):
    result = [0] * (maximum + 1)
    result[0] = 1
    for part in range(1, maximum + 1, 2):
        for energy in range(maximum, part - 1, -1):
            result[energy] += result[energy - part]
    return result


def channel_case(modulus, length, weight, profile):
    partner_weight = length * modulus - weight
    total = weight + partner_weight
    angles = [index for index in range(total + 1) if index & ~total == 0]
    energies = [energy for shift, _, _ in profile for angle in angles
                for energy in (weight + shift - partner_weight + angle,
                               weight + partner_weight - shift - angle)]
    maximum = max(0, *energies)
    dimensions = folded_dimensions(modulus, length, maximum)
    partitions = partition_counts(maximum)
    assert all(a <= b for a, b in zip(dimensions, partitions))

    def size(table, energy):
        return table[energy] if energy >= 0 else 0

    pieces = []
    for shift, r, c in profile:
        row_energies = [weight + shift - partner_weight + angle for angle in angles]
        column_energies = [weight + partner_weight - shift - angle for angle in angles]
        row = sum(size(dimensions, energy) for energy in row_energies)
        column = sum(size(dimensions, energy) for energy in column_energies)
        assert r <= row and c <= column, (modulus, length, shift, r, c, row, column)
        pieces.append({
            "shift": shift, "actual_ranks": [r, c], "folded_channels": [row, column],
            "partition_channels": [sum(size(partitions, energy) for energy in row_energies),
                                   sum(size(partitions, energy) for energy in column_energies)],
        })
    return {"mode_cutoff": modulus, "length": length, "state_energy": weight,
            "operator_weight": weight, "partner_weight": partner_weight,
            "angle_indices": angles, "pieces": pieces}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "pass",
        "scope": "exact ambient high-mode folded-space dimensions and valid two-cut channels; no improved factoring exponent",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "verification": verify(), "dimension_examples": [],
    }
    for modulus, length in ((8, 2), (16, 3), (16, 4), (32, 3), (32, 4), (32, 8), (64, 16)):
        energy = modulus * length
        count = folded_dimensions(modulus, length, energy)[energy]
        report["dimension_examples"].append({
            "mode_cutoff": modulus, "length_bound": length, "energy": energy,
            "folded_dimension": count, "regular_partition_dimension": partition_counts(energy)[energy],
        })
    evidence_path = Path(__file__).resolve().parents[1] / "evidence/complement-shell-ranks.json"
    if evidence_path.exists():
        profiles = json.loads(evidence_path.read_text())["cases"]
        report["complete_piece_channel_checks"] = [
            channel_case(case["mode_cutoff"], case["length"], case["operator_weight"],
                         case["nonzero_profile"]) for case in profiles]
        report["complete_profile_evidence_sha256"] = hashlib.sha256(evidence_path.read_bytes()).hexdigest()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"verification": report["verification"],
                      "dimension_examples": report["dimension_examples"]}), flush=True)


if __name__ == "__main__":
    main()
