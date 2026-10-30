#!/usr/bin/env python3
"""Verify exact Johnson ranks and complete folded-Gram source factorizations."""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations, product
import argparse
import hashlib
import json
from math import comb, sqrt
from pathlib import Path
import random

from check_complement_shell_crossings import shell
from check_homogeneous_central_symbols import binomial_two
from check_modular_crossing import rank_two
from check_residue_folded_channels import (
    component_dimension, explicit_folded_rows,
)
from complement_shell_tensor import folded_profile, multiply, pair_occupations
from folded_gram_ranks import (
    binomial_parity, central_block, complete_component_rank, complete_gram_channels,
    folded_gram_rank, frobenius_weights, johnson_binomial_rank,
    residue_occupancy_counts, tensor_power_blocks, unweighted_complete_gram_channels,
)


def verify_elementary_split():
    # Polynomials in a,b,c,d use four base-4 exponent digits. Total degree
    # never exceeds three, so ordinary carryless multiplication is exact.
    a, b, c, d = (1 << exponent for exponent in (1, 4, 16, 64))
    g = ((a, b), (c, d))
    determinant = multiply(a, d) ^ multiply(b, c)
    tensor = [[0] * 8 for _ in range(8)]
    for row, col in product(range(8), repeat=2):
        value = 1
        for bit in range(3):
            value = multiply(value, g[row >> bit & 1][col >> bit & 1])
        tensor[row][col] = value
    blocks = [[0] * 8 for _ in range(8)]
    for row, col in product(range(4), repeat=2):
        first = g[row & 1][col & 1]
        second = g[row >> 1 & 1][col >> 1 & 1]
        blocks[row][col] = multiply(first, multiply(second, second))
    for offset in (4, 6):
        for row, col in product(range(2), repeat=2):
            blocks[offset + row][offset + col] = multiply(determinant, g[row][col])
    columns = (
        1 << 0, (1 << 1) | (1 << 2) | (1 << 4),
        (1 << 3) | (1 << 5) | (1 << 6), 1 << 7,
        (1 << 1) | (1 << 2), (1 << 5) | (1 << 6),
        (1 << 2) | (1 << 4), (1 << 3) | (1 << 5),
    )
    assert rank_two(list(columns)) == 8
    checked = 0
    for row, col in product(range(8), repeat=2):
        left = right = 0
        for k in range(8):
            if columns[col] >> k & 1:
                left ^= tensor[row][k]
            if columns[k] >> row & 1:
                right ^= blocks[k][col]
        assert left == right, (row, col)
        checked += 1
    return {"generic_polynomial_intertwiner_entries": checked,
            "constant_change_of_basis_columns": list(columns), "basis_rank": 8}


def verify_characters(maximum=64):
    for power in range(maximum + 1):
        character = [0] * (power + 1)
        for highest, twist, multiplicity in tensor_power_blocks(power):
            local = [1]
            for weight in frobenius_weights(highest):
                following = local + [0] * weight
                for degree, value in enumerate(local):
                    following[degree + weight] += value
                local = following
            for degree, value in enumerate(local):
                character[degree + twist] += multiplicity * value
        assert character == [comb(power, degree) for degree in range(power + 1)]
    return maximum + 1


def dense_johnson_rank(p, offset, order, sign=1):
    masks = [sum(1 << k for k in subset) for subset in combinations(range(2 * p), p)]
    rows = [sum(binomial_parity(offset + sign * (a & b).bit_count(), order) << j
                for j, b in enumerate(masks)) for a in masks]
    return rank_two(rows)


def verify_johnson():
    positive = 0
    for p, i, q in product(range(7), range(4), range(9)):
        expected = dense_johnson_rank(p, q - 1, q - i)
        assert folded_gram_rank(i, p, q) == expected
        assert johnson_binomial_rank(p, q - 1, q - i) == expected
        assert expected <= component_dimension(i, p, q)
        positive += 1
    rng = random.Random(20261030119)
    negative = 0
    for _ in range(160):
        length = rng.randrange(1, 10)
        p = rng.randrange(min(5, length) + 1)
        i = rng.randrange(length - p + 1)
        shift = rng.randrange(-length, length + 1)
        order = length + shift - 2 * i - p
        expected = dense_johnson_rank(p, length - i, order, -1)
        assert complete_component_rank(length, i, p, shift) == expected
        negative += 1
    return {"unmerged_dense_Johnson_ranks": positive,
            "merged_dense_Johnson_ranks": negative}


def verify_expanded_grams():
    checked = 0
    for modulus, length in ((4, 1), (8, 2), (8, 3)):
        for energy in range(0, 41, 2):
            vectors, _ = explicit_folded_rows(modulus, length, energy)
            gram = [sum(((a & b).bit_count() & 1) << j for j, b in enumerate(vectors))
                    for a in vectors]
            expected = 0
            for i, p, residue_energy, count in residue_occupancy_counts(modulus, length):
                winding = energy - residue_energy
                if winding >= 0 and winding % modulus == 0:
                    expected += count * folded_gram_rank(i, p, winding // modulus)
            assert rank_two(gram) == expected
            checked += 1
    return checked


def complement_mask(mask, modes):
    return sum(((mask >> i) & 1) << (modes - i - 1) for i in range(modes))


def incoming_generators(parameter, state, modes, complemented):
    """Equality-forced incoming coordinates, independent of winding."""
    common = state & (parameter if complemented else complement_mask(parameter, modes))
    annihilated = common
    result = defaultdict(list)
    while True:
        unchanged = state ^ annihilated
        created = (complement_mask(parameter ^ annihilated, modes) if complemented
                   else parameter ^ complement_mask(annihilated, modes))
        double = unchanged & created
        single = unchanged ^ created
        result[double, single].append(unchanged & ~created)
        if not annihilated:
            break
        annihilated = (annihilated - 1) & common
    return result


def merged_gram_polynomial(first, second, length, reverse=False):
    result = 0
    for (double, single), rows in first.items():
        i = double.bit_count()
        k = i + single.bit_count() // 2
        for a in rows:
            for b in second.get((double, single), ()):
                result ^= binomial_two(length - i - (a & b).bit_count()) << (k + i)
    if reverse:
        result = sum(((result >> degree) & 1) << (2 * length - degree)
                     for degree in range(2 * length + 1))
    return result


def verify_complete_entries():
    records = []
    for modulus, length, weight in ((8, 2, 8), (16, 3, 23), (16, 4, 32)):
        labels = shell(modulus, length, weight)
        occupations = [pair_occupations(parts, modulus) for parts in labels]
        masks = [sum(1 << ((part - 1) // 2) for part in parts) for parts in labels]
        d = len(labels)
        first = {(x, y): incoming_generators(masks[x], masks[y], modulus // 2, False)
                 for x, y in product(range(d), repeat=2)}
        second = {(x, y): incoming_generators(masks[x], masks[y], modulus // 2, True)
                  for x, y in product(range(d), repeat=2)}
        for x, y, u, v in product(range(d), repeat=4):
            expected = folded_profile(occupations[x], occupations[y], occupations[u], occupations[v])
            row = merged_gram_polynomial(first[x, y], first[v, u], length)
            column = merged_gram_polynomial(second[x, u], second[v, y], length, True)
            assert row == column == expected
        records.append({"mode_cutoff": modulus, "length": length, "weight": weight,
                        "labels": d, "complete_entry_polynomials_in_both_cuts": d ** 4})
    return records


def shell_dimension(modulus, length, weight):
    counts = [[0] * (weight + 1) for _ in range(length + 1)]
    counts[0][0] = 1
    for part in range(1, modulus, 2):
        for size in range(length, 0, -1):
            for energy in range(weight, part - 1, -1):
                counts[size][energy] += counts[size - 1][energy - part]
    return counts[length][weight]


def summarize_channels(modulus, length, weight, actual=None):
    dimension = shell_dimension(modulus, length, weight)
    profile = complete_gram_channels(modulus, length, weight)
    if actual is not None:
        channels = {shift: (r, c) for shift, r, c in profile}
        for shift, r, c in actual:
            first, second = channels[shift // modulus]
            assert r <= first and c <= second
    capped = [[shift * modulus, min(dimension ** 2, r), min(dimension ** 2, c)]
              for shift, r, c in profile if r and c]
    return {
        "mode_cutoff": modulus, "length": length, "operator_weight": weight,
        "partner_weight": length * modulus - weight, "labels": dimension,
        "uncapped_merged_Gram_channels": [[shift * modulus, r, c] for shift, r, c in profile],
        "capped_merged_Gram_channels": capped,
        "rank_sum_upper_bound": sum(sqrt(r * c) for _, r, c in capped),
        "finite_seed_threshold": dimension ** sqrt(2),
        "actual_complete_ranks_checked": actual is not None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {
        "status": "pass",
        "scope": "exact polynomial-time Johnson ranks and complete two-cut Gram channels; no improved factoring exponent",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "rank_module_sha256": hashlib.sha256(Path(__file__).with_name("folded_gram_ranks.py").read_bytes()).hexdigest(),
        "elementary_GL2_split": verify_elementary_split(),
        "tensor_character_checks": verify_characters(),
        "Johnson_rank_checks": verify_johnson(),
        "explicit_high_mode_Gram_checks": verify_expanded_grams(),
        "complete_entry_checks": verify_complete_entries(),
    }
    print(json.dumps({key: report[key] for key in (
        "tensor_character_checks", "Johnson_rank_checks", "explicit_high_mode_Gram_checks")}), flush=True)
    source = Path(__file__).resolve().parents[1] / "evidence/complement-shell-ranks.json"
    profiles = json.loads(source.read_text())["cases"]
    report["complete_profile_evidence_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    report["complete_profile_bounds"] = [
        summarize_channels(case["mode_cutoff"], case["length"], case["operator_weight"],
                           case["nonzero_profile"]) for case in profiles]
    report["dense_length_examples"] = [
        summarize_channels(modulus, length, modulus * length // 2)
        for modulus, length in ((32, 8), (64, 16))]
    report["large_Johnson_examples"] = [
        {"half_universe": p, "matrix_dimension": comb(2 * p, p),
         "largest_explicit_block": max(len(central_block(m)[0])
                                       for m, _, _ in tensor_power_blocks(2 * p)),
         "ranks_at_winding_powers_of_two": [[q, folded_gram_rank(0, p, q)]
                                           for q in (1, 2, 4, 8, 16, 32, 64) if q <= p]}
        for p in (8, 16, 32, 64)]
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"dense_length_examples": report["dense_length_examples"]}), flush=True)


if __name__ == "__main__":
    main()
