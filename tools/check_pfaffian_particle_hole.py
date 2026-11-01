#!/usr/bin/env python3
"""Verify the generic particle-hole/exterior identity and exact cut ranks."""

from functools import lru_cache
from itertools import combinations, product
import argparse
import hashlib
import json
from pathlib import Path

from check_complement_shell_crossings import shell
from complement_shell_tensor import multiply, pair_occupations
from explore_complement_edge_pencils import (
    all_length_labels, compile_code, evaluate_local,
    rank_profile as parity_rank_profile,
)
from pfaffian_particle_hole import (
    CHARGE, PARTICLE_HOLE, local_exterior_table, one_particle_matrix,
    particle_hole_mask, rank_profile,
)


def minor(matrix, first, second):
    if len(first) != len(second):
        return 0
    @lru_cache(None)
    def recurse(position, columns):
        if position == len(first):
            return 1
        result = 0
        candidates = columns
        while candidates:
            bit = candidates & -candidates
            j = bit.bit_length() - 1
            result ^= multiply(matrix[first[position]][second[j]],
                               recurse(position + 1, columns ^ bit))
            candidates ^= bit
        return result
    return recurse(0, (1 << len(second)) - 1)


def verify_generic_identity():
    # Base-16 monomial exponents: every individual exponent is at most eight,
    # so this is an injective encoding of all multivariate identities below.
    p, q, r = 1 << 1, 1 << 16, 1 << 256
    table = tuple(evaluate_local((p, q, r, 1, 1, 1)))
    assert table == local_exterior_table(p, q, r)
    g = one_particle_matrix(p, q, r)
    delta = p ^ r
    powers = [1]
    for _ in range(4):
        powers.append(multiply(powers[-1], delta))
    checked = unequal = charge = column = 0
    reflected = local_exterior_table(q, p, r)
    flip = (0, 2, 1, 3)
    for x, y, u, v in product(range(4), repeat=4):
        value = table[((x * 4 + y) * 4 + u) * 4 + v]
        first, second = particle_hole_mask(x, y), particle_hole_mask(v, u)
        a = [i for i in range(4) if first >> i & 1]
        b = [i for i in range(4) if second >> i & 1]
        if len(a) != len(b):
            assert value == 0
            unequal += 1
        elif not a:
            assert value == delta
        else:
            assert multiply(value, powers[len(a) - 1]) == minor(g, a, b)
        checked += 1
        if CHARGE[x] + CHARGE[y] != CHARGE[v] + CHARGE[u]:
            assert value == 0
            charge += 1
        assert value == reflected[((x * 4 + flip[u]) * 4 + flip[y]) * 4 + v]
        column += 1
    return {
        "generic_multivariate_local_table_entries": 256,
        "exterior_minor_identities": checked,
        "vanishing_unequal_exterior_degrees": unequal,
        "signed_charge_zero_entries": charge,
        "column_cut_parameter_swap_identities": column,
        "particle_hole_permutation": list(PARTICLE_HOLE),
    }


def verify_grade_and_energy_maps():
    checked = 0
    modulus = 8
    lows = (1, 3)
    for x, y in product(product(range(4), repeat=2), repeat=2):
        grade_x = sum(mask.bit_count() for mask in x)
        grade_y = sum(mask.bit_count() for mask in y)
        energy_x = sum(low * (a & 1) + (modulus - low) * (a >> 1)
                       for low, a in zip(lows, x))
        energy_y = sum(low * (a & 1) + (modulus - low) * (a >> 1)
                       for low, a in zip(lows, y))
        masks = [particle_hole_mask(a, b) for a, b in zip(x, y)]
        counts = [sum((mask >> i) & 1 for mask in masks) for i in range(4)]
        assert grade_x == 2 + counts[2] - counts[1]
        assert grade_y == 2 + counts[3] - counts[0]
        assert energy_x == sum(lows) + sum(
            (modulus - low) * ((mask >> 2) & 1) - low * ((mask >> 1) & 1)
            for low, mask in zip(lows, masks))
        assert energy_y == sum(lows) + sum(
            (modulus - low) * ((mask >> 3) & 1) - low * (mask & 1)
            for low, mask in zip(lows, masks))
        checked += 1
    return checked


def verify_profiles(include_large=False):
    records = []
    for pairs, length in ((2, 2), (3, 3), (4, 2), (4, 4)):
        labels = all_length_labels(pairs, length)
        result = rank_profile(labels)
        # A separate parity-block implementation evaluates Pfaffian matchings.
        reference = parity_rank_profile(compile_code(labels), (6, 3, 2, 1, 1, 1))
        assert result["profile"] == reference["profile"]
        result.update({"complement_pairs": pairs, "length": length,
                       "independent_parity_block_comparison": True})
        records.append(result)
    labels = [pair_occupations(label, 16) for label in shell(16, 3, 23)]
    result = rank_profile(labels)
    reference = parity_rank_profile(compile_code(labels), (6, 3, 2, 1, 1, 1))
    assert result["profile"] == reference["profile"]
    result.update({"complement_pairs": 4, "length": 3, "energy": 23,
                   "independent_parity_block_comparison": True})
    records.append(result)
    if include_large:
        result = rank_profile(all_length_labels(5, 5))
        result.update({"complement_pairs": 5, "length": 5,
                       "independent_parity_block_comparison": False})
        records.append(result)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--include-large", action="store_true")
    args = parser.parse_args()
    report = {
        "status": "pass",
        "scope": "generic exact exterior-power reduction and actual complete cut ranks; no improved factoring exponent",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "module_sha256": hashlib.sha256(Path(__file__).with_name("pfaffian_particle_hole.py").read_bytes()).hexdigest(),
        "generic_identity": verify_generic_identity(),
        "grade_and_energy_map_checks": verify_grade_and_energy_maps(),
    }
    print(json.dumps({"generic_identity": report["generic_identity"],
                      "grade_and_energy_map_checks": report["grade_and_energy_map_checks"]}), flush=True)
    report["complete_profiles"] = verify_profiles(args.include_large)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report["complete_profiles"]), flush=True)


if __name__ == "__main__":
    main()
