#!/usr/bin/env python3
"""Exact finite checks of the uniform eighth-jet identity over F_2.

The all-size statement is proved by the finite unmatched-field analysis in
homogeneous-central-eighth-jet.md; these computations check its formulas.
"""

from __future__ import annotations

import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path
import random

from check_homogeneous_central_profiles import python_profile, original_profile
from check_homogeneous_central_symbols import binomial_two
from check_modular_crossing import determinant_operator_two, rank_two, regular_partitions


CROSS_PAIRS = frozenset(((1, 1), (1, 3), (1, 5), (1, 7), (3, 5)))
SAME_PAIRS = frozenset(((1, 7), (3, 5)))


def vacuum_pairing(a, b):
    if len(a) != len(b):
        return 0
    rows = [sum(int(n & m == 1) << i for i, m in enumerate(b)) for n in a]
    return int(rank_two(rows) == len(a))


def gamma_eight(a, b):
    """Sum of full-contraction minors after the permitted unmatched pairs."""
    vertices = [(0, n) for n in a] + [(1, m) for m in b]
    result = 0
    for i, j in combinations(range(len(vertices)), 2):
        side, n = vertices[i]
        other, m = vertices[j]
        permitted = SAME_PAIRS if side == other else CROSS_PAIRS
        if tuple(sorted((n, m))) not in permitted:
            continue
        remaining_a = tuple(value for k, (s, value) in enumerate(vertices)
                            if s == 0 and k != i and k != j)
        remaining_b = tuple(value for k, (s, value) in enumerate(vertices)
                            if s == 1 and k != i and k != j)
        result ^= vacuum_pairing(remaining_a, remaining_b)
    return result


def charge_eight(state):
    return sum(part % 8 in (1, 7) for part in state) & 1


def translate_low(polynomial, bound=8):
    result = 0
    mask = (1 << (bound + 1)) - 1
    while polynomial:
        bit = polynomial & -polynomial
        result ^= binomial_two(bit.bit_length() - 1) & mask
        polynomial ^= bit
    return result


def residual_low_jet(a, u, b, y):
    h, j, energy = sum(b), sum(a), sum(y)
    assert j <= h and sum(u) == energy
    polynomial = python_profile(a, u, b, y)
    scalar = vacuum_pairing(a, b) * int(u == y)
    assert polynomial.bit_count() & 1 == scalar
    residual = polynomial ^ (scalar << (energy + j))
    return translate_low(residual), residual


def generalized_binomial_two(exponent, degree):
    if exponent >= 0:
        return int(degree & ~exponent == 0)
    return math.comb(-exponent + degree - 1, degree) & 1


def quadratic_zero_mode_jet(n, m, state, same_slot=False, bound=8):
    """Independent two-field normal-order formula on a fixed state."""
    a, b = (n - 1) // 2, (m - 1) // 2
    result = 0
    for part in state:
        c = (part - 1) // 2
        ca, aa = int(a & ~c == 0), int(a & c == 0)
        cb, ab = int(b & ~c == 0), int(b & c == 0)
        if same_slot:
            result ^= ca * ab ^ aa * cb
        else:
            for degree in range(bound + 1):
                coefficient = (ca * ab * generalized_binomial_two(part, degree)
                               ^ aa * cb * generalized_binomial_two(-part, degree))
                result ^= coefficient << degree
    return result


def verify():
    pair_checks = 0
    # Each table entry gives the order left after its regularizing power.
    table = {(1, 1): 6, (1, 3): 4, (3, 1): 4,
             (1, 5): 2, (5, 1): 2, (3, 3): None,
             (1, 7): 0, (7, 1): 0, (3, 5): 0, (5, 3): 0}
    for energy in range(0, 34, 2):
        for state in regular_partitions(2, energy):
            for (n, m), order in table.items():
                jet = quadratic_zero_mode_jet(n, m, state)
                if order is None:
                    assert jet == 0
                else:
                    assert jet & ((1 << order) - 1) == 0
                    assert (jet >> order) & 1 == charge_eight(state)
                pair_checks += 1
            for n, m in ((1, 3), (1, 5), (1, 7), (3, 5)):
                value = quadratic_zero_mode_jet(n, m, state, same_slot=True)
                assert value == (charge_eight(state) if (n, m) in SAME_PAIRS else 0)
                pair_checks += 1

    entry_checks = 0
    original_checks = 0
    for h in range(1, 13):
        for j in range(h + 1):
            if (h - j) & 1:
                continue
            for a, b in product(regular_partitions(2, j), regular_partitions(2, h)):
                assert vacuum_pairing(a, b) == determinant_operator_two(a, (), b)
                for energy in (0, 4, 8, 10, 12):
                    for u, y in product(regular_partitions(2, energy), repeat=2):
                        jet, _ = residual_low_jet(a, u, b, y)
                        assert jet == ((gamma_eight(a, b) * int(u == y) * charge_eight(y)) << 8)
                        entry_checks += 1
                        if h <= 4 and energy <= 4:
                            assert python_profile(a, u, b, y) == original_profile(a, u, b, y)
                            original_checks += 1

    rng = random.Random(820261030)
    large_checks = 0
    for _ in range(200):
        j, h = rng.choice(((8, 24), (16, 16), (24, 40), (30, 34), (32, 32), (40, 40)))
        energy = rng.choice((16, 24, 32, 40))
        a = rng.choice(regular_partitions(2, j))
        b = rng.choice(regular_partitions(2, h))
        y = rng.choice(regular_partitions(2, energy))
        u = y if rng.randrange(2) else rng.choice(regular_partitions(2, energy))
        jet, _ = residual_low_jet(a, u, b, y)
        assert jet == ((gamma_eight(a, b) * int(u == y) * charge_eight(y)) << 8)
        large_checks += 1

    records = []
    for h in (8, 16, 24, 32, 40):
        labels = regular_partitions(2, h)
        rows = [sum(gamma_eight(a, b) << i for i, b in enumerate(labels)) for a in labels]
        charge_rank = sum(charge_eight(y) for y in labels)
        records.append({"weight": h, "labels": len(labels),
                        "gamma_eight_rank": rank_two(rows), "state_charge_rank": charge_rank,
                        "complete_T8_rank_in_either_cut": rank_two(rows) * charge_rank})

    a, b, state = (1, 3, 5, 23), (1, 3, 13, 15), (1, 5, 9, 17)
    jet, residual = residual_low_jet(a, state, b, state)
    assert jet == 1 << 8
    full_translation = translate_low(residual, residual.bit_length())
    assert (full_translation & -full_translation).bit_length() - 1 == 8
    odd_jet, odd_residual = residual_low_jet((1,), (1,), (1,), (1,))
    assert odd_jet & 255 != 0
    return {"status": "pass", "scope": "all-size eighth-jet formula supported by exact finite checks; no improved factoring exponent",
            "quadratic_zero_mode_checks": pair_checks, "complete_small_entry_checks": entry_checks,
            "original_integer_operator_checks": original_checks, "sampled_larger_entry_checks": large_checks,
            "complete_T8_rank_records": records,
            "order_eight_witness": {"operator": list(a), "partner": list(b), "state": list(state),
                                    "gamma_eight": 1, "charge_eight": 1, "residual_order": 8},
            "odd_energy_negative_control": {"operator": [1], "partner": [1], "state": [1],
                                            "low_residual_jet_bits": odd_jet},
            "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = verify()
    encoded = json.dumps(report, sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")
