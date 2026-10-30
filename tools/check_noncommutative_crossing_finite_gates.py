#!/usr/bin/env python3
"""Exact universal, orbit-support, and constant-augmentation rank-sum gates."""

from itertools import product
import argparse
import hashlib
import json
import random
from pathlib import Path

from check_noncommutative_crossing_group_degeneration import (
    BilinearGroup, word_polynomial, both_ranks, combine,
)
from check_modular_crossing import rank_two


class Cyclic:
    def __init__(self, order):
        self.d = order
    def mul(self, a, b):
        return (a + b) % self.d
    def inverse(self, a):
        return -a % self.d


class Dihedral:
    def __init__(self, rotations):
        self.rotations = rotations
        self.d = 2 * rotations
    def mul(self, a, b):
        x, s = a % self.rotations, a // self.rotations
        y, t = b % self.rotations, b // self.rotations
        return (x + (-1 if s else 1) * y) % self.rotations + self.rotations * (s ^ t)
    def inverse(self, a):
        x, s = a % self.rotations, a // self.rotations
        return ((1 if s else -1) * x) % self.rotations + self.rotations * s


def support_rank_one_audit(group, name, audit_units=True):
    d = group.d
    abelian = all(group.mul(a, b) == group.mul(b, a)
                  for a, b in product(range(d), repeat=2))
    elementary = all(group.mul(a, a) == 0 for a in range(d))
    support_solutions, orbit_solutions = 0, 0
    # Every NONZERO rank-one-column orbit tensor is forced to one of
    # these d pairs of permutations by the proof's y=e equation.
    for constant in range(d):
        f = [group.mul(group.inverse(x), constant) for x in range(d)]
        g = [group.mul(group.inverse(y), constant) for y in range(d)]
        supported = all(
            g[y] == group.mul(group.mul(x, group.inverse(y)), f[x])
            for x, y in product(range(d), repeat=2))
        if not supported:
            continue
        support_solutions += 1
        consistent = True
        # A nonzero kernel entry must retain the support for every free y.
        for S in range(d):
            a = group.mul(group.inverse(S), constant)
            for y in range(d):
                x, u = group.mul(S, y), group.mul(y, a)
                v = group.mul(S, u)
                if u != f[x] or v != g[y]:
                    consistent = False
        orbit_solutions += consistent
    assert support_solutions == (d if abelian else 0)
    assert orbit_solutions == (d if elementary else 0)
    units = 0
    for kernel in range(1 << d) if audit_units else ():
        if not kernel.bit_count() & 1:
            continue
        rows = []
        for y in range(d):
            row = 0
            for u in range(d):
                a = group.mul(group.inverse(y), u)
                if kernel >> a & 1:
                    row |= 1 << u
            rows.append(row)
        assert rank_two(rows) == d
        units += 1
    return {"group": name, "order": d, "abelian": abelian,
            "elementary_abelian": elementary,
            "rank_one_support_candidates": support_solutions,
            "rank_one_orbit_candidates": orbit_solutions,
            "augmentation_one_convolution_units_checked": units}


def arbitrary_control_matrices(group, theta, weights):
    d = group.d
    kernels = [word_polynomial(group, theta[group.word_index[S]], weights,
                               tuple(range(group.width)))
               for S in range(d)]
    assert all(__import__("functools").reduce(int.__xor__, kernel.values(), 0) == 1
               for kernel in kernels)
    matrices = [[0] * (d * d) for _ in range(sum(weights) + 1)]
    for S, y in product(range(d), repeat=2):
        x = group.mul(S, y)
        for a, polynomial in kernels[S].items():
            u = group.mul(y, a)
            v = group.mul(S, u)
            while polynomial:
                bit = polynomial & -polynomial
                matrices[bit.bit_length() - 1][x * d + y] ^= 1 << (u * d + v)
                polynomial ^= bit
    endpoint = combine(matrices, (1 << len(matrices)) - 1)
    assert endpoint == [1 << (y * d + x) for x, y in product(range(d), repeat=2)]
    return matrices


def arbitrary_control_audit(group, name):
    rng = random.Random(20261030)
    maps = [list(range(8)), [0] * 8, [7] * 8,
            [i & 3 for i in range(8)], [(i * i) & 7 for i in range(8)]]
    maps += [[rng.randrange(8) for _ in range(8)] for _ in range(11)]
    count = 0
    minimum = 64
    for weights in ((1, 1, 1), (1, 1, 2), (1, 2, 4)):
        for theta in maps:
            matrices = arbitrary_control_matrices(group, theta, weights)
            for mask in range(1, 1 << len(matrices), 2):
                row, col = both_ranks(combine(matrices, mask), 3)
                assert row == 64 and col >= 2
                minimum = min(minimum, col)
                count += 1
    return {"group": name, "arbitrary_control_maps": len(maps),
            "weight_assignments": 3, "unit_recombinations_checked": count,
            "minimum_column_rank_seen": minimum,
            "all_unit_row_ranks": 64}


def exact_integer_gates():
    assert 99 * 99 > 2 * 70 * 70
    small = []
    for d in range(2, 6):
        left = (1 << 140) * (d * d - 1) ** 70
        right = d ** 198
        assert left > right
        small.append({"labels": d, "lower_power": str(left),
                      "target_upper_power": str(right), "strictly_excluded": True})
    non_elementary = []
    for d, first_num, second_num, den, target in (
        (6, 8, 5, 1, 13), (7, 48, 32, 5, 16), (8, 56, 39, 5, 19)
    ):
        assert 2 * (d * d - 1) * den * den > first_num * first_num
        assert (d * d - 2) * den * den > second_num * second_num
        assert first_num + second_num == target * den
        assert target ** 70 > d ** 99
        non_elementary.append({
            "labels": d, "rank_sum_strictly_above": target,
            "target_strictly_below": target,
            "first_radical_lower": [first_num, den],
            "second_radical_lower": [second_num, den],
            "lower_integer_power": str(target ** 70),
            "target_upper_integer_power": str(d ** 99),
        })
    assert 2 * 25 > 7 * 7
    assert 62 * 25 > 39 * 39
    assert 19 ** 70 > 8 ** 99
    return {"sqrt_two_upper": [99, 70],
            "universal_d2_through_d5": small,
            "non_elementary_orbit_d6_through_d8": non_elementary,
            "order_eight": {"rank_sum_strictly_above": 19,
                           "target_strictly_below": 19,
                           "nineteen_power_seventy": str(19 ** 70),
                           "eight_power_ninety_nine": str(8 ** 99)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    d8 = BilinearGroup(2, 1, ((0, 1), (0, 0)))
    q8 = BilinearGroup(2, 1, ((1, 1), (0, 1)))
    c4c2 = BilinearGroup(2, 1, ((1, 0), (0, 0)))
    elementary = BilinearGroup(3, 0, ((0, 0, 0), (0, 0, 0), (0, 0, 0)))
    records = [support_rank_one_audit(group, name) for name, group in (
        ("D8", d8), ("Q8", q8), ("C4_times_C2", c4c2),
        ("C8", Cyclic(8)), ("C2_cubed", elementary))]
    records += [support_rank_one_audit(group, name, audit_units=False)
                for name, group in (("S3", Dihedral(3)), ("D10", Dihedral(5)),
                                    ("C3", Cyclic(3)), ("C9", Cyclic(9)))]
    result = {
        "status": "Exact integer gates, convolution units, and orbit-support audits passed.",
        "scope": "Universal d<=5 and non-elementary-group orbit d<=8 gates over every field; stronger constant-augmentation 2-group gate in characteristic two.",
        "integer_certificates": exact_integer_gates(),
        "rank_one_orbit_and_unit_audits": records,
        "arbitrary_control_audits": [arbitrary_control_audit(d8, "D8"),
                                    arbitrary_control_audit(q8, "Q8")],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
