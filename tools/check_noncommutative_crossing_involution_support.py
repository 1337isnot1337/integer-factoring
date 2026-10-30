#!/usr/bin/env python3
"""Audit the involution-free diagonal minor of every group orbit kernel."""

from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import random

from check_noncommutative_crossing_abelian_midpoints import rank_mod


class TableGroup:
    def __init__(self, labels, multiply):
        self.labels = tuple(labels)
        self.d = len(self.labels)
        index = {label: i for i, label in enumerate(self.labels)}
        self.table = [[index[multiply(a, b)] for b in self.labels] for a in self.labels]
        self.inv = [next(b for b in range(self.d)
                         if self.table[a][b] == self.table[b][a] == 0)
                    for a in range(self.d)]
        assert all(self.table[0][a] == self.table[a][0] == a for a in range(self.d))
        assert all(self.table[self.table[a][b]][c] == self.table[a][self.table[b][c]]
                   for a, b, c in product(range(self.d), repeat=3))
    def mul(self, a, b):
        return self.table[a][b]


def cyclic(order):
    return TableGroup(range(order), lambda a, b: (a + b) % order)


def dihedral(rotations):
    return TableGroup(product(range(2), range(rotations)),
                      lambda a, b: (a[0] ^ b[0],
                                    (a[1] + (-1 if a[0] else 1) * b[1]) % rotations))


def quaternion(half_rotation_order):
    n = half_rotation_order
    return TableGroup(product(range(2), range(2 * n)),
                      lambda a, b: (a[0] ^ b[0],
                                    (a[1] + (-1 if a[0] else 1) * b[1] + n * a[0] * b[0]) % (2 * n)))


def direct_product(first, second):
    return TableGroup(product(range(first.d), range(second.d)),
                      lambda a, b: (first.mul(a[0], b[0]), second.mul(a[1], b[1])))


def maximum_involution_free(group):
    d = group.d
    involutions = [a for a in range(1, d) if group.mul(a, a) == 0]
    neighbors = []
    for a in range(d):
        neighbors.append(sum(1 << b for b in range(d)
                             if a != b and group.mul(group.mul(a, group.inv[b]),
                                                     group.mul(a, group.inv[b])) == 0))
    best = 0
    def recurse(candidates, selected):
        nonlocal best
        if selected.bit_count() + candidates.bit_count() <= best.bit_count():
            return
        if not candidates:
            best = selected
            return
        choices = [i for i in range(d) if candidates >> i & 1]
        vertex = max(choices, key=lambda i: (neighbors[i] & candidates).bit_count())
        rest = candidates & ~(1 << vertex)
        recurse(rest & ~neighbors[vertex], selected | (1 << vertex))
        recurse(rest, selected)
    recurse((1 << d) - 1, 0)
    subset = [a for a in range(d) if best >> a & 1]
    assert all(a == b or group.mul(group.mul(a, group.inv[b]),
                                   group.mul(a, group.inv[b])) != 0
               for a, b in product(subset, repeat=2))
    return subset, involutions


def entry(group, kernel, x, y, u, v):
    S = group.mul(x, group.inv[y])
    if S != group.mul(v, group.inv[u]):
        return 0
    return kernel[S][group.mul(group.inv[y], u)]


def dense_column(group, kernel):
    d = group.d
    result = [[0] * (d * d) for _ in range(d * d)]
    for S, y, a in product(range(d), repeat=3):
        x = group.mul(S, y)
        u = group.mul(y, a)
        v = group.mul(S, u)
        result[x * d + u][v * d + y] = kernel[S][a]
    return result


def audit_kernel(group, subset, kernel, prime):
    d = group.d
    S0, a0 = next((S, a) for S, a in product(range(d), repeat=2) if kernel[S][a])
    normalized = [[kernel[group.mul(S0, S)][group.mul(a, a0)]
                   for a in range(d)] for S in range(d)]
    value = normalized[0][0]
    assert value
    for x, y in product(subset, repeat=2):
        assert entry(group, normalized, x, y, x, y) == (value if x == y else 0)
    # Recheck the local relabelling on every supported normalized position.
    for S, y, a in product(range(d), repeat=3):
        x = group.mul(S, y)
        u = group.mul(y, a)
        v = group.mul(S, u)
        assert entry(group, normalized, x, y, u, v) == entry(
            group, kernel, group.mul(S0, x), y, group.mul(u, a0),
            group.mul(group.mul(S0, v), a0))
    original_rank = rank_mod(dense_column(group, kernel), prime)
    normalized_rank = rank_mod(dense_column(group, normalized), prime)
    assert original_rank == normalized_rank >= len(subset)
    return original_rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    elementary = TableGroup(range(8), lambda a, b: a ^ b)
    odd_nonabelian = TableGroup(product(range(3), range(7)),
                               lambda a, b: ((a[0] + b[0]) % 3,
                                             (a[1] + pow(2, a[0], 7) * b[1]) % 7))
    examples = (
        ("C4", cyclic(4), 2), ("C9", cyclic(9), 9),
        ("C2_cubed", elementary, 1), ("D8", dihedral(4), 2),
        ("Q8", quaternion(2), 4), ("Q16", quaternion(4), 8),
        ("D10", dihedral(5), 5), ("D12", dihedral(6), 3),
        ("S3_times_C3", direct_product(dihedral(3), cyclic(3)), 9),
        ("nonabelian_order_21", odd_nonabelian, 21),
    )
    rng = random.Random(20261030)
    reports = []
    for name, group, expected in examples:
        subset, involutions = maximum_involution_free(group)
        assert len(subset) == expected
        ranks = []
        for prime in (2, 3):
            for _ in range(2):
                kernel = [[rng.randrange(prime) for _ in range(group.d)] for _ in range(group.d)]
                ranks.append({"field": prime, "column_rank": audit_kernel(group, subset, kernel, prime)})
        reports.append({"group": name, "order": group.d,
                        "nonidentity_involutions": involutions,
                        "maximum_involution_free_size": len(subset),
                        "maximum_subset": subset, "kernel_rank_checks": ranks})
    assert 3 ** 18 > 2 ** 5 * 5 ** 10
    assert 17 ** 2 > 2 * 12 ** 2
    result = {
        "status": "All group laws, optimal subsets, normalized diagonal minors, and dense ranks passed.",
        "scope": "Arbitrary finite-group orbit kernels over every field; finite arithmetic audits over F2 and F3.",
        "examples": reports,
        "rank_floor_target_certificate": [3 ** 18, 2 ** 5 * 5 ** 10],
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
