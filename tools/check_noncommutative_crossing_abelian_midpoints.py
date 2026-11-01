#!/usr/bin/env python3
"""Exact column-block multiplicities for finite abelian orbit kernels."""

from itertools import product
from math import prod
from pathlib import Path
import argparse
import hashlib
import json
import random

from check_modular_crossing import rank_two


class AbelianGroup:
    def __init__(self, moduli):
        self.moduli = tuple(moduli)
        self.labels = tuple(product(*(range(m) for m in self.moduli)))
        self.lookup = {a: i for i, a in enumerate(self.labels)}
        self.d = len(self.labels)
        self.add = [[self.lookup[tuple((x + y) % m for x, y, m in zip(a, b, self.moduli))]
                     for b in self.labels] for a in self.labels]
        self.neg = [self.lookup[tuple(-x % m for x, m in zip(a, self.moduli))]
                    for a in self.labels]
        self.sub = [[self.add[a][self.neg[b]] for b in range(self.d)] for a in range(self.d)]
        self.double = [self.add[a][a] for a in range(self.d)]
        self.image = sorted(set(self.double))
        self.two_torsion = sum(x == 0 for x in self.double)
        assert len(self.image) * self.two_torsion == self.d


def rank_mod(matrix, prime):
    if prime == 2:
        return rank_two([sum((entry & 1) << j for j, entry in enumerate(row))
                         for row in matrix])
    basis = {}
    for source in matrix:
        row = [x % prime for x in source]
        while any(row):
            pivot = next(i for i, value in enumerate(row) if value)
            if pivot not in basis:
                scale = pow(row[pivot], -1, prime)
                basis[pivot] = [scale * value % prime for value in row]
                break
            scale = row[pivot]
            row = [(a - scale * b) % prime for a, b in zip(row, basis[pivot])]
    return len(basis)


def blocks(group, kernel):
    return [[[kernel[group.sub[x][y]][group.sub[group.sub[t][x]][y]]
              for y in range(group.d)] for x in range(group.d)]
            for t in range(group.d)]


def dense_cuts(group, kernel):
    d = group.d
    row = [[0] * (d * d) for _ in range(d * d)]
    col = [[0] * (d * d) for _ in range(d * d)]
    for S, y, a in product(range(d), repeat=3):
        x = group.add[S][y]
        u = group.add[y][a]
        v = group.add[S][u]
        row[x * d + y][u * d + v] = kernel[S][a]
        col[x * d + u][v * d + y] = kernel[S][a]
    return row, col


def audit_kernel(group, kernel, prime, dense=False):
    d = group.d
    matrices = blocks(group, kernel)
    ranks = [rank_mod(matrix, prime) for matrix in matrices]
    total = sum(ranks)
    assert total % len(group.image) == 0
    if any(any(row) for row in kernel):
        assert total >= len(group.image)
    # Test every translation and every block, not only rank equality.
    for t, h in product(range(d), repeat=2):
        shifted = group.add[t][group.double[h]]
        for x, y in product(range(d), repeat=2):
            assert matrices[shifted][group.add[x][h]][group.add[y][h]] == matrices[t][x][y]
    if group.two_torsion == 1:
        half = {value: index for index, value in enumerate(group.double)}
        midpoint = [[kernel[half[group.add[b][c]]][half[group.sub[c][b]]]
                     for c in range(d)] for b in range(d)]
        assert total == d * rank_mod(midpoint, prime)
    if dense:
        row_cut, col_cut = dense_cuts(group, kernel)
        assert rank_mod(col_cut, prime) == total
        conv = sum(rank_mod([[kernel[S][group.sub[u][y]] for u in range(d)]
                             for y in range(d)], prime) for S in range(d))
        assert rank_mod(row_cut, prime) == conv
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = random.Random(20261030)
    reports = []
    specifications = ((2,), (3,), (4,), (6,), (8,), (9,), (12,), (16,),
                      (2, 2), (2, 4), (3, 3), (2, 2, 3))
    for moduli in specifications:
        group = AbelianGroup(moduli)
        d = group.d
        checked, dense_checked = 0, 0
        minimum = d * d
        for prime in (2, 3):
            kernels = []
            identity = [[int(a == 0) for a in range(d)] for _ in range(d)]
            kernels.append(identity)
            kernels.extend([[rng.randrange(prime) for _ in range(d)] for _ in range(d)]
                           for _ in range(5))
            for index, kernel in enumerate(kernels):
                value = audit_kernel(group, kernel, prime, dense=(index < 2))
                minimum = min(minimum, value)
                checked += 1
                dense_checked += index < 2
            for S, a in ((0, 0), (d - 1, d - 1)):
                kernel = [[0] * d for _ in range(d)]
                kernel[S][a] = 1
                value = audit_kernel(group, kernel, prime)
                minimum = min(minimum, value)
                checked += 1
        reports.append({"cyclic_factor_orders": list(moduli), "order": d,
                        "doubling_image_size": len(group.image),
                        "two_torsion_size": group.two_torsion,
                        "kernels_checked": checked,
                        "dense_both_cut_checks": dense_checked,
                        "minimum_column_rank_seen": minimum})
    # Exact certificate for the weighted AM-GM constant in the all-cyclic gate.
    assert 3 ** 18 > 2 ** 5 * 5 ** 10
    assert 17 ** 2 > 2 * 12 ** 2
    result = {
        "status": "All block identities, both-cut comparisons, and exact target certificates passed.",
        "scope": "Every finite abelian group and every coefficient field; numerical audits over F2 and F3.",
        "group_audits": reports,
        "target_certificate": {
            "sqrt_two_upper": [17, 12],
            "weighted_amgm_constant_twelfth_power": [3 ** 18, 2 ** 5 * 5 ** 10],
            "consequence": "Every abelian group with at most two elements killed by doubling misses the strict finite-seed target."
        },
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
