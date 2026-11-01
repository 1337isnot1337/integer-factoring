#!/usr/bin/env python3
"""Exhaust every F_2 link tensor in the stated three-piece boundary family."""
import argparse
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path


def rank4(mask):
    pivots = [0] * 4
    answer = 0
    for i in range(4):
        row = mask >> (4*i) & 15
        while row:
            p = row.bit_length()-1
            if pivots[p]:
                row ^= pivots[p]
            else:
                pivots[p] = row
                answer += 1
                break
    return answer


def rank2(mask):
    a, b = mask & 3, mask >> 2
    return int(bool(a or b)) + int(bool(a and b and a != b))


def radical_sum_at_least(p, q, bound):
    remainder = bound*bound-p-q
    return remainder <= 0 or 4*p*q >= remainder*remainder


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    columns = [0] * 65536
    permutation = [1 << (8*x+4*u+2*y+v)
                   for x, y, u, v in product(range(2), repeat=4)]
    for mask in range(1, 65536):
        bit = mask & -mask
        columns[mask] = columns[mask ^ bit] ^ permutation[bit.bit_length()-1]
    ranks = [rank4(mask) for mask in range(65536)]
    costs = [ranks[mask]*ranks[columns[mask]] for mask in range(65536)]
    expected = {1: [(0, 4), (1, 1), (4, 0)],
                2: [(0, 16), (1, 9), (4, 4), (9, 1), (16, 0)]}
    frontiers = []
    inequalities = 0
    for target_rank in (1, 2):
        target = sum(1 << (8*x+4*y+2*y+x)
                     for x in range(2) for y in range(target_rank))
        pairs = {}
        for tensor in range(65536):
            pair = costs[tensor ^ target], costs[tensor]
            pairs.setdefault(pair, tensor)
            assert radical_sum_at_least(*pair, 2*target_rank)
            inequalities += 1
        frontier = sorted(pair for pair in pairs if not any(
            other != pair and other[0] <= pair[0] and other[1] <= pair[1]
            for other in pairs))
        assert frontier == expected[target_rank]
        for p, q in frontier:
            assert isqrt(p)**2 == p and isqrt(q)**2 == q
            assert isqrt(p)+isqrt(q) == 2*target_rank
        frontiers.append({"target_matrix_rank": target_rank,
                          "target_tensor_mask": target,
                          "frontier": [{"rank_products": pair,
                                        "witness_tensor_mask": pairs[pair]}
                                       for pair in frontier]})
    # Every nontrivial G gives two nonzero terminal matrices, and the
    # exact rank triangle bound needed for the second row of K.
    identity = 9
    terminal_cases = []
    for g in range(1, 16):
        if g == identity:
            continue
        rank_g, rank_complement = rank2(g), rank2(g ^ identity)
        assert rank_g+rank_complement >= 2
        terminal_cases.append([g, rank_g, rank_complement])
    report = {"status": "passed", "field": "F_2", "external_labels": 2,
              "tensors_per_target_rank": 65536,
              "exact_radical_inequalities": inequalities,
              "frontiers": frontiers,
              "all_nontrivial_terminal_matrices": terminal_cases,
              "exact_optimal_spectral_base": 4,
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope": "Specified three-piece right-boundary family and word-rank spectral bound only; no general tensor-rank or factoring lower bound."}
    rendered = json.dumps(report, indent=2, sort_keys=True)+"\n"
    print(rendered, end="")
    if args.output:
        args.output.write_text(rendered)


if __name__ == "__main__":
    main()
