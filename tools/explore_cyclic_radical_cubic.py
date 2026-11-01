#!/usr/bin/env python3
"""Exact all-size cut profiles for cyclic-output cubic symplectic codes."""

import argparse
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import isqrt, sqrt
from pathlib import Path

from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local, minimum_decomposition
from check_modular_crossing import rank_two


LOCAL = evaluate_local((2, 16, 1, 1, 1, 1))


@lru_cache(None)
def canonical_entry(pattern, parity):
    x, y, u, v, fx, fu = pattern
    size = 1 + max(pattern)

    def term(a, b):
        sites = {}
        for leg, single, double in ((0, x, fx), (1, y, a), (2, u, fu), (3, v, b)):
            assert single != double
            sites.setdefault(single, [0, 0, 0, 0])[leg] = 1
            sites.setdefault(double, [0, 0, 0, 0])[leg] = 3
        value = 1
        for aa, bb, cc, dd in sites.values():
            value = multiply(value, LOCAL[((aa * 4 + bb) * 4 + cc) * 4 + dd])
        return value

    value = 0
    for a in range(size):
        if a == y:
            continue
        for b in range(size):
            if b != v:
                value ^= term(a, b)
    # Outside positions are indistinguishable. Distinct ordered outside
    # pairs occur n(n-1) times and cancel; equal pairs occur n times.
    if (parity - size) % 2:
        for a in range(size):
            if a != y:
                value ^= term(a, size)
        for b in range(size):
            if b != v:
                value ^= term(size, b)
        value ^= term(size, size)
    return value


def entry(pairs, x, y, u, v):
    raw = (x, y, u, v, (x + 1) % pairs, (u + 1) % pairs)
    names = {}
    pattern = tuple(names.setdefault(value, len(names)) for value in raw)
    return canonical_entry(pattern, pairs % 2)


def polynomial_gcd(first, second):
    while second:
        while first and first.bit_length() >= second.bit_length():
            first ^= second << (first.bit_length() - second.bit_length())
        first, second = second, first
    return first


def circulant_rank(first_row, size):
    return size - (polynomial_gcd(first_row, (1 << size) | 1).bit_length() - 1)


def block_polynomials(pairs, distance, cut):
    indices = ((0, distance), (distance, 0))
    if cut == 0:
        return [[entry(pairs, x, y, u, v) for v, u in indices] for x, y in indices]
    return [[entry(pairs, x, y, u, v) for v, y in indices] for x, u in indices]


def binary_profile(pairs):
    assert pairs >= 3 and pairs % 2
    diagonal = [[entry(pairs, 0, 0, j, j) for j in range(pairs)],
                [entry(pairs, 0, j, 0, j) for j in range(pairs)]]
    kinds = [(pairs, [block_polynomials(pairs, 1, cut) for cut in range(2)])]
    if pairs > 3:
        kinds.append((pairs * (pairs - 3) // 2,
                      [block_polynomials(pairs, 2, cut) for cut in range(2)]))
    all_values = diagonal[0] + [value for count, cuts in kinds for row in cuts[0] for value in row]
    active = [degree for degree in range(13) if any(value >> degree & 1 for value in all_values)]
    signatures = [sum(((value >> degree) & 1) << i for i, value in enumerate(all_values)) for degree in active]
    pivots, chosen = {}, []
    for degree, signature in zip(active, signatures):
        value = signature
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                chosen.append(degree)
                break
    assert chosen[0] == 0 and 3 not in active and 12 not in active
    ranks = []
    for mask in range(1 << len(chosen)):
        functional = sum(1 << degree for i, degree in enumerate(chosen) if mask >> i & 1)
        pair = []
        for cut in range(2):
            first = sum(((value & functional).bit_count() % 2) << j for j, value in enumerate(diagonal[cut]))
            rank = circulant_rank(first, pairs)
            for count, blocks in kinds:
                rows = [sum(((value & functional).bit_count() % 2) << j for j, value in enumerate(row))
                        for row in blocks[cut]]
                rank += count * rank_two(rows)
            pair.append(rank)
        ranks.append(pair)
    assert ranks[1] == [pairs * pairs] * 2
    cost, witness = minimum_decomposition([sqrt(r * c) for r, c in ranks], target=1)
    lower, _ = minimum_decomposition([isqrt(r * c) for r, c in ranks], target=1)
    return {'pairs': pairs, 'labels': pairs, 'length': 3,
            'basis_degrees': chosen, 'basis_monomials': [[a % 4, a // 4, 3 - a % 4 - a // 4] for a in chosen],
            'all_binary_ranks': ranks, 'cost': cost, 'witness': witness,
            'integer_lower': lower, 'target': pairs ** sqrt(2),
            'exclusion': lower ** 70 > pairs ** 99,
            'diagonal_polynomial_first_rows': diagonal,
            'off_diagonal_polynomial_blocks': kinds}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs', default='3,5,7,9,15,31,63,127,255,511,1023')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = {'scope': 'Exact binary coefficient spaces for cyclic-output cubic radical codes; no claimed exponent improvement',
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'cases': []}
    for pairs in map(int, args.pairs.split(',')):
        case = binary_profile(pairs)
        report['cases'].append(case)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
        print(json.dumps({k: case[k] for k in ('pairs', 'basis_degrees', 'cost', 'integer_lower', 'target', 'exclusion')}), flush=True)


if __name__ == '__main__':
    main()
