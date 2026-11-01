#!/usr/bin/env python3
"""Independent polynomial checks for syndrome-preserving cubic output maps."""

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from explore_circulant_cubic_outputs import profile
from pfaffian_particle_hole import local_exterior_table


H, K, A, B, C = 771, 582, 531, 257, 325
TABLE = local_exterior_table(2, 16, 1)


def literal(pairs, x, y, u, v, a, b):
    """Contract every input term literally, using the exterior-minor table."""
    result = 0
    for c in range(pairs):
        if c == y:
            continue
        for e in range(pairs):
            if e == v:
                continue
            value = 1
            for site in range(pairs):
                occupations = [int(site == single) + 3 * int(site == double)
                               for single, double in ((x, a), (y, c), (u, b), (v, e))]
                aa, bb, cc, dd = occupations
                value = multiply(value, TABLE[((aa * 4 + bb) * 4 + cc) * 4 + dd])
                if not value:
                    break
            result ^= value
    return result


def polynomial_checks():
    records = []
    for pairs in (3, 5, 7, 9):
        count = 0
        for a in range(1, pairs):
            for b in range(1, pairs):
                assert literal(pairs, 0, 0, 0, 0, a, b) == H
                value = literal(pairs, 0, 1, 0, 1, a, b)
                assert value == literal(pairs, 0, 1, 0, 1, b, a)
                if a == b:
                    assert value == K
                count += 2
            for b in range(pairs):
                if b == 1:
                    continue
                assert literal(pairs, 0, 0, 1, 1, a, b) == 0
                special_a, special_b = a == 1, b == 0
                expected = C if (special_a and special_b) or a == b else A if special_a or special_b else B
                assert literal(pairs, 0, 1, 1, 0, a, b) == expected
                assert literal(pairs, 1, 0, 0, 1, b, a) == expected
                count += 3
        records.append({'pairs': pairs, 'generic_polynomial_comparisons': count})
    return records


def all_ranks(basis):
    working = [0] * len(basis[0])
    ranks = [0] * (1 << len(basis))
    previous = 0
    for step in range(1, len(ranks)):
        gray = step ^ (step >> 1)
        change = (gray ^ previous).bit_length() - 1
        working = [a ^ b for a, b in zip(working, basis[change])]
        ranks[gray] = rank_two(working)
        previous = gray
    return ranks


def dense_checks():
    records = []
    for pairs, mask in ((3, 1), (5, 1), (5, 7), (7, 7), (9, 31)):
        report = profile(pairs, mask)
        weights = [[int(i != j and (mask >> ((j - i) % pairs - 1)) & 1)
                    for j in range(pairs)] for i in range(pairs)]
        assert all(sum(row) % 2 == 1 for row in weights)
        crossing = {}
        for i in range(pairs):
            for j in range(i + 1, pairs):
                first = weights[i][j] ^ weights[j][i]
                second = (weights[i][j] * weights[j][i]) ^ (sum(a * b for a, b in zip(weights[i], weights[j])) % 2)
                crossing[i, j] = B ^ ((A ^ B) if first else 0) ^ ((C ^ B) if second else 0)
        chosen = report['basis_degrees']
        rows, columns = [[0] * (pairs * pairs) for _ in chosen], [[0] * (pairs * pairs) for _ in chosen]
        for x, y, u, v in product(range(pairs), repeat=4):
            value = 0
            if x == y == u == v:
                value = H
            elif x == v and y == u:
                value = crossing[min(x, y), max(x, y)]
            elif x == u and y == v:
                value = K
            for i, degree in enumerate(chosen):
                if value >> degree & 1:
                    rows[i][x * pairs + y] ^= 1 << (v * pairs + u)
                    columns[i][x * pairs + u] ^= 1 << (v * pairs + y)
        ranks = list(map(list, zip(all_ranks(rows), all_ranks(columns))))
        digest = hashlib.sha256(json.dumps(ranks, separators=(',', ':')).encode()).hexdigest()
        assert digest == report['rank_table_sha256']
        records.append({'pairs': pairs, 'output_mask': mask, 'dense_entries': pairs ** 4,
                        'complete_binary_combinations': len(ranks), 'rank_table_sha256': digest})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = {'status': 'pass',
              'scope': 'All odd dimensions: polynomial block formula and all-field rank-charge theorem for the stated syndrome-preserving outputs; finite dense controls',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'literal_exterior_polynomial_checks': polynomial_checks(),
              'independent_dense_rank_checks': dense_checks(),
              'polynomial_encodings_p_stride4_q': {'H': H, 'K': K, 'A': A, 'B': B, 'C': C}}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps(report, sort_keys=True))


if __name__ == '__main__':
    main()
