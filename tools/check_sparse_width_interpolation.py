#!/usr/bin/env python3
"""Finite tensor/rank checks for the sparse-width interpolation proof."""

import argparse
from bisect import bisect_right
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_modular_recombination import crossing, permute, ranks
from hamming_pfaffian import matrices
from check_partition_channels import rank_mod

ROOT = Path(__file__).resolve().parents[1]


def blocks(length, available):
    result = []
    while length:
        width = available[bisect_right(available, length)-1]
        result.append(width)
        length -= width
    return result


def pieces(width, case):
    source = matrices(width)
    d = 1 << width
    result = []
    for piece in case['witness']:
        value = [0]*(d*d)
        for term in piece['terms']:
            image = permute(source[term['shift']], d, term['permutation'])
            value = [a ^ b for a, b in zip(value, image)]
        assert ranks(value, d) == (piece['row_rank'], piece['column_rank'])
        result.append(value)
    return result


def kronecker(left, right):
    rows = []
    for a in left:
        for b in right:
            rows.append([x*y for x in a for y in b])
    return rows


def dense(rows, size):
    return [[(row >> j) & 1 for j in range(size)] for row in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    geometric = [1 << j for j in range(20)]
    polynomial = [1, 2, 4, 16, 256, 65536, 4294967296]
    schedule_cases = 0
    for length in list(range(1, 4097))+[65535, 65536, 65537, 1000000]:
        for available, power in [(geometric, 1), (polynomial, 2)]:
            chosen = blocks(length, available)
            assert sum(chosen) == length and all(b in available for b in chosen)
            if power == 1:
                assert len(chosen) <= length.bit_length()
            else:
                # Here K=2,s=2, so (4) implies k <= 2sqrt(2L).
                assert len(chosen)**2 <= 8*length
            schedule_cases += 1
    evidence = json.loads((ROOT/'evidence/hamming-pfaffian.json').read_text())
    cases = {case['width']: case for case in evidence['full_span_cases']}
    left, right = pieces(1, cases[1]), pieces(2, cases[2])
    total = [0]*64
    tested = 0
    for a, b in product(left, right):
        tensor = [0]*64
        for x, y, u, v in product(range(8), repeat=4):
            entry = ((a[(x//4)*2+y//4] >> ((u//4)*2+v//4)) & 1) & (
                (b[(x%4)*4+y%4] >> ((u%4)*4+v%4)) & 1)
            tensor[x*8+y] |= entry << (u*8+v)
        ra, ca = ranks(a, 2)
        rb, cb = ranks(b, 4)
        assert ranks(tensor, 8) == (ra*rb, ca*cb)
        for permutation, expected in [((0, 1, 2, 3), ra*rb), ((0, 2, 1, 3), ca*cb)]:
            da = dense(permute(a, 2, permutation), 4)
            db = dense(permute(b, 4, permutation), 16)
            assert rank_mod(kronecker(da, db), 2) == expected
        total = [a ^ b for a, b in zip(total, tensor)]
        tested += 1
    assert tested == 18 and total == crossing(8)
    result = {'status': 'pass; sparse-width theorem is conditional on source rank bounds',
              'greedy_schedule_cases': schedule_cases, 'product_pieces_checked': tested,
              'crossing_entries_checked': 4096,
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
