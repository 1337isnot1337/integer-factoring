#!/usr/bin/env python3
"""Exact finite raw-coefficient ranks on non-affine quadratic Boolean codes."""

import argparse
from functools import lru_cache
import json
from math import log2
from pathlib import Path
import time
from quadratic_crossing import (
    binary_rank, binomial_profile_from_histogram, expanded_profile,
    label_polynomial, truth_word,
)


@lru_cache(maxsize=None)
def profile_for_histogram(histogram):
    return expanded_profile(binomial_profile_from_histogram(histogram))


def matrices(variables, width):
    dimension = 1 << width
    length = 1 << variables
    all_bits = (1 << length) - 1
    words = [truth_word(variables, label_polynomial(variables, width, i))
             for i in range(dimension)]
    pair_masks = []
    for a in words:
        row = []
        for b in words:
            row.append((all_bits ^ (a | b), (~a) & b, a & (~b), a & b))
        pair_masks.append(row)
    rows, columns = {}, {}
    start = time.monotonic()
    for x in range(dimension):
        for y in range(dimension):
            first = pair_masks[x][y]
            row_id = x * dimension + y
            for u in range(dimension):
                column_row = x * dimension + u
                for v in range(dimension):
                    second = pair_masks[u][v]
                    histogram = tuple((a & b).bit_count() for a in first for b in second)
                    polynomial = profile_for_histogram(histogram)
                    assert polynomial.bit_count() % 2 == int(x == v and y == u)
                    row_bit = 1 << (u * dimension + v)
                    column_bit = 1 << (y * dimension + v)
                    while polynomial:
                        bit = polynomial & -polynomial
                        degree = bit.bit_length() - 1
                        rows.setdefault(degree, [0] * dimension**2)[row_id] |= row_bit
                        columns.setdefault(degree, [0] * dimension**2)[column_row] |= column_bit
                        polynomial ^= bit
        if x % max(1, dimension // 8) == 0:
            print(json.dumps({'x': x, 'dimension': dimension,
                              'seconds': round(time.monotonic() - start, 3),
                              'histograms': profile_for_histogram.cache_info().currsize}), flush=True)
    return rows, columns


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variables', type=int, default=3)
    parser.add_argument('--width', type=int, default=5)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    rows, columns = matrices(args.variables, args.width)
    profiles = []
    for degree in sorted(rows):
        r, c = binary_rank(rows[degree]), binary_rank(columns[degree])
        profiles.append({'degree': degree, 'row_rank': r, 'column_rank': c,
                         'rank_product': r * c})
        print(json.dumps(profiles[-1]), flush=True)
    result = {'variables': args.variables, 'binary_width': args.width,
              'labels': 1 << args.width, 'evaluation_positions': 1 << args.variables,
              'raw_coefficient_profiles': profiles,
              'maximum_rank_product': max(p['rank_product'] for p in profiles),
              'finite_maximum_rho': max(log2(p['rank_product']) for p in profiles) / (4 * args.width),
              'seconds': time.monotonic() - start,
              'scope': 'Finite literal coefficient ranks; no asymptotic or optimized-span claim'}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
