#!/usr/bin/env python3
"""Exact F4 rational-place rank screen of a frozen affine-path snapshot.

Finite ranks are only a proposal filter. They are not complete local Smith
profiles or a pressure certificate. Every retained word can be reconstructed
from the frozen exact frame catalogue before a full-cost follow-up.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import time

from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import scale
from search_affine_coupled_reflection_frames import affine_seed, catalogue
from square_zero_coupled_reflection_frames import rank


def evaluate(value, point):
    power, result = 1, 0
    for exponent in range(max(value[0].bit_length(), value[1].bit_length())):
        coefficient = ((value[0] >> exponent) & 1) + 2 * ((value[1] >> exponent) & 1)
        result ^= F4_PRODUCT[coefficient][power]
        power = F4_PRODUCT[power][point]
    return result


def rows(values, cut, point):
    output = [(0, 0)] * 16
    for index, polynomial in enumerate(values):
        x, y, u, v = ((index // 4 ** leg) % 4 for leg in range(4))
        row, column = (4 * x + y, 4 * u + v) if cut == 0 else (4 * x + u, 4 * y + v)
        value = evaluate(polynomial, point)
        a, b = output[row]
        output[row] = a ^ ((value & 1) << column), b ^ ((value >> 1) << column)
    return output


def bitplane_rank(matrix):
    pivots = {}
    for a, b in matrix:
        while a | b:
            bit = ((a | b) & -(a | b)).bit_length() - 1
            coefficient = ((a >> bit) & 1) + 2 * ((b >> bit) & 1)
            if bit in pivots:
                first, second = scale(pivots[bit], coefficient)
                a, b = a ^ first, b ^ second
            else:
                pivots[bit] = scale((a, b), (0, 1, 3, 2)[coefficient])
                break
    return len(pivots)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('snapshot', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--points', type=int, nargs='+', default=[2, 3])
    args = parser.parse_args()
    if any(point not in range(4) or point == 1 for point in args.points):
        parser.error('Use F4 places other than the crossing endpoint one.')
    raw = args.snapshot.read_bytes()
    snapshot = json.loads(raw)
    base = tensor(4)
    frames = catalogue(base)
    metadata = [{key: value for key, value in frame.items() if key != 'delta'} for frame in frames]
    assert metadata == snapshot['frame_catalogue']
    directory = Path(__file__).parent
    for name, expected in snapshot['source_sha256'].items():
        assert sha256((directory / name).read_bytes()).hexdigest() == expected
    caches = {(point, cut): (rows(base['values'], cut, point),
                             [rows(frame['delta'], cut, point) for frame in frames])
              for point in args.points for cut in (0, 1)}
    retained, histograms, independent_controls = [], Counter(), 0
    started = time.monotonic()
    for index, record in enumerate(snapshot['screen_records']):
        ranks = []
        for point in args.points:
            pair = []
            for cut in (0, 1):
                base_rows, delta_rows = caches[point, cut]
                matrix = list(base_rows)
                for frame, coefficient in record['word']:
                    matrix = [(a ^ c, b ^ d) for (a, b), (c, d) in
                              zip(matrix, [scale(row, coefficient) for row in delta_rows[frame]])]
                result = bitplane_rank(matrix)
                if index < 12:
                    dense = [[((a >> j) & 1) + 2 * ((b >> j) & 1) for j in range(16)] for a, b in matrix]
                    assert rank(dense) == result
                    independent_controls += 1
                pair.append(result)
            ranks.append(pair)
        signature = tuple(value for pair in ranks for value in pair)
        histograms[signature] += 1
        if min(signature) < 16:
            retained.append({'snapshot_index': index, 'word': record['word'],
                             'seed_sha256': record['seed_sha256'], 'natural_rho': record['natural_rho'],
                             'ranks': ranks})
    retained.sort(key=lambda row: (min(value for pair in row['ranks'] for value in pair),
                                   sum(value for pair in row['ranks'] for value in pair)))
    reconstruction_entries = 0
    for record in retained[:24]:
        seed = affine_seed(base, frames, record['word'])
        assert sha256(json.dumps(seed['values'], separators=(',', ':')).encode()).hexdigest() == record['seed_sha256']
        for point, observed in zip(args.points, record['ranks']):
            for cut in (0, 1):
                assert bitplane_rank(rows(seed['values'], cut, point)) == observed[cut]
                reconstruction_entries += 16 ** 2
    names = (Path(__file__).name, 'search_affine_coupled_reflection_frames.py',
             'reflection_laurent_physical_frames.py', 'reflection_dual_rail_gaussian.py',
             'physical_split_crossings.py', 'square_zero_coupled_reflection_frames.py')
    output = {'scope': 'Exact finite-place ranks for the frozen stated sample; no local-profile, optimality, or pressure conclusion.',
              'snapshot_sha256': sha256(raw).hexdigest(), 'snapshot_bytes': len(raw),
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
              'points': args.points, 'screened_seeds': len(snapshot['screen_records']),
              'singular_seeds': len(retained), 'independent_rank_controls': independent_controls,
              'full_coefficient_reconstruction_entries': reconstruction_entries,
              'rank_signature_fields': [[point, cut] for point in args.points for cut in (0, 1)],
              'histogram': [{'ranks': list(key), 'count': count} for key, count in sorted(histograms.items())],
              'retained': retained, 'elapsed_seconds': time.monotonic() - started}
    args.output.write_text(json.dumps(output, separators=(',', ':')) + '\n')
    print(json.dumps({'screened': output['screened_seeds'], 'singular': len(retained),
                      'best_ranks': retained[0]['ranks'] if retained else None,
                      'seconds': output['elapsed_seconds']}))


if __name__ == '__main__':
    main()
