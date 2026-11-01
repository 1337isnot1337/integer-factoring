#!/usr/bin/env python3
"""Exact quadratic-over-F4 place ranks for a frozen affine-path sample.

One representative of each F16/F4 conjugate pair is evaluated.  Field
restriction doubles every matrix rank; this is independently checked
against direct F16 elimination.  Ranks alone are not pressure bounds.
"""
from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import time

from check_fused_multipoint_profiles import field_rank
from f16_identity_operator_profiles import MUL
from fused_six_vertex_crossings import power
from physical_split_crossings import F4_PRODUCT
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import scale
from screen_affine_finite_places import bitplane_rank
from search_affine_coupled_reflection_frames import affine_seed, catalogue


@lru_cache(maxsize=None)
def evaluate(value, point):
    result = 0
    for component, multiplier in zip(value, (1, 2)):
        while component:
            bit = component & -component
            result ^= MUL[multiplier][power(point, bit.bit_length() - 1)]
            component ^= bit
    return result


def restricted_rows(values, cut, point):
    result = [(0, 0)] * 32
    for index, value in enumerate(values):
        x, y, u, v = ((index // 4 ** leg) % 4 for leg in range(4))
        row, column = (4 * x + y, 4 * u + v) if cut == 0 else (4 * x + u, 4 * y + v)
        coefficient = evaluate(tuple(value), point)
        a, b = coefficient & 3, coefficient >> 2
        block = ((a, F4_PRODUCT[2][b]), (b, a ^ b))
        for i in (0, 1):
            low, high = result[2 * row + i]
            for j in (0, 1):
                low ^= (block[i][j] & 1) << (2 * column + j)
                high ^= (block[i][j] >> 1) << (2 * column + j)
            result[2 * row + i] = low, high
    return result


def direct_rank(values, cut, point):
    matrix = [[0] * 16 for _ in range(16)]
    for index, value in enumerate(values):
        x, y, u, v = ((index // 4 ** leg) % 4 for leg in range(4))
        row, column = (4 * x + y, 4 * u + v) if cut == 0 else (4 * x + u, 4 * y + v)
        matrix[row][column] = evaluate(tuple(value), point)
    return field_rank(matrix)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('snapshot', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.snapshot.read_bytes()
    snapshot = json.loads(raw)
    directory = Path(__file__).parent
    for name, expected in snapshot['source_sha256'].items():
        if sha256((directory / name).read_bytes()).hexdigest() != expected:
            raise ArithmeticError('A frozen proposal source changed: ' + name)
    base = tensor(4)
    frames = catalogue(base)
    assert [{k: v for k, v in frame.items() if k != 'delta'} for frame in frames] == snapshot['frame_catalogue']
    orbits = [sorted((point, power(point, 4))) for point in range(4, 16) if point < power(point, 4)]
    points = [orbit[0] for orbit in orbits]
    caches = {(point, cut): (restricted_rows(base['values'], cut, point),
                             [restricted_rows(frame['delta'], cut, point) for frame in frames])
              for point in points for cut in (0, 1)}
    records, histogram, controls = [], Counter(), 0
    started = time.monotonic()
    for index, row in enumerate(snapshot['screen_records']):
        ranks = []
        for point in points:
            pair = []
            for cut in (0, 1):
                original, deltas = caches[point, cut]
                matrix = list(original)
                for frame, coefficient in row['word']:
                    matrix = [(a ^ c, b ^ d) for (a, b), (c, d) in
                              zip(matrix, [scale(value, coefficient) for value in deltas[frame]])]
                result = bitplane_rank(matrix)
                assert result % 2 == 0
                pair.append(result // 2)
            ranks.append(pair)
        signature = tuple(value for pair in ranks for value in pair)
        histogram[signature] += 1
        if min(signature) < 16:
            records.append({'snapshot_index': index, 'word': row['word'], 'seed_sha256': row['seed_sha256'],
                            'natural_rho': row['natural_rho'], 'ranks': ranks})
        if index < 8:
            seed = affine_seed(base, frames, row['word'])
            for point, pair in zip(points, ranks):
                for cut in (0, 1):
                    assert direct_rank(seed['values'], cut, point) == pair[cut]
                    controls += 1
    records.sort(key=lambda row: (min(v for pair in row['ranks'] for v in pair),
                                  sum(v for pair in row['ranks'] for v in pair)))
    names = (Path(__file__).name, 'screen_affine_finite_places.py', 'search_affine_coupled_reflection_frames.py',
             'reflection_laurent_physical_frames.py', 'check_fused_multipoint_profiles.py',
             'f16_identity_operator_profiles.py', 'fused_six_vertex_crossings.py')
    result = {'scope': 'Exact finite ranks at all quadratic closed places represented in F16; no complete-profile or pressure claim.',
              'snapshot_sha256': sha256(raw).hexdigest(), 'screened': len(snapshot['screen_records']),
              'conjugate_pairs': orbits, 'representative_points': points,
              'independent_F16_rank_controls': controls, 'singular_candidates': len(records),
              'rank_histogram': [{'signature': list(key), 'count': value} for key, value in sorted(histogram.items())],
              'retained': records, 'elapsed_seconds': time.monotonic() - started,
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names}}
    args.output.write_text(json.dumps(result, separators=(',', ':')) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('retained', 'rank_histogram', 'source_sha256')}))


if __name__ == '__main__':
    main()
