#!/usr/bin/env python3
"""Replay the complete declared F4 line and affine-coset control packets."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
import lzma
from pathlib import Path

from check_fused_f4_envelope import dense_rank
from check_fused_scalar_projection import TARGET
from f4_lattice_linear import add, scale
from fused_f4_envelope import Envelope, normalize
from joint_energy_conjugations import Space


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    archive = root/'evidence/fused-f4-finite-controls.json.xz'
    packet = json.loads(lzma.decompress(archive.read_bytes()))
    seeds = [0, TARGET]+[r['coordinates'] for r in packet['predecessor_data']['complement_pairs']]
    expected = sorted({normalize((a, b)) for a, b in product(seeds, repeat=2)}-{(0, 0)})
    lines = packet['lines']
    assert lines['manifest']['registered_count'] == lines['results']['visited'] == len(expected) == 484
    assert [tuple(r['coordinates']) for r in lines['results']['records']] == expected
    envelope, line_ranks = Envelope(), 0
    for record in lines['results']['records']:
        point = tuple(record['coordinates'])
        ranks = envelope.ranks(point)
        assert ranks == tuple(record['ranks'])
        line_ranks += 2
        for orbit in record['two_piece_orbit']:
            complement = add((TARGET, 0), scale(point, orbit['scalar']))
            assert envelope.ranks(complement) == tuple(orbit['complement_ranks'])
            line_ranks += 2
    cosets = packet['cosets']
    manifest, result = cosets['manifest'], cosets['results']
    directions = manifest['binary_direction_basis']
    assert len(Space(directions).basis) == len(directions) == 6
    space = Space(directions)
    def canonical(value):
        for p in sorted(space.pivots, reverse=True):
            if value >> p & 1:
                value ^= space.pivots[p][0]
        return value
    representatives = sorted({canonical(value) for value in seeds[2:]})
    assert representatives == manifest['affine_coset_representatives'] and len(representatives) == 15
    expected_points = []
    for representative in representatives:
        for coefficients in product(range(4), repeat=6):
            point = (representative, 0)
            for value, coefficient in zip(directions, coefficients):
                point = add(point, scale((value, 0), coefficient))
            expected_points.append(point)
    assert len(expected_points) == len(set(expected_points)) == 61440
    text = ''.join(f'{a} {b}\n' for a, b in expected_points)
    assert sha256(text.encode()).hexdigest() == manifest['candidate_text_sha256']
    records = [json.loads(line) for line in cosets['native_rows'].splitlines()]
    records = [r for r in records if r['status'] == 'seed_control']
    assert len(records) == result['visited'] == manifest['registered'] == 61440
    assert [tuple(r['coordinates']) for r in records] == expected_points
    histogram = Counter(tuple(r['ranks']) for r in records)
    assert histogram == Counter({tuple(r['ranks']): r['count'] for r in result['histogram']})
    assert len(histogram) == 696 and min(a*b for a, b in histogram) == result['minimum_product'] == 162
    independent = result['independently_replayed_rank_class_representatives']
    assert {tuple(r['ranks']) for r in independent} == set(histogram)
    record_lookup = {tuple(r['coordinates']): tuple(r['ranks']) for r in records}
    for record in independent:
        point, ranks = tuple(record['coordinates']), tuple(record['ranks'])
        assert record_lookup[point] == ranks
        for side, (low, high) in enumerate(envelope.matrices(point)):
            matrix = [[((low >> (81*r+c)) & 1) | (((high >> (81*r+c)) & 1) << 1)
                       for c in range(81)] for r in range(81)]
            assert dense_rank(matrix) == ranks[side]
    report = {'status': 'pass', 'projective_line_points': 484,
              'complete_affine_coset_points': 61440, 'coset_rank_classes': 696,
              'independent_full_matrix_ranks': 1392, 'line_and_complement_block_ranks': line_ranks,
              'archive_sha256': sha256(archive.read_bytes()).hexdigest(),
              'scope': 'Checks both finite coverage ledgers, all line/complement ranks and every '
              'coset rank-class representative with independent full dense elimination. Does not '
              'repeat every native coset rank or certify the full F4 envelope. Score remains 1/2.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
