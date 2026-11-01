#!/usr/bin/env python3
"""Check the completed registered affine-zero F4 threshold census ledger.

Reconstructs every chunk's zero/normalization counts independently by
linear algebra. Does not repeat the billions of native rank tests.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import lzma
from pathlib import Path
from random import Random
from tempfile import TemporaryDirectory

from check_fused_scalar_projection import DUAL, coefficient_space, partitions, write_blocks
from fused_f4_coverage import affine_zero_counts
from joint_energy_conjugations import Space


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = Random(0xF4C0A195C2)
    for n in range(1, 5):
        for _ in range(16):
            directions = [rng.getrandbits(9) for _ in range(n)]
            mask, free = rng.getrandbits(2*n), rng.randrange(2*n+1)
            length = 1 << free
            start = rng.randrange(1 << (2*n-free))*length
            actual = {'visited': length, 'normalized_threshold_tests': 0, 'zeros': 0, 'scalar_orbit_skips': 0}
            for i in range(start, start+length):
                code, low, high = (i ^ (i >> 1)) ^ mask, 0, 0
                for j, direction in enumerate(directions):
                    if code >> j & 1:
                        low ^= direction
                    if code >> (j+n) & 1:
                        high ^= direction
                if not low | high:
                    actual['zeros'] += 1
                else:
                    pivot = (low | high).bit_length()-1
                    if ((low >> pivot) & 1) | (((high >> pivot) & 1) << 1) == 1:
                        actual['normalized_threshold_tests'] += 1
                    else:
                        actual['scalar_orbit_skips'] += 1
            assert actual == affine_zero_counts(start, start+length, directions, mask)
    root = Path(__file__).resolve().parents[1]
    archive = root/'evidence/fused-f4-zero-census.json.xz'
    packet = json.loads(lzma.decompress(archive.read_bytes()))
    registration, summary = packet['registration'], packet['summary']
    manifest = registration['manifest']
    total, directions, mask = manifest['registered'], manifest['directions'], manifest['gray_xor_mask']
    assert total == 1 << 34 and len(directions) == 17
    assert len(Space(directions).basis) == 17
    assert len(Space([manifest['offset']]+directions).basis) == 18
    assert (manifest['offset'] & DUAL).bit_count() % 2 == 1
    assert all((v & DUAL).bit_count() % 2 == 0 for v in directions)
    original = Space(manifest['original_binary_spanning_points'])
    declared = Space([manifest['offset']]+directions)
    assert len(original.basis) == 18
    for v in manifest['original_binary_spanning_points']:
        declared.binary_coordinates(v)
    for v in [manifest['offset']]+directions:
        original.binary_coordinates(v)
    assert registration['affine_value'] == summary['affine_value'] == 0
    assert registration['threshold'] == summary['threshold'] == 55
    for name, expected in manifest['source_sha256'].items():
        assert sha256(registration['as_run_sources'][name].encode()).hexdigest() == expected
    cuts, _ = coefficient_space(ordered=True)
    with TemporaryDirectory(prefix='fused-f4-zero-input-') as temporary:
        path = Path(temporary)/'blocks.bin'
        write_blocks(path, cuts, partitions())
        assert sha256(path.read_bytes()).hexdigest() == manifest['block_input_sha256']
        assert path.read_bytes() == (root/'evidence/fused-f4-subspace-input.bin').read_bytes()
    position, counts, chunk_controls = 0, dict(visited=0, normalized_threshold_tests=0, zeros=0, scalar_orbit_skips=0), []
    for record in packet['chunks']:
        data, retained = record['summary'], record['retained_text']
        encoded = (json.dumps(data, indent=2, sort_keys=True)+'\n').encode()
        assert sha256(encoded).hexdigest() == record['summary_sha256']
        assert sha256(retained.encode()).hexdigest() == data['retained_sha256']
        assert retained == '' and data['retained'] == 0
        assert data['status'] == 'complete' and data['start'] == position
        independent = affine_zero_counts(position, data['stop'], directions, mask)
        assert data['visited'] == independent['visited']
        assert data['zeros'] == independent['zeros']
        assert data['scalar_orbit_skips'] == independent['scalar_orbit_skips']
        assert data['rejected'] == independent['normalized_threshold_tests']
        for key in counts:
            counts[key] += independent[key]
        chunk_controls.append({'start': position, 'stop': data['stop'], **independent})
        position = data['stop']
    assert position == total and len(packet['chunks']) == summary['chunks'] == 512
    assert summary['status'] == 'complete' and summary['visited'] == summary['registered'] == total
    assert counts['normalized_threshold_tests'] == summary['rejected'] == (total-1)//3
    assert counts['zeros'] == summary['zeros'] == 1
    assert counts['scalar_orbit_skips'] == summary['scalar_orbit_skips'] == 2*(total-1)//3
    assert summary['retained'] == 0 and not summary['histogram']
    report = {'status': 'pass', 'covered_coordinate_positions': total,
              'exact_native_threshold_tests': counts['normalized_threshold_tests'],
              'independent_chunk_normalization_counts': chunk_controls,
              'scalar_orbit_skips': counts['scalar_orbit_skips'], 'zero_vectors': 1,
              'nonzero_rank_product_floor_in_declared_affine_zero_space': 56,
              'small_independent_bruteforce_count_controls': 64,
              'archive_sha256': sha256(archive.read_bytes()).hexdigest(),
              'scope': 'Checks the entire registered ledger and each chunk\'s normalization/zero '
              'counts by independent binary linear algebra. Native exact threshold rejections '
              'are not all independently recomputed. The result concerns the specified 17-dimensional '
              'F4 kernel in W18, not the full 33-dimensional envelope or other fields; no score change.',
              'source_sha256': {n: sha256((root/n).read_bytes()).hexdigest() for n in
                                ['tools/check_fused_f4_zero_census.py', 'tools/fused_f4_coverage.py',
                                 'tools/check_fused_scalar_projection.py', 'tools/joint_energy_conjugations.py']}}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'covered_coordinate_positions', 'exact_native_threshold_tests',
                       'nonzero_rank_product_floor_in_declared_affine_zero_space')}))


if __name__ == '__main__':
    main()
