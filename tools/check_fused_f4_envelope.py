#!/usr/bin/env python3
"""Independent dense F4 rank and literal coefficient controls for the envelope."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random
import subprocess
from tempfile import TemporaryDirectory

from check_fused_scalar_projection import DUAL, TARGET, dense_binary_rank, write_blocks
from f4_lattice_linear import INV, MUL, add, eliminate, encode, scale
from fused_f4_envelope import Envelope, normalize
from joint_energy_conjugations import Space


def dense_rank(matrix):
    # Ordinary lowest-column Gaussian elimination; no charge decomposition.
    matrix = [row[:] for row in matrix]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next((r for r in range(rank, len(matrix)) if matrix[r][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        factor = INV[matrix[rank][column]]
        matrix[rank] = [MUL[factor][entry] for entry in matrix[rank]]
        for r in range(rank+1, len(matrix)):
            factor = matrix[r][column]
            if factor:
                matrix[r] = [x ^ MUL[factor][y] for x, y in zip(matrix[r], matrix[rank])]
        rank += 1
    return rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--scanner', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    for values in product(range(4), repeat=4):
        matrix = [list(values[:2]), list(values[2:])]
        assert dense_rank(matrix) == len(eliminate([encode(row) for row in matrix]))
    envelope, rng = Envelope(), Random(0xF41601A195C2)
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root/'evidence/fused-f4-subspace-manifest.json').read_text())
    original = Space(manifest['original_binary_spanning_points'])
    declared = Space([manifest['offset']]+manifest['directions'])
    assert len(original.basis) == len(declared.basis) == 18
    assert len(manifest['directions']) == 17 and manifest['registered'] == 1 << 34
    for value in manifest['original_binary_spanning_points']:
        declared.binary_coordinates(value)
    for value in [manifest['offset']]+manifest['directions']:
        original.binary_coordinates(value)
    assert manifest['gray_xor_mask'] < 1 << 34
    assert (manifest['offset'] & DUAL).bit_count() % 2 == 1
    assert all((v & DUAL).bit_count() % 2 == 0 for v in manifest['directions'])
    with TemporaryDirectory(prefix='fused-f4-input-') as directory:
        rebuilt = Path(directory)/'blocks.bin'
        write_blocks(rebuilt, envelope.cuts, envelope.groups)
        assert rebuilt.read_bytes() == (root/'evidence/fused-f4-subspace-input.bin').read_bytes()
        assert sha256(rebuilt.read_bytes()).hexdigest() == manifest['block_input_sha256']
    points = [(TARGET, 0)]
    pivot = (DUAL & -DUAL).bit_length()-1
    for _ in range(64):
        low, high = rng.getrandbits(33), rng.getrandbits(33)
        low ^= (1 ^ ((low & DUAL).bit_count() % 2)) << pivot
        high ^= ((high & DUAL).bit_count() % 2) << pivot
        points.append((low, high))
    for _ in range(16):
        low = rng.getrandbits(33)
        low ^= (1 ^ ((low & DUAL).bit_count() % 2)) << pivot
        points.append((low, 0))
    records = []
    for point in points:
        ranks = envelope.ranks(point)
        for side, (low, high) in enumerate(envelope.matrices(point)):
            matrix = [[((low >> (81*r+c)) & 1) | (((high >> (81*r+c)) & 1) << 1)
                       for c in range(81)] for r in range(81)]
            assert dense_rank(matrix) == ranks[side]
            if not high:
                assert dense_binary_rank(low) == ranks[side]
        assert all(envelope.ranks(scale(point, scalar)) == ranks for scalar in (1, 2, 3))
        assert all(normalize(scale(point, scalar)) == normalize(point) for scalar in (1, 2, 3))
        records.append({'coordinates': point, 'ranks': ranks})
    projections = []
    for scalar in (1, 5, 9, 13):
        pieces = [envelope.projected_piece(i, j, scalar) for i, j in product(range(3), repeat=2)]
        projections.append({'scalar': scalar, 'pieces': pieces,
                            'complete_decomposition_ranks': envelope.decomposition(pieces)})
    native = []
    if args.binary:
        with TemporaryDirectory(prefix='fused-f4-controls-') as directory:
            directory = Path(directory)
            write_blocks(directory/'blocks.bin', envelope.cuts, envelope.groups)
            (directory/'points.txt').write_text(''.join(f'{a} {b}\n' for a, b in points))
            run = subprocess.run([str(args.binary.resolve()), str(directory/'blocks.bin'),
                                  str(directory/'points.txt'), '0', '1'],
                                 check=True, capture_output=True, text=True)
            native = [json.loads(line) for line in run.stdout.splitlines()]
            controls = [row for row in native if row['status'] == 'seed_control']
            assert len(controls) == len(records)
            assert all(tuple(a['coordinates']) == tuple(b['coordinates']) and
                       tuple(a['ranks']) == tuple(b['ranks']) for a, b in zip(controls, records))
    scan_controls = []
    if args.scanner:
        directions, mask = [816, 257557, 384425], 42
        assert all((DUAL & value).bit_count() % 2 == 0 for value in directions)
        with TemporaryDirectory(prefix='fused-f4-scan-controls-') as directory:
            directory = Path(directory)
            write_blocks(directory/'blocks.bin', envelope.cuts, envelope.groups)
            for offset, affine, threshold in ((106954755, 1, 6561), (106954755, 1, 501), (0, 0, 55)):
                text = f'{offset} 3\n'+'\n'.join(map(str, directions))+'\n'
                (directory/'subspace.txt').write_text(text)
                expected, zeros, skips = [], 0, 0
                for position in range(64):
                    code = (position ^ (position >> 1)) ^ mask
                    point = (offset, 0)
                    for j, value in enumerate(directions):
                        point = add(point, scale((value, 0), ((code >> j) & 1) | (((code >> (j+3)) & 1) << 1)))
                    if point == (0, 0):
                        zeros += 1
                    elif not affine and normalize(point) != point:
                        skips += 1
                    else:
                        ranks = envelope.ranks(point)
                        if ranks[0]*ranks[1] <= threshold:
                            expected.append({'position': position, 'coordinates': list(point), 'ranks': list(ranks)})
                run = subprocess.run([str(args.scanner.resolve()), str(directory/'blocks.bin'),
                                      str(directory/'subspace.txt'), str(directory/'retained.jsonl'),
                                      '0', '64', str(mask), str(threshold), str(affine)],
                                     check=True, capture_output=True, text=True)
                summary = json.loads(run.stdout)
                actual = [json.loads(line) for line in (directory/'retained.jsonl').read_text().splitlines()]
                assert actual == expected and summary['visited'] == 64
                assert summary['retained'] == len(expected) and summary['zeros'] == zeros
                assert summary['scalar_orbit_skips'] == skips
                assert summary['rejected']+len(expected)+zeros+skips == 64
                scan_controls.append({'affine_value': affine, 'threshold': threshold, 'summary': summary,
                                      'retained': actual})
    names = ['tools/check_fused_f4_envelope.py', 'tools/fused_f4_envelope.py',
             'tools/search_fused_f4_envelope.cpp', 'tools/check_fused_scalar_projection.py',
             'tools/scan_fused_f4_subspace.cpp',
             'tools/f4_lattice_linear.py', 'tools/joint_energy_conjugations.py']
    report = {'status': 'pass', 'all_two_by_two_f4_matrices': 256,
              'independent_full_81_by_81_ranks': 2*len(records),
              'exact_nonbinary_points': 64, 'records': records,
              'field_projection_decompositions': projections, 'native_controls': native,
              'complete_small_gray_scan_controls': scan_controls,
              'registered_subspace_dimension': 18, 'registered_vectors_per_affine_value': 1 << 34,
              'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest() for name in names},
              'scope': 'Representation, rank and literal decomposition controls only. '
              'No exhaustive F4 envelope exclusion or improved factoring exponent.'}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'all_two_by_two_f4_matrices', 'independent_full_81_by_81_ranks')}))


if __name__ == '__main__':
    main()
