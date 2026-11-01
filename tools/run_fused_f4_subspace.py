#!/usr/bin/env python3
"""Checkpoint the two registered finite F4 rank-product threshold scans."""
from __future__ import annotations

import argparse
from collections import Counter
import fcntl
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import time

from check_fused_f4_envelope import dense_rank
from check_fused_scalar_projection import DUAL
from fused_f4_envelope import Envelope, normalize


def write_json(path, data):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')
    temporary.replace(path)


def candidate(position, manifest, affine):
    code = (position ^ (position >> 1)) ^ manifest['gray_xor_mask']
    low, high = (manifest['offset'] if affine else 0), 0
    n = len(manifest['directions'])
    for j, direction in enumerate(manifest['directions']):
        if code >> j & 1:
            low ^= direction
        if code >> (j+n) & 1:
            high ^= direction
    return low, high


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--blocks', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--affine', type=int, choices=(0, 1), required=True)
    parser.add_argument('--chunk-size', type=int, default=1 << 25)
    parser.add_argument('--max-chunks', type=int)
    args = parser.parse_args()
    if args.chunk_size <= 0 or (args.max_chunks is not None and args.max_chunks <= 0):
        parser.error('Chunk budgets must be positive.')
    args.directory.mkdir(parents=True, exist_ok=True)
    lock = (args.directory/'run.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads(args.manifest.read_text())
    for name, expected in manifest['source_sha256'].items():
        assert sha256((root/name).read_bytes()).hexdigest() == expected, name
    assert sha256(args.binary.read_bytes()).hexdigest() == manifest['binary_sha256']
    assert sha256(args.blocks.read_bytes()).hexdigest() == manifest['block_input_sha256']
    assert manifest['registered'] == 1 << (2*len(manifest['directions']))
    threshold = 501 if args.affine else 55
    frozen = {name: (root/name).read_text() for name in manifest['source_sha256']}
    registration = {'manifest': manifest, 'affine_value': args.affine,
                    'threshold': threshold, 'chunk_size': args.chunk_size,
                    'as_run_sources': frozen}
    if (args.directory/'registration.json').exists():
        assert json.loads((args.directory/'registration.json').read_text()) == registration
    else:
        write_json(args.directory/'registration.json', registration)
    declaration = f"{manifest['offset'] if args.affine else 0} {len(manifest['directions'])}\n"
    declaration += '\n'.join(map(str, manifest['directions']))+'\n'
    (args.directory/'subspace.txt').write_text(declaration)
    envelope, ranks_seen = Envelope(), {}
    counts = Counter(visited=0, retained=0, rejected=0, zeros=0, scalar_orbit_skips=0)
    histogram, position, chunks, elapsed = Counter(), 0, [], 0.

    def consume(summary, retained_path):
        nonlocal position, elapsed
        assert summary['status'] == 'complete' and summary['start'] == position
        assert summary['stop'] > position and summary['visited'] == summary['stop']-position
        assert summary['visited'] == sum(summary[k] for k in
                ('retained', 'rejected', 'zeros', 'scalar_orbit_skips'))
        assert sha256(retained_path.read_bytes()).hexdigest() == summary['retained_sha256']
        records = [json.loads(line) for line in retained_path.read_text().splitlines()]
        assert len(records) == summary['retained']
        previous = position-1
        for record in records:
            p = record['position']
            assert previous < p < summary['stop'] and p >= position
            previous = p
            point, ranks = tuple(record['coordinates']), tuple(record['ranks'])
            assert point == candidate(p, manifest, args.affine)
            assert ((point[0] & DUAL).bit_count() % 2) == args.affine
            assert (point[1] & DUAL).bit_count() % 2 == 0 and point != (0, 0)
            assert args.affine or normalize(point) == point
            assert len(ranks) == 2 and all(type(r) is int and 1 <= r <= 81 for r in ranks)
            assert ranks[0]*ranks[1] <= threshold
            histogram[ranks] += 1
            if ranks not in ranks_seen:
                independent = []
                for low, high in envelope.matrices(point):
                    matrix = [[((low >> (81*r+c)) & 1) | (((high >> (81*r+c)) & 1) << 1)
                               for c in range(81)] for r in range(81)]
                    independent.append(dense_rank(matrix))
                assert tuple(independent) == ranks
                ranks_seen[ranks] = record
        for key in counts:
            counts[key] += summary[key]
        position, elapsed = summary['stop'], elapsed+summary['seconds']
        chunks.append(summary)

    def save():
        complete = position == manifest['registered']
        if complete:
            if args.affine:
                assert counts['zeros'] == counts['scalar_orbit_skips'] == 0
            else:
                assert counts['zeros'] == 1
                assert counts['scalar_orbit_skips'] == 2*(manifest['registered']-1)//3
        report = {'status': 'complete' if complete else 'incomplete', 'registered': manifest['registered'],
                  'affine_value': args.affine, 'threshold': threshold, **counts,
                  'native_seconds': elapsed, 'chunks': len(chunks),
                  'histogram': [{'ranks': key, 'count': value} for key, value in sorted(histogram.items())],
                  'independent_full_rank_class_representatives': list(ranks_seen.values()),
                  'scope': 'Coverage of the declared proper F4 subspace only. Rejections use monotone '
                  'partial rank products. Scalar-orbit skips in affine zero preserve both ranks. '
                  'No full-envelope, all-field or improved factoring claim.'}
        write_json(args.directory/'summary.json', report)
        return report

    for path in sorted(args.directory.glob('chunk-*.json')):
        consume(json.loads(path.read_text()), path.with_suffix('.retained.jsonl'))
    save()
    fresh = 0
    while position < manifest['registered'] and (args.max_chunks is None or fresh < args.max_chunks):
        stop = min(position+args.chunk_size, manifest['registered'])
        prefix = args.directory/f'chunk-{position:012d}'
        retained = prefix.with_suffix('.retained.jsonl')
        temporary = retained.with_suffix('.tmp')
        run = subprocess.run([str(args.binary.resolve()), str(args.blocks.resolve()),
                              str((args.directory/'subspace.txt').resolve()), str(temporary.resolve()),
                              str(position), str(stop), str(manifest['gray_xor_mask']),
                              str(threshold), str(args.affine)], check=True, capture_output=True, text=True)
        summary = json.loads(run.stdout)
        summary['retained_sha256'] = sha256(temporary.read_bytes()).hexdigest()
        consume(summary, temporary)
        temporary.replace(retained)
        write_json(prefix.with_suffix('.json'), summary)
        report = save()
        print(json.dumps({key: report[key] for key in
                          ('status', 'affine_value', 'visited', 'registered', 'retained', 'chunks', 'native_seconds')}), flush=True)
        fresh += 1


if __name__ == '__main__':
    main()
