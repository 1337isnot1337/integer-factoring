#!/usr/bin/env python3
"""Stop exact paired-frame closure at its proved 5175-dimensional ceiling."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
from hashlib import sha256
import json
from pathlib import Path
import time

from close_independent_paired_affine_span import operators, image, insert
from dense_coupled_reflection_frames import save
from reflection_dual_rail_gaussian import tensor
from solve_companion_affine_degree import coefficient, full_vector, unpack
from check_degree16_reflection_crossing import paired_traces


def canonical_basis():
    traceless = []
    for i in range(4):
        for j in range(4):
            if i != j:
                traceless.append([(i, j)])
    for i in range(1, 4):
        traceless.append([(0, 0), (i, i)])
    rows = []
    for first in traceless:
        for second in traceless:
            entries = [x + 4 * y + 16 * u + 64 * v
                       for x, v in first for y, u in second]
            for power in range(23):
                low = sum((1 << (power + 25 * index)) | (1 << (power + 2 + 25 * index))
                          for index in entries)
                rows.append((low, 0))
    assert len(rows) == 5175
    assert len({a.bit_length() - 1 for a, b in rows}) == 5175
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('basis', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--basis-output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    started_utc = datetime.now(timezone.utc).isoformat()
    raw = args.basis.read_bytes()
    original = json.loads(gzip.decompress(raw))
    for name, expected in original['source_sha256'].items():
        assert sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() == expected
    rows = {row['pivot']: tuple(int(x, 16) for x in row['bitplanes']) for row in original['basis']}
    assert len(rows) == original['dimension'] == len(original['basis'])
    for pivot, vector in rows.items():
        assert 0 <= pivot < 6400
        assert (vector[0] | vector[1]).bit_length() - 1 == pivot
        assert coefficient(vector, pivot) == 1
    # Highest polynomial-coordinate pivots expose the missing directions
    # quickly; the exact ambient ceiling is independent of this order.
    queue = [full_vector(tensor(4)['values'])] + [
        vector for _, vector in sorted(rows.items(), reverse=True)
    ]
    cursor, images = 0, 0
    ops = operators()
    while len(rows) < 5175:
        assert cursor < len(queue), 'The exact closure did not reach the ambient dimension.'
        for operator in ops:
            result = insert(rows, image(queue[cursor], operator))
            images += 1
            if result is not None:
                queue.append(result[1])
                assert len(rows) <= 5175
        cursor += 1
    # Every independent stored direction lies in the stated ambient space.
    odd_mask = sum(1 << i for i in range(1, 25, 2))
    for vector in rows.values():
        values = unpack(vector)
        assert all(not (a.bit_count() % 2 or b.bit_count() % 2) for a, b in values)
        assert all(not ((a & odd_mask).bit_count() % 2 or (b & odd_mask).bit_count() % 2)
                   for a, b in values)
        assert all(not any(a or b for a, b in matrix) for matrix in paired_traces(values))
    canonical = canonical_basis()
    # Exact membership of the sparse canonical basis in the constructed
    # row span is implied by equal dimensions and the ambient checks above.
    basis = {
        'scope': 'canonical complete degree24 paired-trace affine space, proved attained by physical-frame closure',
        'dimension': 5175, 'ambient_field_coordinates': 6400,
        'packing': 'field coefficient at entry*25+degree in two F4 bitplanes',
        'basis': [{'pivot': a.bit_length() - 1, 'bitplanes': [hex(a), '0x0']} for a, b in canonical],
        'input_basis_sha256': sha256(raw).hexdigest(),
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, 'close_independent_paired_affine_span.py',
                                       'solve_companion_affine_degree.py',
                                       'check_degree16_reflection_crossing.py',
                                       'reflection_dual_rail_gaussian.py',
                                       'dense_coupled_reflection_frames.py')},
    }
    encoded = json.dumps(basis, separators=(',', ':')).encode()
    compressed = gzip.compress(encoded, mtime=0)
    args.basis_output.parent.mkdir(parents=True, exist_ok=True)
    args.basis_output.write_bytes(compressed)
    certificate = {
        'scope': 'exact full paired-frame span saturation; no pressure certificate',
        'started_utc': started_utc, 'finished_utc': datetime.now(timezone.utc).isoformat(),
        'input_basis_sha256': sha256(raw).hexdigest(), 'initial_dimension': original['dimension'],
        'exact_constructed_dimension': len(rows), 'ambient_dimension': 225 * 23,
        'queue_vectors_processed': cursor, 'operator_images_processed': images,
        'zero_endpoint_and_derivative_entries_checked': len(rows) * 256,
        'zero_paired_trace_entries_checked': len(rows) * 32,
        'canonical_basis_dimension': len(canonical),
        'canonical_basis_sha256': sha256(compressed).hexdigest(),
        'canonical_basis_bytes': len(compressed), 'elapsed_seconds': time.monotonic() - started,
        'source_sha256': basis['source_sha256'],
    }
    save(args.output, certificate)
    print(json.dumps({key: value for key, value in certificate.items() if key != 'source_sha256'}))


if __name__ == '__main__':
    main()
