#!/usr/bin/env python3
"""Rebuild or inspect the exact seven- and eight-letter fused endpoint seeds."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from fused_six_vertex_crossings import construct
from f16_identity_operator_profiles import evaluate_one


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--rebuild', action='store_true',
                        help='Repeat the full-qubit invariant-subspace construction, about a minute at eight letters.')
    args = parser.parse_args()
    directory = Path(__file__).parent
    evidence = directory.parent / 'evidence'
    source_names = ('fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py', 'f4_lattice_linear.py')
    source_hashes = {name: sha256((directory / name).read_bytes()).hexdigest() for name in source_names}
    saved = json.loads((evidence / 'fused-larger-endpoints.json').read_text())
    assert saved['source_sha256'] == source_hashes
    records = []
    for old in saved['cases']:
        started = time.monotonic()
        q = old['alphabet']
        cache = json.loads((evidence / f'fused-q{q}-matrix.json').read_text())
        matrix = [[{int(k): value for k, value in polynomial.items()} for polynomial in row] for row in cache]
        encoded = json.dumps(matrix, sort_keys=True, separators=(',', ':')).encode()
        assert sha256(encoded).hexdigest() == old['tensor_sha256']
        assert all(evaluate_one(polynomial) == int(i == j)
                   for i, row in enumerate(matrix) for j, polynomial in enumerate(row))
        assert max(k for row in matrix for polynomial in row for k in polynomial) == old['degree']
        if args.rebuild:
            source = construct(old['m'], old['theta'])
            assert source['invariant'] and source['endpoint_identity']
            assert source['matrix'] == matrix
            assert source['degree'] == old['degree']
            assert source['common_degree'] == old['cancelled_common_degree']
        records.append({'m': old['m'], 'theta': old['theta'], 'alphabet': q,
                        'degree': old['degree'], 'endpoint_entries': q ** 4,
                        'tensor_sha256': old['tensor_sha256'],
                        'full_invariant_subspace_rebuilt': args.rebuild,
                        'elapsed_seconds': time.monotonic() - started})
    result = {'scope': 'Exact cached seed/endpoint replay, with optional full invariant-subspace rebuild; no profile or score claim.',
              'source_sha256': source_hashes,
              'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'cases': records}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'alphabets': [row['alphabet'] for row in records],
                      'endpoint_entries': sum(row['endpoint_entries'] for row in records),
                      'full_invariant_subspace_rebuilt': args.rebuild}))


if __name__ == '__main__':
    main()
