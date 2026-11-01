#!/usr/bin/env python3
"""Independent literal controls and finite root-of-unity fusion certificates."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

from compressed_fused_six_vertex import construct_compressed, rectangular_row
from f16_identity_operator_profiles import MUL, evaluate_one, multiply
from fused_six_vertex_crossings import construct, sym, gate, power, pscale
from check_fused_multipoint_profiles import field_rank


def digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def literal_row_control(m, theta, shift):
    compressed = rectangular_row(m, theta, shift)
    checked = 0
    for a, b in product(range(2), range(m + 1)):
        vector = {a | (word << 1): {0: coefficient}
                  for word, coefficient in sym(m, theta, b).items()}
        for j in range(m):
            vector = gate(vector, j, theta, power(theta, 2 * (shift - j)))
        reconstructed = {}
        for c, d in product(range(m + 1), range(2)):
            coefficient = compressed['matrix'].get((c, d, a, b), {})
            raw = multiply(coefficient, compressed['common_factor'])
            for word, scalar in sym(m, theta, c).items():
                value = pscale(raw, scalar)
                if value:
                    reconstructed[word | (d << m)] = value
        assert reconstructed == vector
        checked += 2 ** (m + 1)
    return checked


def endpoint_rank(endpoint, q):
    total = 0
    for charge in range(2 * q - 1):
        indices = [q * a + b for a, b in product(range(q), repeat=2) if a + b == charge]
        total += field_rank([[endpoint[i][j] for j in indices] for i in indices])
    return total


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-m', type=int, default=18)
    args = parser.parse_args()
    if args.max_m < 7:
        raise ValueError('The published control set includes m=7.')
    directory = Path(__file__).parent
    evidence = directory.parent / 'evidence'
    row_controls = []
    for m, theta in [(m, 4) for m in range(1, 8)] + [(2, 2), (3, 2), (4, 2)]:
        for shift in sorted({0, m - 1}):
            row_controls.append({'m': m, 'theta': theta, 'shift': shift,
                                 'literal_entries': literal_row_control(m, theta, shift)})
    literal_controls = []
    for m, theta in [(1, 4), (2, 2), (2, 4), (3, 2), (3, 4), (4, 2), (4, 4), (5, 4), (6, 4), (7, 4)]:
        compressed = construct_compressed(m, theta)
        if m <= 5:
            original = construct(m, theta)
            assert original['invariant']
            assert original['matrix'] == compressed['matrix']
            assert original['degree'] == compressed['degree']
            assert original['common_degree'] == compressed['common_degree']
            provenance = 'Repeated full 2m-bit construction in this run.'
        else:
            cached = json.loads((evidence / f'fused-q{m+1}-matrix.json').read_text())
            original_matrix = [[{int(k): v for k, v in p.items()} for p in row] for row in cached]
            assert original_matrix == compressed['matrix']
            provenance = 'Complete comparison with the separately certified literal matrix cache.'
        literal_controls.append({'m': m, 'theta': theta, 'alphabet': m + 1,
                                 'matrix_entries': (m + 1) ** 4,
                                 'matrix_sha256': digest(compressed['matrix']),
                                 'provenance': provenance})
    cases = []
    for m, theta in [(2, 2), (3, 2), (4, 2)] + [(m, 4) for m in range(1, args.max_m + 1)]:
        started = time.monotonic()
        source = construct_compressed(m, theta, normalize_joint=True)
        q = m + 1
        assert source['invariant'] and source['endpoint_identity']
        assert all(evaluate_one(value) == int(i == j)
                   for i, row in enumerate(source['matrix']) for j, value in enumerate(row))
        rank = endpoint_rank(source['endpoint_before_normalization'], q)
        assert rank == q * q
        record = {key: value for key, value in source.items()
                  if key not in ('matrix', 'endpoint_before_normalization', 'endpoint_inverse')}
        record.update({'raw_endpoint_rank': rank,
                       'raw_endpoint_sha256': digest(source['endpoint_before_normalization']),
                       'endpoint_inverse_sha256': digest(source['endpoint_inverse']),
                       'normalized_matrix_sha256': digest(source['matrix']),
                       'endpoint_entries': q ** 4,
                       'elapsed_seconds': time.monotonic() - started})
        cases.append(record)
        print(json.dumps({'m': m, 'theta': theta, 'q': q, 'degree': source['degree'],
                          'joint_normalized': source['joint_endpoint_normalization'],
                          'seconds': record['elapsed_seconds']}), flush=True)
    sources = ('compressed_fused_six_vertex.py', Path(__file__).name,
               'fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py',
               'f4_lattice_linear.py', 'check_fused_multipoint_profiles.py')
    result = {
        'scope': 'Exact finite compressed fusion and explicit joint endpoint normalization; no improved profile or factoring score claimed.',
        'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in sources},
        'literal_row_controls': row_controls,
        'literal_full_matrix_controls': literal_controls,
        'cases': cases,
        'totals': {'literal_row_entries': sum(row['literal_entries'] for row in row_controls),
                   'literal_full_matrix_entries': sum(row['matrix_entries'] for row in literal_controls),
                   'rectangular_symmetry_checks': sum(row['rectangular_symmetry_checks'] for row in cases),
                   'output_symmetry_checks': sum(row['output_symmetry_checks'] for row in cases),
                   'strict_endpoint_entries': sum(row['endpoint_entries'] for row in cases)}}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['totals']))


if __name__ == '__main__':
    main()
