#!/usr/bin/env python3
"""Exact independent checks of the explicit degree16 affine crossing."""
from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path

from check_dense_coupled_reflection_frames import coefficients, apply_matrix, multiply
from close_independent_paired_affine_span import operators, image
from dense_coupled_reflection_frames import digest
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize
from solve_companion_affine_degree import full_vector, scaled, unpack


def paired_traces(values):
    first, second = [], []
    for y in range(4):
        for u in range(4):
            a = b = 0
            for x in range(4):
                low, high = values[x + 4 * y + 16 * u + 64 * x]
                a, b = a ^ low, b ^ high
            first.append((a, b))
    for x in range(4):
        for v in range(4):
            a = b = 0
            for y in range(4):
                low, high = values[x + 4 * y + 16 * y + 64 * v]
                a, b = a ^ low, b ^ high
            second.append((a, b))
    return first, second


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--crossing', type=Path,
                        default=Path('evidence/degree16-reflection-crossing.json'))
    parser.add_argument('--basis', type=Path,
                        default=Path('evidence/independent-paired-affine-degree-basis.json.gz'))
    parser.add_argument('--record', type=Path,
                        default=Path('evidence/independent-paired-affine-degree-prefix.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    data = json.loads(args.crossing.read_text())
    record = json.loads(args.record.read_text())
    basis_raw = args.basis.read_bytes()
    basis = json.loads(gzip.decompress(basis_raw))
    assert sha256(args.record.read_bytes()).hexdigest() == data['source_record_sha256']
    assert sha256(basis_raw).hexdigest() == data['source_basis_sha256']
    assert data['source_basis_sha256'] == record['output_basis_sha256']
    values = list(map(tuple, data['values']))
    assert len(values) == 256 and all(a >= 0 and b >= 0 for a, b in values)
    assert digest(values) == data['seed_sha256']
    assert max((a | b).bit_length() - 1 for a, b in values) == 16
    for index, (a, b) in enumerate(values):
        x, y, u, v = (index // 4 ** leg % 4 for leg in range(4))
        assert (a.bit_count() % 2, b.bit_count() % 2) == (int(x == v and y == u), 0)
    trace = 1 | (1 << 16)
    for matrix in paired_traces(values):
        assert matrix == [(trace if i == j else 0, 0) for i in range(4) for j in range(4)]
    base = tensor(4)
    a, b = full_vector(base['values'])
    for index, scalar in data['basis_witness']:
        vector = tuple(int(value, 16) for value in basis['basis'][index]['bitplanes'])
        x, y = scaled(vector, scalar)
        a, b = a ^ x, b ^ y
    assert unpack((a, b)) == values
    assert len(data['basis_witness']) == data['basis_witness_terms'] == 1168
    # These trace checks verify the claimed affine constraints on every
    # retained basis vector, independently of the search's degree solve.
    for row in basis['basis']:
        vector = tuple(int(value, 16) for value in row['bitplanes'])
        entry_values = unpack(vector)
        assert all(not (x.bit_count() % 2 or y.bit_count() % 2) for x, y in entry_values)
        assert all(not any(a or b for a, b in matrix) for matrix in paired_traces(entry_values))
    # Independent exponent-dictionary actions check the fast packed
    # transvection operators in both shift directions and both physical pairs.
    packed = full_vector(base['values'])
    ops = operators()
    checked = 0
    for pair, (leg, dual) in enumerate(((0, 3), (1, 2))):
        for target, source in ((0, 1), (3, 2)):
            linear = next(o for o in ops if (o['pair'], o['target'], o['source'], o['kind'])
                          == (pair, target, source, 'linear'))
            quadratic = next(o for o in ops if (o['pair'], o['target'], o['source'], o['kind'])
                             == (pair, target, source, 'quadratic'))
            for scalar in (1, 2):
                matrix = [[0] * 4 for _ in range(4)]
                matrix[target][source] = scalar
                reference = [coefficients(value) for value in base['values']]
                reference = apply_matrix(reference, matrix, leg, 4, identity=True)
                reference = apply_matrix(reference, list(map(list, zip(*matrix))), dual, 4, identity=True)
                la, lb = scaled(image(packed, linear), scalar)
                qa, qb = scaled(image(packed, quadratic), multiply(scalar, scalar))
                actual = unpack((packed[0] ^ la ^ qa, packed[1] ^ lb ^ qb))
                assert reference == [coefficients(value) for value in actual]
                checked += 256
    bell, _ = bell_normalize(base)
    double = [0, 0]
    for x in range(4):
        for u in range(4):
            a, b = bell['values'][x + 4 * x + 16 * u + 64 * u]
            double[0] ^= a
            double[1] ^= b
    assert double == [sum(1 << power for power in (8, 12, 20, 24)), 0]
    result = {
        'status': 'passed',
        'scope': 'exact strict crossing, degree16, sharp paired-trace floor, basis reconstruction and independent physical-action controls; no pressure certificate',
        'crossing_sha256': sha256(args.crossing.read_bytes()).hexdigest(),
        'tensor_sha256': data['seed_sha256'],
        'exact_degree': 16,
        'endpoint_entries': 256,
        'paired_trace_polynomial_entries': 32,
        'basis_vectors_with_zero_endpoint_and_paired_traces': len(basis['basis']),
        'affine_witness_terms': len(data['basis_witness']),
        'independent_physical_action_entries': checked,
        'original_coupled_bell_double_trace_bitplanes': double,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, 'solve_companion_affine_degree.py',
                                       'close_independent_paired_affine_span.py',
                                       'check_dense_coupled_reflection_frames.py')},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'source_sha256'}))


if __name__ == '__main__':
    main()
