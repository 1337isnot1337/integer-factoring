#!/usr/bin/env python3
"""Complete singular-place certificates for two primitive finite-fiber cores.

The exact fraction-free determinant is checked independently by enough
evaluations in F256 to prove a polynomial identity. Every local Smith order
is checked again by truncated block-Toeplitz ranks.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from check_frobenius_rational_trace import field_of_degree
from check_fused_multipoint_profiles import toeplitz_rank
from check_reflection_laurent_frames import cut_matrix
from f16_weighted_local_oracle import F16WeightedLocalOracle
from reflection_laurent_physical_frames import digest


EXPECTED = {
    'infinity4_primitive_degree20_deflated': {
        'degree': 10,
        'tensor_sha256': 'cb796f650714d5b8f4952c0b87bb0097337582b66e9674bc16e80bae22e580ed',
        'row_determinant_order': 112,
    },
    'infinity4_primitive_degree22_deflated': {
        'degree': 11,
        'tensor_sha256': '30fa0b9badc2f0e7853faab8686db7d902a835f15754cfa9d1c67f2af62a2d47',
        'row_determinant_order': 120,
    },
}


def finite_field_tables():
    field = field_of_degree(8)
    multiplication = [[field.mul(a, b) for b in range(256)]
                      for a in range(256)]
    inverse = [0] + [field.power(a, 254) for a in range(1, 256)]
    omega = next(a for a in range(2, 256)
                 if multiplication[a][a] ^ a ^ 1 == 0)
    embedding = (0, 1, omega, omega ^ 1)
    return field, multiplication, inverse, embedding


def numeric_determinant(matrix, multiplication, inverse):
    """Ordinary field elimination, independent of the polynomial oracle."""
    work = [row[:] for row in matrix]
    result = 1
    for j in range(len(work)):
        pivot = next((i for i in range(j, len(work)) if work[i][j]), None)
        if pivot is None:
            return 0
        work[j], work[pivot] = work[pivot], work[j]
        value = work[j][j]
        result = multiplication[result][value]
        scalar = inverse[value]
        work[j] = [multiplication[scalar][v] for v in work[j]]
        for i in range(j + 1, len(work)):
            scalar = work[i][j]
            if scalar:
                work[i] = [a ^ multiplication[scalar][b]
                           for a, b in zip(work[i], work[j])]
    return result


def evaluate_matrix(matrix, point, multiplication, embedding, degree):
    powers = [1]
    for _ in range(degree):
        powers.append(multiplication[powers[-1]][point])
    result = []
    for row in matrix:
        output = []
        for polynomial in row:
            value = 0
            for exponent, coefficient in polynomial.items():
                value ^= multiplication[embedding[coefficient]][powers[exponent]]
            output.append(value)
        result.append(output)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=Path('evidence/deflated-finite-fiber-input.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    input_bytes = args.input.read_bytes()
    inputs = json.loads(input_bytes)
    field, multiplication, inverse, embedding = finite_field_tables()
    cases = []
    determinant_evaluations = toeplitz_controls = endpoint_entries = 0
    for entry in inputs['cases']:
        expected = EXPECTED[entry['label']]
        seed = entry['seed']
        degree = seed['degree']
        assert degree == expected['degree']
        assert digest(seed['values']) == expected['tensor_sha256']
        assert seed['q_x'] == seed['q_u'] == 4
        # Direct endpoint reconstruction avoids relying on a stored flag.
        for index, (low, high) in enumerate(seed['values']):
            x, y, u, v = ((index // 4 ** leg) % 4 for leg in range(4))
            assert (low.bit_count() & 1) == int(x == v and y == u)
            assert high.bit_count() % 2 == 0
            endpoint_entries += 1
        records = []
        for cut in (0, 1):
            matrix = cut_matrix(seed, cut)
            oracle = F16WeightedLocalOracle(matrix, degree)
            zero = oracle.profile([0] * 16, [0] * 16, (0, 1))
            infinity = oracle.profile([0] * 16, [0] * 16, None)
            exponent = expected['row_determinant_order'] if cut == 0 else 0
            assert zero['determinant'] == {exponent: 1}
            assert infinity['determinant'] == {exponent: 1}
            # A 16 by 16 matrix of degree at most d has determinant degree
            # at most 16d. Agreement at 16d+1 distinct field points proves
            # the claimed identity over F256, and hence over embedded F4.
            count = 16 * degree + 1
            assert count <= 256
            for point in range(count):
                evaluated = evaluate_matrix(matrix, point, multiplication,
                                            embedding, degree)
                observed = numeric_determinant(evaluated, multiplication, inverse)
                assert observed == field.power(point, exponent)
                determinant_evaluations += 1
            local_records = []
            for name, local, profile in (
                ('zero', matrix, zero['orders']),
                ('homogeneous_infinity',
                 [[{degree - k: c for k, c in value.items()} for value in row]
                  for row in matrix], infinity['orders']),
            ):
                checks = []
                # All lengths through the maximum order plus one determine
                # the entire multiset, including the zero channels.
                for length in range(1, max(profile) + 2):
                    actual = toeplitz_rank(local, length)
                    target = sum(max(length - order, 0) for order in profile)
                    assert actual == target
                    checks.append([length, actual])
                    toeplitz_controls += 1
                local_records.append({'place': name, 'orders': list(profile),
                                      'independent_toeplitz_ranks': checks})
            assert sum(zero['orders']) == exponent
            assert sum(infinity['orders']) == 16 * degree - exponent
            records.append({'cut': cut, 'determinant': {str(exponent): 1},
                            'determinant_evaluation_count': count,
                            'local_profiles': local_records,
                            'other_finite_singular_places': []})
        cases.append({'label': entry['label'], 'degree': degree,
                      'tensor_sha256': expected['tensor_sha256'],
                      'parent_tensor_sha256': entry['parent_tensor_sha256'],
                      'parent_degree': entry['parent_degree'],
                      'parameter_ramification': entry['ramification'],
                      'cuts': records})
    names = (Path(__file__).name, 'check_frobenius_rational_trace.py',
             'check_carry_fourier_elimination.py',
             'check_fused_multipoint_profiles.py', 'check_reflection_laurent_frames.py',
             'f16_weighted_local_oracle.py', 'check_fused_closed_profiles.py',
             'f16_identity_operator_profiles.py', 'fused_six_vertex_crossings.py',
             'reflection_laurent_physical_frames.py')
    result = {
        'scope': 'Complete singular-place catalogues of two explicit primitive cores; no pressure or factoring improvement.',
        'input_path': str(args.input), 'input_sha256': sha256(input_bytes).hexdigest(),
        'independent_field': {'degree': 8, 'irreducible_polynomial': field.polynomial,
                              'F4_embedding': list(embedding)},
        'cases': cases,
        'counts': {'independent_determinant_evaluations': determinant_evaluations,
                   'independent_block_Toeplitz_ranks': toeplitz_controls,
                   'strict_endpoint_entries': endpoint_entries},
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
        'elapsed_seconds': time.monotonic() - started, 'passed': True,
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('scope', 'counts', 'elapsed_seconds', 'passed')}))


if __name__ == '__main__':
    main()
