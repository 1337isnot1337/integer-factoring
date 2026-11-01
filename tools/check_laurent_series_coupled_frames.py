#!/usr/bin/env python3
"""Exact independent controls for coupled matrix Laurent-series frames."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import random

from check_reflection_laurent_frames import evaluate
from f16_identity_operator_profiles import MUL
from fused_six_vertex_crossings import power
from laurent_series_coupled_frames import construct, from_scalar_terms
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, digest
from square_zero_coupled_reflection_frames import matrix_terms, transpose


def evaluate_action(values, q, leg, matrix):
    stride = q ** leg
    output = [0] * len(values)
    for index in range(len(values)):
        target = index // stride % q
        for source, coefficient in enumerate(matrix[target]):
            if coefficient:
                output[index] ^= MUL[coefficient][values[index + (source - target) * stride]]
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(312371)
    records, coefficients_checked, involution_entries, evaluated_entries = [], 0, 0, 0
    for interval, samples in ((4, 12), (8, 3)):
        base = tensor(interval)
        q = base['q_x']
        bell, factors = bell_normalize(base)
        for sample in range(samples):
            coefficients = {}
            for pole in range(1, 3 if interval == 4 else 5):
                matrix = [[0] * q for _ in range(q)]
                for i in range(q // 2):
                    for j in range(q // 2, q):
                        if pole <= 2 or (i.bit_count() + j.bit_count()) % 2 == 0:
                            matrix[i][j] = rng.randrange(4)
                coefficients[pole] = matrix
            changed = construct(base, coefficients, restore=False, bell_data=(bell, factors))
            restored = construct(changed, coefficients, restore=False, bell_data=(changed, factors))
            assert restored['values'] == bell['values']
            involution_entries += q ** 4
            coefficients_checked += q ** 4
            if interval == 4 and sample < 4:
                for point in (2, 3, 4, 8):
                    matrix = [[int(i == j) for j in range(q)] for i in range(q)]
                    for pole, coefficient in coefficients.items():
                        scalar = power(point, -pole)
                        for i in range(q):
                            for j in range(q):
                                matrix[i][j] ^= MUL[scalar][coefficient[i][j]]
                    values = [evaluate(value, point) for value in bell['values']]
                    for leg in range(4):
                        values = evaluate_action(values, q, leg, transpose(matrix) if leg % 2 else matrix)
                    assert values == [evaluate(value, point) for value in changed['values']]
                    evaluated_entries += q ** 4
            records.append({'interval': interval, 'sample': sample, 'coefficients': coefficients,
                            'degree': changed['degree'], 'tensor_sha256': digest(changed['values'])})
        # A second, independent expansion in the scalar f checks all f^r
        # cross terms against the sequential four-leg construction.
        generator = [[0] * q for _ in range(q)]
        generator[0][q // 2] = 1
        terms = matrix_terms(bell, generator)
        for scalars in ({1: 1, 2: 1}, {1: 2, 2: 3}, {1: 3, 2: 2}):
            matrices = {pole: [[MUL[c][x] for x in row] for row in generator] for pole, c in scalars.items()}
            direct = construct(base, matrices, restore=False, bell_data=(bell, factors))
            expanded = from_scalar_terms(bell, terms, scalars)
            assert direct['values'] == expanded['values']
            coefficients_checked += q ** 4
    sources = (Path(__file__).name, 'laurent_series_coupled_frames.py', 'square_zero_coupled_reflection_frames.py',
               'reflection_laurent_physical_frames.py', 'reflection_dual_rail_gaussian.py')
    result = {'scope': 'Exact polynomiality, degree, endpoint, involution, and independent evaluation controls; no pressure claim.',
              'random_series_frames': len(records), 'coefficient_comparison_entries': coefficients_checked,
              'involution_entries': involution_entries, 'independent_F16_evaluations': evaluated_entries,
              'records': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in sources}}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('records', 'source_sha256')}))


if __name__ == '__main__':
    main()
