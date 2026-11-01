#!/usr/bin/env python3
"""Independent weighted-profile and gradient controls for physical bases."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import random

import numpy as np

from fused_physical_basis_oracle import FusedPhysicalBasisOracle, identity, quantum_pascal, vandermonde
from fused_six_vertex_crossings import power


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    generator = random.Random(311927)
    records = []
    for m, theta in ((2, 2), (3, 4)):
        q = m + 1
        bases = [('quantum_pascal_evaluation', quantum_pascal(q, power(theta, 2)), vandermonde(range(q))),
                 ('evaluation_common', vandermonde(range(q)), vandermonde(range(q))),
                 ('pascal_common', quantum_pascal(q, 1), quantum_pascal(q, 1))]
        for name, first, second in bases:
            oracle = FusedPhysicalBasisOracle(m, theta, first, second)
            comparisons = 0
            largest_profile_error = 0.0
            for point in ((0, 1), (2, 1), None):
                for cut in (0, 1):
                    heights = [Fraction(generator.randrange(-9, 10), 7) for _ in range(4 * q)]
                    scale = Fraction(5, 3)
                    _, _, _, native = oracle.log_profile_with_gradient(cut, point, heights, scale)
                    exact = np.asarray([float(value) for value in oracle.profile(cut, point, heights, scale)])
                    error = float(np.max(np.abs(native - exact)))
                    assert error < 2e-8
                    largest_profile_error = max(largest_profile_error, error)
                    comparisons += 1
            rng = np.random.default_rng(1100 + q)
            heights = rng.normal(0, 0.37, 4 * q)
            scale, epsilon = 0.7, 1e-6
            _, gradient, scale_gradient, _ = oracle.log_profile_with_gradient(1, (0, 1), heights, scale)
            errors = []
            for j in range(4 * q):
                plus, minus = heights.copy(), heights.copy()
                plus[j] += epsilon
                minus[j] -= epsilon
                numerical = (oracle.log_profile_with_gradient(1, (0, 1), plus, scale)[0]
                             - oracle.log_profile_with_gradient(1, (0, 1), minus, scale)[0]) / (2 * epsilon)
                errors.append(abs(float(numerical - gradient[j])))
            numerical = (oracle.log_profile_with_gradient(1, (0, 1), heights, scale + epsilon)[0]
                         - oracle.log_profile_with_gradient(1, (0, 1), heights, scale - epsilon)[0]) / (2 * epsilon)
            errors.append(abs(float(numerical - scale_gradient)))
            assert max(errors) < 2e-6
            records.append({'m': m, 'theta': theta, 'alphabet': q, 'basis_name': name,
                            'first_basis': first, 'second_basis': second,
                            'complete_exact_profile_comparisons': comparisons,
                            'largest_profile_error': largest_profile_error,
                            'gradient_coordinates': len(errors), 'largest_gradient_error': max(errors),
                            'endpoint_entries': q ** 4})
    directory = Path(__file__).parent
    names = (Path(__file__).name, 'fused_physical_basis_oracle.py', 'compressed_fused_six_vertex.py',
             'f16_weighted_local_oracle.py', 'fused_six_vertex_crossings.py',
             'f16_identity_operator_profiles.py', 'f4_lattice_linear.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp',
             'native_rectangular_smith_f4.cpp', 'native_weighted_smith_f4.py')
    result = {'scope': 'Exact rational F16 profiles versus native restriction of scalars after complete physical basis changes; floating gradients are controls, not an asymptotic certificate.',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in names},
              'cases': records,
              'totals': {'complete_exact_profile_comparisons': sum(row['complete_exact_profile_comparisons'] for row in records),
                         'gradient_coordinates': sum(row['gradient_coordinates'] for row in records),
                         'endpoint_entries': sum(row['endpoint_entries'] for row in records)}}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['totals']))


if __name__ == '__main__':
    main()
