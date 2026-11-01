#!/usr/bin/env python3
"""Exact coefficient and native pressure controls for deflated fiber seeds."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from asymmetric_crossing_pressure import AsymmetricPressure
from deflate_crossing_parameter import deflate
from reflection_laurent_physical_frames import digest
from search_affine_coupled_reflection_frames import balanced


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    original = args.input.read_bytes()
    data = json.loads(original)
    cases = []
    for item in data['cases']:
        seed = item['seed']
        assert digest(seed['values']) == item['tensor_sha256']
        reduced = deflate(seed)
        degree = reduced['degree']
        ramification = reduced['parameter_ramification_removed']
        before = np.asarray(item['initial_params'], dtype=float)
        after = before.copy()
        after[-1] *= ramification
        cost1 = balanced(AsymmetricPressure(seed, 1), before)
        cost2 = balanced(AsymmetricPressure(reduced, 1), after)
        assert abs(cost1[0]-cost2[0]) < 2e-10
        cases.append({'label': item['label'] + '_deflated', 'seed': reduced,
                      'tensor_sha256': digest(reduced['values']), 'parent_tensor_sha256': item['tensor_sha256'],
                      'parent_degree': seed['degree'], 'degree': degree, 'ramification': ramification,
                      'initial_params': after.tolist(), 'initial_balanced_rho': cost2[0],
                      'pressure_identity_error': abs(cost1[0]-cost2[0]),
                      'initial_types': item['initial_types']})
    names = (Path(__file__).name, 'deflate_crossing_parameter.py', 'asymmetric_crossing_pressure.py',
             'search_affine_coupled_reflection_frames.py', 'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    result = {'scope': __doc__, 'parent_input_sha256': sha256(original).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'passed': True, 'cases': cases}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps([{key: item[key] for key in ('label', 'parent_degree', 'degree', 'ramification',
                                                 'initial_balanced_rho', 'pressure_identity_error')} for item in cases]))


if __name__ == '__main__':
    main()
