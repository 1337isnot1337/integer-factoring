#!/usr/bin/env python3
"""Exact controls for affine operator bit-plane compilation."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
import random

from compiled_f4_identity_operators import compile_family, combine_kernel
from coupled_identity_operator_filters import identity_kernel, states
from f4_identity_operator_filters import zero_profiles, infinity_profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    controls = []
    for n, P, modes, geometry in ((3, 1, 2, 'initial'), (8, 3, 2, 'mixed'),
                                  (4, 1, 3, 'terminal'), (8, 1, 3, 'initial')):
        alphabet = states(2, modes)
        sector = [state for state in alphabet if sum(state) == 1]
        perturbations = [{(a, b): {0: 1, -1-(i%2): 1}}
                         for i, (a, b) in enumerate(product(sector, repeat=2))]
        constant = identity_kernel(2, modes)
        compiled = compile_family(n, P, modes, constant, perturbations, geometry)
        rng = random.Random(100*n+10*P+modes)
        for sample in range(12):
            colors = tuple(rng.randrange(4) for perturbation in perturbations)
            kernel = combine_kernel(constant, perturbations, colors)
            actual = compiled(colors)
            expected = (zero_profiles(n, P, modes, kernel, geometry),
                        infinity_profiles(n, P, modes, kernel, geometry))
            assert actual == expected
            controls.append({'n': n, 'ramification': P, 'modes': modes,
                             'geometry': geometry, 'colors': colors})
    constant = identity_kernel(2, 2)
    perturbations = [{((1, 0), (0, 1)): {0: 1, 600: 1}}]
    compiled = compile_family(2, 1, 2, constant, perturbations, 'initial')
    kernel = combine_kernel(constant, perturbations, (2,))
    assert compiled((2,)) == (zero_profiles(2, 1, 2, kernel), infinity_profiles(2, 1, 2, kernel))
    result = {'status': 'Exact complete affine operator compilation controls; no pressure claim.',
              'two_cut_controls': controls, 'late_second_plane_degree': 600}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'complete_two_cut_comparisons': len(controls)+1,
                      'late_second_plane_degree': 600}))


if __name__ == '__main__':
    main()
