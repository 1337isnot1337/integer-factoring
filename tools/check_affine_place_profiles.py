#!/usr/bin/env python3
"""Independent exact controls for affine reflection multipoint profiles."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from affine_place_profile_oracle import AffinePlaceOracle, bitplanes, polynomial, primitive_seed
from f16_identity_operator_profiles import multiply
from reflection_dual_rail_gaussian import tensor
from search_affine_coupled_reflection_frames import affine_seed, catalogue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    base = tensor(4)
    frames = catalogue(base)
    seeds = [('base', base), ('column_11', affine_seed(base, frames, [(56, 2), (58, 2)])),
             ('row_11', affine_seed(base, frames, [(18, 1), (66, 1)]))]
    rng = np.random.default_rng(312343)
    records, endpoint_entries, gradient_coordinates = [], 0, 0
    for name, seed in seeds:
        primitive, content = primitive_seed(seed)
        assert primitive['values'] == seed['values'] and content == {0: 1}
        # A nonconstant common factor with nontrivial endpoint scalar must
        # recover the exact original tensor after full normalization.
        factor = {0: 1, 1: 3}
        multiplied = dict(seed, values=[bitplanes(multiply(polynomial(value), factor))
                                        for value in seed['values']], degree=seed['degree'] + 1)
        # Divide by factor(1)=2 to keep the strict endpoint before the test.
        from reflection_laurent_physical_frames import scale
        multiplied['values'] = [scale(value, 3) for value in multiplied['values']]
        recovered, found = primitive_seed(multiplied)
        assert recovered['values'] == seed['values'] and max(found) == 1
        oracle = AffinePlaceOracle(seed)
        endpoint_entries += 4 ** 4
        for point in ((0, 1), (2, 1), (3, 1), None, (4, 1)):
            for cut in (0, 1):
                heights = [Fraction(int(n), 7) for n in rng.integers(-8, 9, 16)]
                scale_value = Fraction(5, 7)
                exact = oracle.profile(cut, point, heights, scale_value)
                native = oracle.log_profile_with_gradient(cut, point, np.asarray(heights, float), float(scale_value))
                assert np.allclose(np.asarray(exact, float), native[3], atol=2e-10, rtol=2e-12)
                unweighted = oracle.log_profile_with_gradient(cut, point, np.zeros(16), 1.0)[3]
                assert all(abs(float(order) - round(float(order))) < 1e-10 for order in unweighted)
                records.append({'seed': name, 'point': point, 'cut': cut,
                                'heights': list(map(str, heights)), 'scale': str(scale_value),
                                'weighted_orders': list(map(str, exact)),
                                'unweighted_orders': [round(float(order)) for order in unweighted]})
        # Non-tied real heights check the derivative used by the optimizer.
        heights = rng.normal(0, .31, 16)
        value, derivative, ds, _ = oracle.log_profile_with_gradient(0, (2, 1), heights, .73)
        epsilon = 1e-6
        for index in range(17):
            plus, minus = heights.copy(), heights.copy()
            sp = sm = .73
            if index < 16:
                plus[index] += epsilon
                minus[index] -= epsilon
                predicted = derivative[index]
            else:
                sp += epsilon
                sm -= epsilon
                predicted = ds
            observed = (oracle.log_profile_with_gradient(0, (2, 1), plus, sp)[0]
                        - oracle.log_profile_with_gradient(0, (2, 1), minus, sm)[0]) / (2 * epsilon)
            assert abs(observed - predicted) < 2e-6
            gradient_coordinates += 1
    sources = (Path(__file__).name, 'affine_place_profile_oracle.py', 'fused_physical_basis_oracle.py',
               'search_affine_coupled_reflection_frames.py', 'f16_weighted_local_oracle.py',
               'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp', 'native_rectangular_smith_f4.cpp')
    result = {'scope': 'Exact rational-height complete-profile controls; exploratory optimization is not certified.',
              'complete_profiles': len(records), 'weighted_channels': len(records) * 16,
              'endpoint_entries': endpoint_entries, 'gradient_coordinates': gradient_coordinates,
              'normalized_scalar_factor_controls': len(seeds), 'records': records,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in sources}}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('records', 'source_sha256')}))


if __name__ == '__main__':
    main()
