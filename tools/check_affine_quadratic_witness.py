#!/usr/bin/env python3
"""Exact full determinant and quadratic-place profiles of one affine path."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from affine_place_profile_oracle import AffinePlaceOracle, primitive_seed
from check_fused_multipoint_profiles import field_rank, local_polynomial, toeplitz_rank
from check_reflection_laurent_frames import cut_matrix
from f16_identity_operator_profiles import MUL, multiply
from f16_weighted_local_oracle import F16WeightedLocalOracle
from fused_six_vertex_crossings import power
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import digest
from search_affine_coupled_reflection_frames import affine_seed, catalogue


def polynomial_power(value, exponent):
    result = {0: 1}
    while exponent:
        if exponent & 1:
            result = multiply(result, value)
        value = multiply(value, value)
        exponent //= 2
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    base = tensor(4)
    frames = catalogue(base)
    word = [[19, 3], [89, 2]]
    seed = affine_seed(base, frames, word)
    primitive, common = primitive_seed(seed)
    assert common == {0: 1} and primitive['values'] == seed['values']
    expected_hash = 'dd9490bc487c838a1f29363522a9f56ab06f83becdbc4124a85265198da2ef9c'
    assert digest(seed['values']) == expected_hash
    cuts = [cut_matrix(seed, cut) for cut in (0, 1)]
    source_cuts = [cut_matrix(base, cut) for cut in (0, 1)]
    factors = {a: {0: MUL[a][power(a, 4)], 1: a ^ power(a, 4), 2: 1}
               for a in (4, 6, 9, 12)}
    assert all(c in range(4) for f in factors.values() for c in f.values())
    determinants = []
    for cut in (0, 1):
        oracle = F16WeightedLocalOracle(cuts[cut], 24)
        result = oracle.profile([0] * 16, [0] * 16, (0, 1))
        roots = (6, 9) if cut == 0 else (4, 12)
        expected = {160: 1} if cut == 0 else {0: 1}
        for root in roots:
            expected = multiply(expected, polynomial_power(factors[root], 24))
        assert result['determinant'] == expected
        at_infinity = oracle.profile([0] * 16, [0] * 16, None)
        determinants.append({'cut': cut, 'degree': max(expected),
                             'zero_order': 160 if cut == 0 else 0,
                             'complete_zero_orders': list(result['orders']),
                             'complete_homogeneous_infinity_orders': list(at_infinity['orders']),
                             'quadratic_factors': [{'representative': root, 'polynomial': factors[root],
                                                    'multiplicity': 24} for root in roots],
                             'determinant': expected})
    records, toeplitz_controls = [], 0
    for point, cut in ((4, 1), (6, 0), (9, 0), (12, 1)):
        for source_cut in source_cuts:
            evaluated = [[local_polynomial(value, point, 24).get(0, 0) for value in row]
                         for row in source_cut]
            assert field_rank(evaluated) == 16
        local = [[local_polynomial(value, point, 24) for value in row] for row in cuts[cut]]
        orders = F16WeightedLocalOracle(local, 24).profile([0] * 16, [0] * 16, (0, 1))['orders']
        assert list(orders) == [0] * 11 + [4] * 4 + [8]
        ranks = []
        for length in (1, 4, 5, 8, 9):
            observed = toeplitz_rank(local, length)
            expected = sum(max(length - order, 0) for order in orders)
            assert observed == expected
            ranks.append([length, observed])
            toeplitz_controls += 1
        records.append({'conjugates': [point, power(point, 4)], 'cut': cut,
                        'complete_orders': list(orders), 'toeplitz_ranks': ranks})
    # Every nonzero determinant factor is now explicit.  Thus the four
    # recorded conjugate pairs are the entire nonzero finite singular set.
    source_names = (Path(__file__).name, 'affine_place_profile_oracle.py',
                    'search_affine_coupled_reflection_frames.py', 'f16_weighted_local_oracle.py',
                    'check_fused_multipoint_profiles.py', 'fused_six_vertex_crossings.py')
    result = {'scope': 'Exact affine crossing and complete singular-place catalogue; no optimized pressure or factoring claim.',
              'word': word, 'frame_metadata': [{k: v for k, v in frames[index].items() if k != 'delta'}
                                              for index, coefficient in word],
              'tensor_sha256': expected_hash, 'degree': seed['degree'], 'common_factor': common,
              'determinants': determinants, 'quadratic_place_profiles': records,
              'independent_block_Toeplitz_ranks': toeplitz_controls,
              'independent_source_place_ranks': 8,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in source_names}}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'degree': seed['degree'], 'tensor_sha256': expected_hash,
                      'nonzero_finite_singular_points': 8, 'complete_profile_representatives': len(records),
                      'independent_block_Toeplitz_ranks': toeplitz_controls}))


if __name__ == '__main__':
    main()
