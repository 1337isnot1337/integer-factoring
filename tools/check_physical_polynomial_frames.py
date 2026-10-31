#!/usr/bin/env python3
"""Exact polynomial-frame inverses, endpoints, and independent local profiles."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from native_matrix_smith_f4 import MatrixOracle
from physical_polynomial_frames import paired_shear, polynomial, endpoint, apply_word
from physical_split_crossings import flatten
from reflection_dual_rail_gaussian import tensor
from search_reflection_physical_bases import packed_polynomial
from weighted_smith_oracle_f4 import WeightedSmithOracle


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    base = tensor(4)
    constant = dict(base, values=[(int(lo.bit_count() % 2), 0) for lo, hi in base['values']], degree=0)
    shapes = [[(0, 1), (1, 1)], [(1, 2)], [(0, 3), (3, 3)], [(1, 1), (4, 2)]]
    records = []
    counts = {'endpoint_entries': 0, 'exact_inverse_entries': 0, 'constant_crossing_entries': 0,
              'unweighted_finite_profile_channels': 0, 'weighted_profile_channels': 0}
    baseline = []
    for cut in (0, 1):
        matrix, _, _ = flatten(base, cut)
        oracle = WeightedSmithOracle([[packed_polynomial(p) for p in row] for row in matrix])
        baseline.append(oracle.profile([0] * 16, [0] * 16, Fraction(1))[0])
    for j, shape in enumerate(shapes):
        pair, target, source = j % 2, j, (j + 1) % 4
        f = polynomial(shape)
        changed = paired_shear(base, pair, target, source, f)
        restored = paired_shear(changed, pair, target, source, f)
        assert restored['values'] == base['values']
        fixed = paired_shear(constant, pair, target, source, f)
        assert fixed['values'] == constant['values']
        counts['endpoint_entries'] += endpoint(changed)
        counts['exact_inverse_entries'] += len(base['values'])
        counts['constant_crossing_entries'] += len(base['values'])
        controls = []
        for cut in (0, 1):
            matrix, _, _ = flatten(changed, cut)
            exact = WeightedSmithOracle([[packed_polynomial(p) for p in row] for row in matrix])
            orders = exact.profile([0] * 16, [0] * 16, Fraction(1))[0]
            assert orders == baseline[cut]
            counts['unweighted_finite_profile_channels'] += 16
            native = MatrixOracle(matrix, 16)
            rw = [Fraction((3 * i + 1) % 7 - 3, 2) for i in range(16)]
            cw = [Fraction((5 * i + 2) % 9 - 4, 2) for i in range(16)]
            for infinity in (False, True):
                independent = WeightedSmithOracle([[packed_polynomial(p) for p in row] for row in matrix],
                                                  infinity=infinity)
                expected = independent.profile(rw, cw, Fraction(3, 2))[0]
                got = native.profile(rw, cw, scale=1.5, infinity=infinity)['orders']
                assert sorted(got) == sorted(map(float, expected))
                counts['weighted_profile_channels'] += 16
                controls.append({'cut': cut, 'infinity': infinity,
                                 'orders': [[p.numerator, p.denominator] for p in expected]})
        records.append({'pair': pair, 'target': target, 'source': source, 'terms': shape,
                        'degree_before': base['degree'], 'degree_after': changed['degree'],
                        'complete_profiles': controls})
    word = [{'pair': 0, 'target': 0, 'source': 1, 'terms': shapes[0]},
            {'pair': 1, 'target': 2, 'source': 3, 'terms': shapes[2]}]
    assert apply_word(base, word)['values'] == apply_word(base, list(reversed(word)))['values']
    names = ('check_physical_polynomial_frames.py', 'physical_polynomial_frames.py',
             'weighted_smith_oracle_f4.py', 'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    result = {'status': 'pass', 'scope': 'exact polynomial-frame identities and complete local controls',
              'counts': counts, 'disjoint_pair_commutation_entries': 256, 'records': records,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass', **counts}))


if __name__ == '__main__':
    main()
