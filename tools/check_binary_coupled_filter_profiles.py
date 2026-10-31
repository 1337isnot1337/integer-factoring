#!/usr/bin/env python3
"""Compare the binary coupled native backend with the independent jet backend."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path

from binary_coupled_filter_profiles import zero_profiles as native_profiles
from coupled_identity_mode_filters import monomial_filters, zero_profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    results = []
    for n, ramification, modes in ((2, 1, 2), (5, 1, 2), (5, 3, 2),
                                   (3, 1, 3), (5, 1, 3), (5, 3, 3)):
        states = list(product(range(2), repeat=modes))
        exponents = tuple(sum(state) - 2 * int(all(state)) for state in states)
        for polynomial in (False, True):
            alpha = monomial_filters(2, modes, exponents)
            if polynomial:
                alpha[states[-1]] = {0: 1, 1: 1, 2: 1}
            expected = zero_profiles(2, n, ramification, modes, alpha)
            actual = native_profiles(n, ramification, modes, alpha)
            assert [(row['charge'], row['orders']) for row in actual] == [
                (row['charge'], row['orders']) for row in expected]
            results.append({'n': n, 'ramification': ramification, 'modes': modes,
                            'polynomial_joint_weight': polynomial,
                            'matrix_size': 2**(2 * modes), 'charge_blocks': actual})
    output = {'status': 'Exact finite-field backend comparisons; no pressure claim.',
              'profile_comparisons': len(results), 'cases': results}
    if args.output:
        args.output.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({'profile_comparisons': len(results),
                      'individual_orders': sum(row['matrix_size'] for row in results)}))


if __name__ == '__main__':
    main()
