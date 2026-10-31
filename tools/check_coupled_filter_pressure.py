#!/usr/bin/env python3
"""Compare coupled-filter search pressure with full finite profile sums."""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import product
import json
from math import log, sqrt
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

import coupled_identity_extreme_filters as extreme
import coupled_identity_mode_filters as initial
from check_coupled_identity_mode_filters import expand_modes
from search_coupled_identity_filters import prepare_pressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    controls = []
    comparisons = 0
    maximum_error = 0.
    for prime, n, geometry in ((2, 2, 'initial'), (2, 4, 'terminal'), (2, 4, 'mixed'),
                               (3, 2, 'initial'), (3, 3, 'terminal'), (3, 3, 'mixed')):
        exponents = tuple(int(r == s) - 1 for r, s in product(range(prime), repeat=2))
        alpha = initial.monomial_filters(prime, 2, exponents)
        for ramification in (1, 3):
            module = initial if geometry == 'initial' else extreme
            options = {} if geometry == 'initial' else {'geometry': geometry}
            blocks = module.zero_profiles(prime, n, ramification, 2, alpha, **options)
            evaluate = prepare_pressure(prime, n, ramification, geometry, exponents, blocks)
            for charge in (Fraction(0), Fraction(2), Fraction(1, 2)):
                data = module.profile_data(prime, n, ramification, charge, 2, alpha,
                                           blocks=blocks, **options)
                zero = np.asarray([a + b for a in data['zero_active']
                                   for b in expand_modes(prime, data['zero_spectators'])], float)
                infinity = np.asarray([a + b for a in data['infinity_active']
                                       for b in expand_modes(prime, data['infinity_spectators'])], float)
                first = np.asarray(expand_modes(prime, data['paired_physical_first']), float)
                second = np.asarray(expand_modes(prime, data['paired_physical_second']), float)
                for scale in (1 / 16, .3, 1.):
                    ratio = sqrt(2)
                    physical_scale = scale / ramification
                    direct = (logsumexp(-physical_scale * zero)
                              + logsumexp(-physical_scale * infinity)
                              - ratio * (logsumexp(-physical_scale * first / ratio)
                                         + logsumexp(-physical_scale * second / ratio)))
                    candidate = evaluate(log(scale), float(charge))
                    error = abs(float(direct) - candidate)
                    assert error < 2e-10, (prime, n, geometry, ramification, charge, scale, direct, candidate)
                    maximum_error = max(maximum_error, error)
                    comparisons += 1
            controls.append({'prime': prime, 'n': n, 'geometry': geometry,
                             'ramification': ramification})
    result = {'status': 'Floating formula regression only; not a strict pressure certificate.',
              'full_expanded_profile_comparisons': comparisons,
              'maximum_absolute_difference': maximum_error, 'controls': controls}
    encoded = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded)


if __name__ == '__main__':
    main()
