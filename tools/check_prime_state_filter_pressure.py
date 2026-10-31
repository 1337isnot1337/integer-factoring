#!/usr/bin/env python3
"""Numerical controls for the ground-separated prime-filter pressure helper."""
from __future__ import annotations

import argparse
from itertools import product
import json
from math import sqrt
from pathlib import Path
import numpy as np
from scipy.special import logsumexp

from prime_state_filter_profiles import active_profiles
from prime_state_filter_pressure import filter_pressure


def expand(orders, prime):
    result = [0]
    for order in orders:
        result = [value + digit * order for value in result for digit in range(prime)]
    return result


def direct_full_pressure(data, prime, ramification, scale):
    zero = np.array([a + b for a in data['zero_active']
                     for b in expand(data['zero_spectators'], prime)])
    infinity = np.array([a + b for a in data['infinity_active']
                         for b in expand(data['infinity_spectators'], prime)])
    physical = data['paired_physical_first'] + data['paired_physical_second']
    ratio = sqrt(2)
    return float(logsumexp(-zero * scale / ramification)
                 + logsumexp(-infinity * scale / ramification)
                 - ratio * sum(logsumexp(-np.arange(prime) * value * scale / ramification / ratio)
                               for value in physical))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases, worst_direct, worst_clearing = 0, 0.0, 0.0
    for prime, n, ramification, charge, scale in product(
            (3, 5), (1, 2), (1, 3), (-2, 0, 1, 3), (.04, .2, .7)):
        exponents = tuple((k * k + 2 * k) % 5 - 2 for k in range(prime))
        data = active_profiles(prime, n, ramification, charge, exponents)
        direct = direct_full_pressure(data, prime, ramification, scale)
        value = filter_pressure(prime, n, ramification, charge, scale,
                                data['beta_orders'], data['infinity_active'])
        error = abs(value - direct)
        assert error < 1e-10
        worst_direct = max(worst_direct, error)
        clearing = 10**6
        shifted = filter_pressure(
            prime, n, ramification, charge, scale,
            [order + 2 * clearing for order in data['beta_orders']],
            [order - 2 * clearing for order in data['infinity_active']])
        error = abs(shifted - value)
        assert error < 1e-8
        worst_clearing = max(worst_clearing, error)
        cases += 1
    result = {
        'status': 'Floating arithmetic controls, not a rigorous pressure certificate.',
        'complete_profile_partition_comparisons': cases,
        'largest_absolute_profile_comparison_error': worst_direct,
        'common_scalar_clearing_comparisons': cases,
        'clearing_power': 10**6,
        'largest_absolute_clearing_error': worst_clearing,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
