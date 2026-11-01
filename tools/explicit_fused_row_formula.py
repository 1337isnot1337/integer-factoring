#!/usr/bin/env python3
"""Closed bit x symmetric-power fusion rows in the repository's convention."""
from __future__ import annotations

from f16_identity_operator_profiles import MUL, INV, multiply
from fused_six_vertex_crossings import power, pscale


def explicit_row(m, theta, shift=0):
    """Return the degree-one row and its exact removed monic common factor.

    The formula is proved over F2[theta^±1,z]. Theta != 1 ensures the
    displayed degree-m-1 scalar is the complete row gcd.
    """
    if m < 1 or not 1 < theta < 16:
        raise ValueError('Use positive m and nonzero theta different from one.')
    factor = power(theta, 2 * shift)
    common = {0: 1}
    for j in range(1, m):
        common = multiply(common, {0: INV[theta], 1: power(theta, 1 + 2 * shift - 2 * j)})
    leading = common[max(common)]
    monic = pscale(common, INV[leading])
    matrix = {}
    for b in range(m + 1):
        entries = {
            (b, 0, 0, b): {0: power(theta, b - 1),
                            1: MUL[factor][power(theta, 1 - b)]},
            (b, 1, 1, b): {0: power(theta, m - b - 1),
                            1: MUL[factor][power(theta, 1 - m + b)]}}
        if b:
            entries[b - 1, 1, 0, b] = {1: MUL[factor][theta ^ power(theta, 2 * b - 2 * m - 1)]}
        if b < m:
            entries[b + 1, 0, 1, b] = {0: power(theta, m) ^ power(theta, m - 2 * b - 2)}
        for key, value in entries.items():
            matrix[key] = pscale({k: v for k, v in value.items() if v}, leading)
    return {'matrix': matrix, 'common_factor': monic, 'degree': 1}
