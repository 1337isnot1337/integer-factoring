#!/usr/bin/env python3
"""Exact rational dual certificate for the sparse Gaussian bulk model."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


def unit(index, size=21):
    return [Fraction(int(j == index)) for j in range(size)]


def add(*vectors):
    return [sum(row) for row in zip(*vectors)]


def scale(vector, scalar):
    return [scalar * value for value in vector]


def check():
    # B,X,V,Y; alpha 0:4, beta 4:8, gamma 8:12, delta 12:16.
    names = [prefix + '_' + boundary for prefix in ('alpha', 'beta', 'gamma', 'delta')
             for boundary in ('B', 'X', 'V', 'Y')] + ['ell', 'k1', 'k2', 'k3', 'k4']
    variables = {name: unit(i) for i, name in enumerate(names)}
    def vector(**coefficients):
        return add(*[scale(variables[name], Fraction(value)) for name, value in coefficients.items()])
    d1 = vector(ell=2, alpha_X=1, alpha_Y=1, beta_V=1, beta_B=1)
    d2 = vector(ell=1, alpha_V=1, alpha_B=1, beta_X=1, beta_Y=1)
    d3 = vector(gamma_X=1, gamma_B=1, delta_V=1, delta_Y=1)
    d4 = vector(ell=-1, gamma_V=1, gamma_Y=1, delta_X=1, delta_B=1)
    ks = [variables['k' + str(j)] for j in range(1, 5)]
    orders = [value for k, d in zip(ks, (d1, d2, d3, d4))
              for value in (k, add(d, scale(k, -1)))]
    occupancies = tuple(map(Fraction, (1, 1, 1, 1, 1, 0))) + (Fraction(1, 2), Fraction(1, 2))
    weighted = add(*[scale(order, weight) for order, weight in zip(orders, occupancies)])
    assert weighted == add(d1, d2, variables['k3'], scale(d4, Fraction(1, 2)))
    gauge = add(*[variables[name] for name in names[:16]])
    lhs = add(gauge, scale(weighted, -1))
    inequalities = [vector(gamma_X=1, delta_V=1, ell=-1, k3=-1),
                    vector(gamma_B=1, delta_Y=1, k3=-1),
                    vector(gamma_B=1, delta_X=1, ell=-2),
                    vector(gamma_X=1, delta_B=1, ell=-1),
                    vector(gamma_V=1, delta_Y=1),
                    vector(gamma_Y=1, delta_V=1, ell=-1)]
    rhs = scale(add(*inequalities), Fraction(1, 2))
    assert lhs == rhs
    assert all(not value for value in lhs[:8])
    assert Fraction(3, 4) ** 2 > Fraction(1, 2)
    nonzero = {name: str(value) for name, value in zip(names, lhs) if value}
    return {'status': 'pass', 'model': 'constant-species sparse Gaussian bulk only',
            'entropy_occupancies_divided_by_p': list(map(str, occupancies)),
            'dual_inequality_weights': ['1/2'] * 6,
            'dual_identity_coefficients': nonzero,
            'entropy_lower_bound': '5 H(p) + 2 H(p/2) >= 6 H(p)',
            'rho_lower_bound': '3/4',
            'exceeds_half_threshold_squared': str(Fraction(3, 4) ** 2 - Fraction(1, 2)),
            'scope': 'No claim for unbounded index slopes, arbitrary physical-letter gauges, or energy projection.',
            'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = json.dumps(check(), indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report)
    print(report, end='')


if __name__ == '__main__':
    main()
