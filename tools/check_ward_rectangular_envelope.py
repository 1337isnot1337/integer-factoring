#!/usr/bin/env python3
"""Protected blocks for the computed Ward rectangular-bound envelope."""

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log2
from pathlib import Path

from stable_pfaffian_bounds import stable_multiplicities
from ward_rectangular_bounds import ANCILLAS, rectangular_block


def entropy(value):
    return -value * log2(value) - (1 - value) * log2(1 - value) if 0 < value < 1 else 0.0


def fusion_entropy(u, v):
    t = (1 - v) / (2 - u - v)
    value = entropy(t) + .5 * (1 - u - v + t * entropy(u / t) + (1 - t) * entropy(v / (1 - t)))
    defect = 1 - u - v
    h3 = -sum(x * log2(x) for x in (u, v, defect) if x > 0)
    alternate = .5 * (h3 + (1 + defect) * entropy((u + defect) / (1 + defect)) + defect)
    assert abs(value - alternate) < 1e-12
    return value


def protected_blocks():
    blocks = caps = 0
    cases = []
    for pairs in range(3, 13):
        multiplicities = {(m, n, twist): count for m, n, twist, count in stable_multiplicities(pairs)}
        for forms in range(1, pairs // 2 + 1):
            for length in range(2 * forms + 1, pairs + 1):
                excess = length - 2 * forms
                a = excess // (4 if 3 * excess <= 2 * pairs else 3)
                b = a + 2 * forms
                m, n = 2 * excess + 1 - 4 * a, 2 * pairs + 1 - 2 * excess + 2 * a
                specifications = (
                    (m, n, a, (excess, excess, length, length), a),
                    (m, n, a, (length, length, length, length), b),
                    (m, n - 4 * forms, b, (length, length, length, length), a),
                    (m + 2 * forms, n - 2 * forms, a, (excess, length, length, excess), a),
                )
                for first, second, twist, grades, y_degree in specifications:
                    count = multiplicities[first, second, twist]
                    assert count > 0
                    for ancillary_row, ancillary_column in ANCILLAS:
                        profile, _ = rectangular_block(first, second, *grades, pairs, ancillary_row, ancillary_column)
                        assert {(x, y): rank for x, y, rank in profile}.get((0, y_degree)) == 1
                        blocks += 1
                    ra, ry, cv, cu = grades
                    cap = min(comb(2 * pairs, ra) * comb(2 * pairs, ry),
                              comb(2 * pairs, cv) * comb(2 * pairs, cu))
                    assert count <= 16 * (2 * pairs + 1) ** 4 * cap
                    caps += 1
                cases.append((pairs, forms, length))
    return {'parameter_cases': len(cases), 'all_ancilla_rank_one_blocks': blocks,
            'multiplicity_vs_physical_dimension_controls': caps}


def entropy_controls():
    checked = 0
    minimum = 2.0
    for index in range(1, 201):
        excess = index / 200
        x = excess / (4 if excess <= 2 / 3 else 3)
        u, v = excess - 2 * x, 1 - excess + x
        base = fusion_entropy(u, v)
        for step in range(41):
            kappa = (1 - excess) * step / 80
            assert kappa <= v - u + 1e-12
            first = fusion_entropy(u + kappa, v - kappa)
            second = fusion_entropy(u, v - 2 * kappa)
            assert min(first, second) + 1e-12 >= base
            checked += 1
        alpha = excess / 2
        t = 7 * alpha / 5 if excess <= 2 / 3 else alpha
        lower = entropy(t) + .5 * (1 - u - v + t * entropy(u / t) + (1 - t) * entropy(v / (1 - t)))
        assert lower + 1e-12 >= 29 * entropy(alpha) / 20
        minimum = min(minimum, lower / (2 * entropy(alpha)))
    assert 2 ** 114 * 3 ** 9 > 5 ** 50
    assert Fraction(7, 8) ** 4 > Fraction(1, 2)
    assert Fraction(3, 2) ** 3 < 4
    return {'floating_entropy_controls': checked,
            'sample_minimum_protected_rho': minimum,
            'exact_integer_endpoint_and_log_estimates': 3,
            'scope': 'The floating controls illustrate the analytic inequalities; the note supplies their proof.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'protected_blocks': protected_blocks(), 'entropy_controls': entropy_controls()}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
