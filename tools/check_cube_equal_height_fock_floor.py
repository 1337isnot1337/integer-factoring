#!/usr/bin/env python3
"""Exact shared generic minors for the equally weighted multi-mode floor."""
from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import random
import time

from check_cube_gaussian_profile_budget import floor_certificate
from cube_homogeneous_gaussian import ZERO, ONE, monomial, multiply, valuation


def determinant(matrix, rows, columns):
    rows, columns = tuple(rows), tuple(columns)

    @lru_cache(None)
    def visit(position, selected):
        if position == len(rows):
            return ONE
        value = ZERO
        for j, column in enumerate(columns):
            if selected >> j & 1:
                value ^= multiply(matrix[rows[position]][column],
                                  visit(position + 1, selected ^ (1 << j)))
        return value

    return visit(0, (1 << len(columns)) - 1)


def sampled_kernel(modes, degree, rng):
    terms = [(a, b, degree - a - b)
             for a in range(degree + 1) for b in range(degree - a + 1)]
    order = 8 * modes
    result = [[ZERO] * order for _ in range(order)]
    for i in range(order):
        for j in range(i + 1, order):
            value = frozenset(term for term in terms if rng.randrange(2))
            expected = int((i // modes) ^ (j // modes) == 7
                           and i % modes == j % modes)
            if len(value) % 2 != expected:
                value ^= monomial(degree, 0, 0)
            result[i][j] = result[j][i] = value
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='evidence/cube-equal-height-fock-floor.json')
    args = parser.parse_args()
    started = time.monotonic()
    rng = random.Random(2026110128)
    records = []
    entry_controls = minor_controls = orientation_controls = 0
    for modes, degree, samples in ((1, 1, 3), (1, 3, 3), (2, 1, 3),
                                    (2, 2, 3), (3, 1, 3), (4, 1, 2)):
        for sample in range(samples):
            source = sampled_kernel(modes, degree, rng)
            for i, row in enumerate(source):
                for j, value in enumerate(row):
                    assert source[j][i] == value
                    assert all(sum(term) == degree for term in value)
                    expected = int((i // modes) ^ (j // modes) == 7
                                   and i % modes == j % modes)
                    assert len(value) % 2 == expected
                    if i == j:
                        assert not value
                    entry_controls += 1
            rows = list(range(modes))
            columns = [7 * modes + i for i in range(modes)]
            common = determinant(source, rows, columns)
            assert common and len(common) % 2 == 1
            assert all(sum(term) == modes * degree for term in common)
            orders = [valuation(common, axis) for axis in range(3)]
            assert sum(orders) <= modes * degree
            minor_controls += 1
            pairs = []
            for omitted in range(3):
                axes = [axis for axis in range(3) if axis != omitted]
                vertex = 1 << omitted
                rows = [v * modes + i for v in (0, vertex) for i in range(modes)]
                columns = [v * modes + i for v in (7, vertex ^ 7) for i in range(modes)]
                for axis in axes:
                    assert all(not ((i // modes) >> axis & 1) for i in rows)
                    assert all((i // modes) >> axis & 1 for i in columns)
                    orientation_controls += 1
                value = determinant(source, rows, columns)
                assert value and len(value) % 2 == 1
                assert all(sum(term) == 2 * modes * degree for term in value)
                pair_orders = [valuation(value, axis) for axis in axes]
                assert sum(pair_orders) <= 2 * modes * degree
                pairs.append({'axes': axes, 'orders': pair_orders,
                              'determinant_monomials': len(value)})
                minor_controls += 1
            records.append({'modes': modes, 'degree': degree, 'sample': sample,
                            'common_minor_orders': orders, 'pair_minors': pairs})
    certificate = floor_certificate()
    assert len(certificate['positive_integer_coefficients']) == 355
    paths = [Path(__file__), Path(__file__).with_name('cube_homogeneous_gaussian.py'),
             Path(__file__).with_name('check_cube_gaussian_profile_budget.py')]
    report = {'schema': 'cube-equal-height-fock-floor-v1',
              'scope': 'All equally weighted full Fock alphabets; common unweighted profiles only.',
              'generic_shared_minors': minor_controls,
              'entry_endpoint_and_degree_controls': entry_controls,
              'shared_cut_orientations': orientation_controls,
              'kernels': records,
              'all_modes_profile_floor': '3/5',
              'positive_coefficients': len(certificate['positive_integer_coefficients']),
              'minimum_positive_coefficient': min(certificate['positive_integer_coefficients']),
              'new_factoring_score': False,
              'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in paths},
              'elapsed_seconds': time.monotonic() - started}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'kernels': len(records),
                      'generic_shared_minors': minor_controls,
                      'entry_controls': entry_controls,
                      'positive_coefficients': report['positive_coefficients'],
                      'seconds': report['elapsed_seconds']}))


if __name__ == '__main__':
    main()
