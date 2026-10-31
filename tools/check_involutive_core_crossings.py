#!/usr/bin/env python3
"""Exact involutive-core, polynomiality, and typed-endpoint controls."""
from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_residual_quartic_codes import apply, residual_operator, restrict
from check_schur_downset_crossing import sparse_kernel
from involutive_core_crossing import CORE_MAPS, build
from weighted_smith_oracle import polynomial_exact_divide, polynomial_multiply


ROOT = Path(__file__).resolve().parent.parent


def core_gate():
    maps = [[sum(1 << state for state in states) for states in leg] for leg in CORE_MAPS]
    gate = restrict(residual_operator(), maps[2], maps[0], maps[1], maps[3])
    expected = [1 << i for i in range(9)]
    expected[1] ^= 1 << 3
    expected[7] ^= 1 << 5
    assert gate == expected
    assert all(apply(gate, apply(gate, 1 << j)) == 1 << j for j in range(9))
    return gate


def independent_exterior_constructor(indices, ambient):
    n = len(indices)
    kernel = sparse_kernel(indices, ambient)
    matrix = [row[:2 * n] for row in kernel[2 * n:]]
    gaussian = {0: 1}
    for i, j in product(range(2 * n), repeat=2):
        coefficient = matrix[i][j]
        if not coefficient:
            continue
        row_bit = 3 * n + i if i < n else i
        column_bit = 2 * n + j if j < n else j - n
        edge = (1 << row_bit) | (1 << column_bit)
        additions = [(mask | edge, polynomial_multiply(value, coefficient))
                     for mask, value in gaussian.items() if not mask & edge]
        for mask, value in additions:
            gaussian[mask] = gaussian.get(mask, 0) ^ value
            if not gaussian[mask]:
                del gaussian[mask]
    quartic = 1 | (1 << (2 * n - 1)) | (1 << (2 * n)) | (1 << (4 * n - 1))
    numerator = {mask: polynomial_multiply(5, value) for mask, value in gaussian.items()}
    for mask, value in gaussian.items():
        if not mask & quartic:
            numerator[mask | quartic] = numerator.get(mask | quartic, 0) ^ value
    lookups = [{state: core for core, states in enumerate(leg) for state in states} for leg in CORE_MAPS]
    bulk = 1 << (n - 2)
    q = 3 * bulk
    values = [0] * q ** 4
    for mask, value in numerator.items():
        labels = []
        input_low = False
        for leg in range(4):
            physical = (mask >> (n * leg)) & ((1 << n) - 1)
            state = (physical & 1) | ((physical >> (n - 1)) << 1)
            if state not in lookups[leg]:
                break
            labels.append(((physical >> 1) & (bulk - 1)) + bulk * lookups[leg][state])
            input_low |= leg == 2 and state == 1
        else:
            if input_low:
                value ^= value << 1
            x, y, u, v = labels
            values[x + q * y + q * q * u + q ** 3 * v] ^= value
    return [polynomial_exact_divide(value, 5) for value in values]


def two_copy_type_control(data):
    q = data['q']
    endpoint = [value.bit_count() & 1 for value in data['values']]
    types = defaultdict(list)
    for first, second in product(range(q), repeat=2):
        types[tuple(sorted((first, second)))].append((first, second))

    def index(x, y, u, v):
        return x + q * y + q * q * u + q ** 3 * v

    checks = 0
    for words in types.values():
        for x, y, u, v in product(words, repeat=4):
            actual = endpoint[index(x[0], y[0], u[0], v[0])] & endpoint[index(x[1], y[1], u[1], v[1])]
            assert actual == int(x == v and y == u)
            checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--seed-output', type=Path)
    args = parser.parse_args()
    gate = core_gate()
    cases = [((0, 1), 2, 0, 0, 0, None, 'direct'),
             ((0, 1, 2), 4, 0, 0, 0, None, 'direct'),
             ((0, 1, 2, 3), 4, 0, 0, 0, None, 'direct'),
             ((0, 1, 3, 7), 8, 0, 0, 0, None, 'direct'),
             ((0, 1, 2), 4, 0, 0, 0, -1, 'direct'),
             ((0, 1, 2), 4, 2, -3, 1, None, 'direct'),
             ((0, 1, 2), 4, -2, 1, -4, 2, 'direct'),
             ((0, 1), 2, 1, 0, 0, None, 'direct'),
             ((0, 1, 2), 4, 0, 0, 0, None, 'ramified-repair'),
             ((0, 1, 2), 4, 0, 0, 0, None, 'temporal')]
    records = []
    seed = None
    for case in cases:
        data = build(*case)
        if case[0] in ((0, 1, 2), (0, 1, 2, 3)) and case[2:] == (0, 0, 0, None, 'direct'):
            independent = independent_exterior_constructor(case[0], case[1])
            assert independent == data['values']
            data['independent_exterior_coefficients_checked'] = len(independent)
        data['two_copy_full_type_endpoint_checks'] = two_copy_type_control(data)
        if case == ((0, 1, 2), 4, 0, 0, 0, None, 'direct'):
            seed = dict(data)
        values = data.pop('values')
        data['coefficient_sha256'] = sha256(json.dumps(values, separators=(',', ':')).encode()).hexdigest()
        records.append(data)
    report = {'status': 'exact polynomial crossing constructions; no pressure or exponent claim',
              'field': 2, 'core_transfer_binary_rows': gate, 'core_transfer_square_is_identity': True,
              'cases': records,
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest()
                                for name in ('tools/involutive_core_crossing.py',
                                             'tools/check_involutive_core_crossings.py',
                                             'tools/check_residual_quartic_codes.py',
                                             'tools/check_schur_downset_crossing.py',
                                             'tools/weighted_smith_oracle.py')}}
    rendered = json.dumps(report, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(rendered)
    if args.seed_output:
        args.seed_output.write_text(json.dumps(seed, separators=(',', ':')) + '\n')
    print(rendered)


if __name__ == '__main__':
    main()
