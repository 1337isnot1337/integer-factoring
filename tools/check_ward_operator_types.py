#!/usr/bin/env python3
"""Literal tensor controls for the combined Ward type-block ranks."""

import argparse
import hashlib
from itertools import product
import json
from pathlib import Path

from check_modular_crossing import rank_two
from check_noncommutative_crossing_ward_rank_descent import entry, labels
from pfaffian_particle_hole import local_exterior_table, particle_hole_mask
from ward_operator_types import local_operators, operator_score, parity_filter, rectangular_operator_profile


def local_block_controls():
    checked = 0
    for eta_axis in ('x', 'y'):
        pairs, stride = 3, 11
        eta = 1 if eta_axis == 'x' else stride
        table = local_exterior_table(3, 1 | (1 << stride), 1)
        ordinary, operated = [[0] * 16 for _ in range(16)], [[0] * 16 for _ in range(16)]
        for a, y, u, v in product(range(4), repeat=4):
            i, j = particle_hole_mask(a, y), particle_hole_mask(v, u)
            value = table[((a * 4 + y) * 4 + u) * 4 + v]
            acted = parity_filter(value, a.bit_count(), stride)
            if a == 3:
                acted ^= table[((0 * 4 + y) * 4 + u) * 4 + v] << eta
            ordinary[i][j], operated[i][j] = value, acted
        g, h = local_operators(pairs, eta_axis)
        sectors = ((0,), (15,), (1, 2, 4, 8), (7, 11, 13, 14), (3, 5, 9, 6, 10, 12))
        for block, indices in enumerate(sectors):
            for row, i in enumerate(indices):
                for col, j in enumerate(indices):
                    if block < 2:
                        assert ordinary[i][j] == 2 and operated[i][j] == 0
                    else:
                        matrix = g if block < 4 else h
                        assert matrix[row][col] == (ordinary[i][j], operated[i][j])
                    checked += 1
        for i, j in product(range(16), repeat=2):
            if not any(i in sector and j in sector for sector in sectors):
                assert ordinary[i][j] == operated[i][j] == 0
                checked += 1
    return checked


def literal_profile(pairs, length):
    first, second = labels(pairs, length), labels(pairs, length - 2)
    size, preimage = len(first), len(second)
    stride = 2 * pairs + 1
    rows, columns = {}, {}
    entries = 0
    for a, b, c, d in product(range(size), range(preimage), range(size), range(size)):
        x, y, u, v = first[a], second[b], first[c], first[d]
        ordinary = entry(pairs, x, y, u, v)
        value = parity_filter(ordinary, length, stride)
        for site, occupied in enumerate(x):
            if occupied == 3:
                source = x[:site] + (0,) + x[site + 1:]
                value ^= entry(pairs, source, y, u, v) << stride
        value <<= stride
        while value:
            bit = value & -value
            exponent = bit.bit_length() - 1
            key = exponent % stride, exponent // stride
            rows.setdefault(key, [0] * (size * preimage))[a * preimage + b] ^= 1 << (d * size + c)
            columns.setdefault(key, [0] * (size * size))[a * size + c] ^= 1 << (d * preimage + b)
            value ^= bit
        entries += 1
    row_ranks = {key: rank_two(matrix) for key, matrix in rows.items() if rank_two(matrix)}
    col_ranks = {key: rank_two(matrix) for key, matrix in columns.items() if rank_two(matrix)}
    typed_row = rectangular_operator_profile(pairs, length, length - 2, length, length, 'y')
    typed_col = rectangular_operator_profile(pairs, length, length, length, length - 2, 'x')
    assert row_ranks == typed_row['profile']
    assert col_ranks == {(b, a): rank for (a, b), rank in typed_col['profile'].items()}
    return {'pairs': pairs, 'length': length, 'literal_polynomial_entries': entries,
            'row_ranks': [[a, b, rank] for (a, b), rank in sorted(row_ranks.items())],
            'column_ranks': [[a, b, rank] for (a, b), rank in sorted(col_ranks.items())]}


def code_controls():
    existing = json.loads(Path('evidence/ward-charge-filtration.json').read_text())
    checked = []
    for case in existing['literal_code_controls']:
        if case.get('forms', 1) != 1:
            continue
        score = operator_score(case['pairs'], case['length'], case['labels'])
        bounds = {(a, b): (r, c) for a, b, r, c in score['profile']}
        assert all(r <= bounds.get((a, b), (0, 0))[0] and c <= bounds.get((a, b), (0, 0))[1]
                   for a, b, r, c in case['exact_jets'])
        checked.append({'pairs': case['pairs'], 'length': case['length'], 'labels': case['labels'],
                        'combined_cost_bound': score['integer_rank_sum_upper_bound'],
                        'differential_cost_bound': score['previous_differential_integer_rank_sum_upper_bound']})
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in (Path(__file__).name, 'ward_operator_types.py')},
              'local_operator_block_entries': local_block_controls(),
              'literal_full_tensor_controls': [literal_profile(2, 2), literal_profile(3, 3)],
              'restricted_code_controls': code_controls()}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
