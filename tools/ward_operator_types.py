#!/usr/bin/env python3
"""Exact global Ward operator ranks on particle-hole type blocks."""

import argparse
from functools import lru_cache
import hashlib
import json
from math import comb, isqrt, log2
from pathlib import Path

from check_modular_crossing import rank_two
from complement_shell_tensor import multiply
from pfaffian_grade_types import G_GRADES, H_GRADES, H_STATES, type_basis
from pfaffian_particle_hole import local_exterior_table, one_particle_matrix
from screen_ward_code_bounds import radical_dimensions
from ward_differential_bounds import differential_score


def parity_filter(value, external_grade, stride):
    result = 0
    while value:
        bit = value & -value
        exponent = bit.bit_length() - 1
        if (external_grade + exponent % stride + exponent // stride) & 1:
            result ^= bit
        value ^= bit
    return result


@lru_cache(None)
def local_operators(pairs, eta_axis):
    stride = 2 * pairs + 5
    g = one_particle_matrix(3, 1 | (1 << stride), 1)
    table = local_exterior_table(3, 1 | (1 << stride), 1)
    h = tuple(tuple(table[((x * 4 + y) * 4 + u) * 4 + v]
                    for v, u in H_STATES) for x, y in H_STATES)
    eta_shift = 1 if eta_axis == 'x' else stride
    result = []
    for matrix, grades, creation in ((g, G_GRADES, {2: 1}),
                                     (h, H_GRADES, {1: 0, 5: 4})):
        operator = [[parity_filter(value, grades[row][0], stride) for value in values]
                    for row, values in enumerate(matrix)]
        for target, source in creation.items():
            for col, value in enumerate(matrix[source]):
                operator[target][col] ^= value << eta_shift
        result.append(tuple(tuple((matrix[i][j], operator[i][j]) for j in range(len(matrix)))
                            for i in range(len(matrix))))
    return tuple(result)


def rectangular_operator_profile(pairs, ra, ry, cv, cu, eta_axis):
    stride = 2 * pairs + 5
    eta_shift = 1 if eta_axis == 'x' else stride
    g, h = local_operators(pairs, eta_axis)
    total = {}
    largest = entries = 0
    records = []
    for scalars in range(pairs + 1):
        for odd in range(pairs - scalars + 1):
            middle = pairs - scalars - odd
            rows = type_basis(odd, middle, ra - scalars, ry - scalars)
            cols = type_basis(odd, middle, cv - scalars, cu - scalars)
            if not rows or not cols:
                continue
            multiplicity = comb(pairs, scalars) * comb(pairs - scalars, odd) * 2 ** (scalars + odd)
            matrices, cache = (g,) * odd + (h,) * middle, {}
            jets = {}
            for i, row in enumerate(rows):
                for j, col in enumerate(cols):
                    pattern = tuple(sorted(matrix[a][b] for matrix, a, b in zip(matrices, row, col)))
                    value = cache.get(pattern)
                    if value is None:
                        product_value, inserted = 1, 0
                        for ordinary, acted in pattern:
                            inserted = multiply(inserted, ordinary) ^ multiply(product_value, acted)
                            product_value = multiply(product_value, ordinary)
                        value = inserted << (scalars + eta_shift)
                        cache[pattern] = value
                    while value:
                        bit = value & -value
                        exponent = bit.bit_length() - 1
                        key = exponent % stride, exponent // stride
                        jets.setdefault(key, [0] * len(rows))[i] ^= 1 << j
                        value ^= bit
            ranks = {key: rank_two(matrix) for key, matrix in jets.items()}
            for key, rank in ranks.items():
                total[key] = total.get(key, 0) + multiplicity * rank
            largest = max(largest, len(rows), len(cols))
            entries += len(rows) * len(cols)
            records.append({'scalars': scalars, 'odd': odd, 'middle': middle,
                            'multiplicity': multiplicity, 'row_dimension': len(rows),
                            'column_dimension': len(cols),
                            'coefficient_ranks': [[a, b, r] for (a, b), r in sorted(ranks.items()) if r]})
    return {'profile': {key: rank for key, rank in total.items() if rank},
            'largest_type_dimension': largest, 'evaluated_type_entries': entries,
            'types': records}


def operator_score(pairs, length, labels):
    row = rectangular_operator_profile(pairs, length, length - 2, length, length, 'y')
    col = rectangular_operator_profile(pairs, length, length, length, length - 2, 'x')
    previous = differential_score(pairs, length, labels)
    old = {(a, b): (r, c) for a, b, r, c in previous['profile']}
    total, profile = 0, []
    for (a, b), (r, c) in sorted(old.items()):
        r = min(r, row['profile'].get((a, b), 0))
        c = min(c, col['profile'].get((b, a), 0))
        if r or c:
            product_ = r * c
            root = isqrt(product_)
            total += root + int(root * root < product_)
            profile.append([a, b, r, c])
    return {'pairs': pairs, 'length': length, 'labels': labels,
            'integer_rank_sum_upper_bound': total,
            'previous_differential_integer_rank_sum_upper_bound': previous['integer_rank_sum_upper_bound'],
            'rho_of_upper_bound': log2(total) / (2 * log2(labels)),
            'largest_type_dimension': max(row['largest_type_dimension'], col['largest_type_dimension']),
            'evaluated_type_entries': row['evaluated_type_entries'] + col['evaluated_type_entries'],
            'profile': profile,
            'raw_row_operator_ranks': [[a, b, r] for (a, b), r in sorted(row['profile'].items())],
            'raw_column_operator_ranks': [[b, a, r] for (a, b), r in sorted(col['profile'].items())]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append', default=[])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with args.output.open('w') as output:
        for case in args.case or ('3:3', '4:4', '6:4', '8:4'):
            pairs, length = map(int, case.split(':'))
            result = operator_score(pairs, length, radical_dimensions(pairs)[length])
            result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            output.write(json.dumps(result, sort_keys=True) + '\n')
            output.flush()
            print(json.dumps({key: value for key, value in result.items() if key not in ('profile', 'raw_row_operator_ranks', 'raw_column_operator_ranks')}), flush=True)


if __name__ == '__main__':
    main()
