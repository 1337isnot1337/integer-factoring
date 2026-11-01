#!/usr/bin/env python3
"""Independent determinant identities and literal controls for Ward charge bounds."""

import hashlib
import json
from itertools import combinations, permutations, product
from math import comb
from pathlib import Path
import random

from check_modular_crossing import rank_two
from check_noncommutative_crossing_ward_rank_descent import actual_code_profile
from explore_symplectic_boundary_codes import boundary_code, half_orientation_code
from ward_channel_free_bounds import channel_free_score
from ward_differential_bounds import differential_score


def determinant(rows):
    return int(rank_two(rows) == len(rows))


def principal(rows, mask):
    indices = [i for i in range(len(rows)) if mask >> i & 1]
    return [sum(((rows[i] >> j) & 1) << k for k, j in enumerate(indices)) for i in indices]


def symmetric(size, code):
    rows = [0] * size
    bit = 0
    for i in range(size):
        for j in range(i, size):
            if code >> bit & 1:
                rows[i] |= 1 << j
                rows[j] |= 1 << i
            bit += 1
    return rows


def poly_product(first, second):
    result = set()
    for a in first:
        for b in second:
            term = tuple(sorted(a + b))
            if term in result:
                result.remove(term)
            else:
                result.add(term)
    return result


def formal_determinant(matrix, indices):
    result = set()
    for perm in permutations(indices):
        term = {()}
        for row, col in zip(indices, perm):
            term = poly_product(term, matrix[row][col])
        result.symmetric_difference_update(term)
    return result


def formal_translation_checks():
    records = []
    for size in range(1, 7):
        variables = {(i, j): n for n, (i, j) in enumerate(
            (pair for i in range(size) for pair in ((i, j) for j in range(i, size))))}
        count = len(variables)
        first = [[{(variables[min(i, j), max(i, j)],)} for j in range(size)] for i in range(size)]
        second = [[{(count + variables[min(i, j), max(i, j)],)} for j in range(size)] for i in range(size)]
        combined = [[first[i][j] | second[i][j] for j in range(size)] for i in range(size)]
        expected = set()
        for mask in range(1 << size):
            left = [i for i in range(size) if mask >> i & 1]
            right = [i for i in range(size) if not mask >> i & 1]
            expected.symmetric_difference_update(poly_product(
                formal_determinant(first, left), formal_determinant(second, right)))
        actual = formal_determinant(combined, list(range(size)))
        assert actual == expected
        records.append({'size': size, 'independent_scalar_variables': 2 * count,
                        'surviving_formal_monomials': len(actual)})
    return records


def border_checks():
    rng = random.Random(492781)
    checked = 0
    for size in range(1, 6):
        q = size * (size + 1) // 2
        for _ in range(160):
            base = symmetric(size, rng.randrange(1 << q))
            columns = [rng.randrange(1 << size) for _ in range(rng.randrange(size + 2))]
            difference = 0
            for selected in range(1 << len(columns)):
                rows = base.copy()
                for j, column in enumerate(columns):
                    if selected >> j & 1:
                        for i in range(size):
                            if column >> i & 1:
                                rows[i] ^= column
                difference ^= determinant(rows)
            bordered = [row | sum(((column >> i) & 1) << (size + j)
                                  for j, column in enumerate(columns)) for i, row in enumerate(base)]
            bordered += columns
            assert difference == determinant(bordered)
            checked += 1
    return checked


def augmentation_checks():
    rng = random.Random(68053)
    checked = 0
    for size in range(1, 4):
        q = size * (size + 1) // 2
        total = 1 << q
        matrices = [symmetric(size, code) for code in range(total)]
        full = (1 << size) - 1
        for _ in range(5):
            offset, operators = rng.randrange(total), rng.randrange(1 << size)
            for degree in range(size + 1):
                values = []
                for code in range(total):
                    value = 0
                    for subset in range(1 << size):
                        if subset.bit_count() == degree and subset & operators == subset:
                            value ^= determinant(principal(matrices[code ^ offset], full ^ subset))
                    values.append(value)
                coefficients = values.copy()
                for bit in range(q):
                    for mask in range(total):
                        if mask >> bit & 1:
                            coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
                assert not any(value and mask.bit_count() > size - degree
                               for mask, value in enumerate(coefficients))
                for left_degree in range(size + 2):
                    for right_degree in range(size + 2):
                        left = [x for x in range(total) if x.bit_count() >= left_degree]
                        right = [x for x in range(total) if x.bit_count() >= right_degree]
                        rows = [sum((int(not x & y) * coefficients[x | y]) << j
                                    for j, y in enumerate(right)) for x in left]
                        bound = sum(comb(size, h) for h in range(left_degree + degree, size - right_degree + 1))
                        assert rank_two(rows) <= bound
                        checked += 1
    return checked


def literal_controls():
    records = []
    for pairs, length in ((3, 3), (4, 4)):
        full = boundary_code(pairs, length, True)
        for half in (False, True):
            code = half_orientation_code(full) if half else full
            exact = actual_code_profile(code)
            candidates = [differential_score(pairs, length, code['labels']),
                          channel_free_score(pairs, length, 1, code['labels'], pairs % 2)]
            for bound in candidates:
                lookup = {(a, b): (r, c) for a, b, r, c in bound['profile']}
                assert all(r <= lookup.get(key, (0, 0))[0] and c <= lookup.get(key, (0, 0))[1]
                           for key, (r, c) in exact.items())
            records.append({'pairs': pairs, 'length': length, 'labels': code['labels'],
                            'half': half, 'exact_jets': [[a, b, r, c] for (a, b), (r, c) in sorted(exact.items())],
                            'differential_cost_upper': candidates[0]['integer_rank_sum_upper_bound'],
                            'channel_free_cost_upper': candidates[1]['integer_rank_sum_upper_bound']})
    overlaps = json.loads(Path('evidence/overlapping-symplectic-code-profiles.json').read_text())
    for case in overlaps['cases']:
        k = case['forms']
        gram = [0] * k
        for column in case['form_column_types']:
            for i in range(k):
                if column >> i & 1:
                    gram[i] ^= column
        bound = channel_free_score(case['pairs'], case['length'], k, case['labels'], rank_two(gram))
        lookup = {(a, b): (r, c) for a, b, r, c in bound['profile']}
        assert all(r <= lookup.get((a, b), (0, 0))[0] and c <= lookup.get((a, b), (0, 0))[1]
                   for a, b, r, c in case['literal_centered_jet_ranks'])
        records.append({'pairs': case['pairs'], 'length': case['length'], 'labels': case['labels'],
                        'forms': k, 'joint_order': bound['joint_order_used'],
                        'exact_jets': case['literal_centered_jet_ranks'],
                        'channel_free_cost_upper': bound['integer_rank_sum_upper_bound']})
    return records


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'formal_symmetric_determinant_translation': formal_translation_checks(),
              'rank_one_bordered_determinant_checks': border_checks(),
              'operator_coefficient_augmentation_rank_checks': augmentation_checks(),
              'literal_code_controls': literal_controls(),
              'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in (Path(__file__).name, 'ward_channel_free_bounds.py',
                                            'ward_differential_bounds.py', 'ward_rectangular_bounds.py')}}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
