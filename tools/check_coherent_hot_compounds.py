#!/usr/bin/env python3
"""Independent small-map checks for coherent occupation-block crossings."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path

from check_block_smith_crossings import toeplitz_rank
from check_cancellation_module_automata import pmul
from coherent_hot_compounds import height, one_particle, profiles, support, validate_maps
from coherent_inverse_systems import derivative, inverse_system
from eight_letter_reflection_retraction import displayed_maps, homogeneous_restriction


def reverse_map(vector):
    return sum(1 << int(f'{word:04b}'[::-1], 2) for word in support(vector))


def determinant(matrix, row_mask, col_mask):
    rows, cols = support(row_mask), support(col_mask)
    assert len(rows) == len(cols)
    result = 0
    for order in permutations(cols):
        value = 1
        for i, j in zip(rows, order):
            value = pmul(value, matrix[i][j])
        result ^= value
    return result


def compressed_entry(matrix, modes, p, l, cut, row_pair, col_pair):
    x, y = row_pair
    v, u = col_pair
    row_vectors = [a | (b << modes) for a, b in product(support(l[x]),
                        support(l[y]) if cut == 'zero' else support(p[y]))]
    col_vectors = [a | (b << modes) for a, b in product(
                        support(p[u]) if cut == 'zero' else support(l[u]), support(p[v]))]
    value = 0
    for r, c in product(row_vectors, col_vectors):
        value ^= determinant(matrix, r, c)
    return value


def check():
    cases = []
    for mixed in (False, True):
        p, l, _ = displayed_maps(mixed)
        cases.append((f'eight_letter_{mixed}', tuple(map(reverse_map, p)), tuple(map(reverse_map, l))))
    for extra in (False, True):
        tops = [(1 << 6) ^ (1 << 9)] + ([1 << 5] if extra else [])
        code = inverse_system(tops, 4)
        cases.append((f'inverse_system_extra_{extra}', tuple(code['embedding']), tuple(code['retraction'])))
    counts, records = Counter(), []
    for name, p, l in cases:
        computed = profiles(p, l, 4, keep_blocks=True)
        hh = [height(support(value)[0]) for value in p]
        literal = homogeneous_restriction(tuple(map(reverse_map, p)), tuple(map(reverse_map, l)), hh)
        assert computed['zero_orders_t'] == list(literal['zero_orders'])
        assert computed['infinity_orders_t'] == list(literal['infinity_orders'])
        counts['literal_cut_profile_comparisons'] += 2
        counts['literal_endpoint_entries'] += literal['literal_endpoint_entries']
        counts['literal_infinity_restoration_entries'] += literal['literal_infinity_restoration_entries']
        counts['deletion_space_controls'] += computed['deletion_checks']
        one = one_particle(4)
        for cut in computed['cuts']:
            for block in cut['blocks']:
                matrix = block.pop('matrix')
                chosen = sorted({0, len(matrix) // 2, len(matrix) - 1})
                for i, j in product(chosen, repeat=2):
                    assert matrix[i][j] == compressed_entry(one, 4, p, l, cut['cut'],
                                                            block['pairs'][i], block['pairs'][j])
                    counts['independent_complete_coefficient_checks'] += 1
                for precision in sorted({1, 4, 8, max(block['orders_t']) + 1}):
                    assert toeplitz_rank(matrix, precision) == sum(max(0, precision - order)
                                                                  for order in block['orders_t'])
                    counts['independent_quotient_rank_checks'] += 1
                counts['native_endpoint_entries'] += block['endpoint_entries']
        computed.pop('seconds')
        records.append({'name': name, 'result': computed,
                        'literal_raw_degree': literal['actual_polynomial_degree']})
    for modes, charge, energy in ((4, 2, 8), (8, 4, 32)):
        cell = [w for w in range(1 << modes) if w.bit_count() == charge and height(w) == energy]
        for mask in range(1, 1 << len(cell)):
            top = sum(1 << word for i, word in enumerate(cell) if mask >> i & 1)
            code = inverse_system([top], modes)
            p, l = code['embedding'], code['retraction']
            assert all((a & b).bit_count() % 2 == int(i == j)
                       for i, a in enumerate(l) for j, b in enumerate(p))
            chosen = set(code['pivot_words'])
            for word in chosen:
                for bit in range(modes):
                    if word >> bit & 1:
                        assert word ^ (1 << bit) in chosen
            for deleted in range(1 << modes):
                value = derivative(top, deleted)
                rebuilt = 0
                for column, row in zip(p, l):
                    if (row & value).bit_count() % 2:
                        rebuilt ^= column
                assert rebuilt == value
                counts['inverse_derivative_reconstructions'] += 1
            counts['complete_inverse_system_forms'] += 1
    invalid = [([1, 1 << 3], [1, 1 << 3], 4),
               ([1, 1 << 1], [1 << 1, 1], 4),
               ([1, (1 << 6) ^ (1 << 5)], [1, 1 << 6], 4),
               ([1, (1 << 7) ^ (1 << 16)], [1, 1 << 7], 5)]
    for args in invalid:
        try:
            validate_maps(*args)
        except ValueError:
            counts['invalid_map_rejections'] += 1
        else:
            raise AssertionError('An invalid coherent map passed the explicit hypotheses.')
    names = ('coherent_hot_compounds.py', 'coherent_hot_compounds.cpp', 'coherent_inverse_systems.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp',
             'eight_letter_reflection_retraction.py', Path(__file__).name)
    return {'status': 'pass', 'scope': __doc__, 'counts': dict(counts), 'records': records,
            'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                              for name in names}, 'factoring_exponent_improved': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = check()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'records'}, indent=2))


if __name__ == '__main__':
    main()
