#!/usr/bin/env python3
"""Independent path enumeration and complete-profile checks for closed blocks."""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_block_smith_crossings import smith_profile, toeplitz_rank, reverse_polynomial
from closed_auxiliary_blocks import tensor, flatten, local_entry, product_masks
from native_rectangular_smith_f4 import RectangularOracle


def literal_product(polynomials):
    coefficients = [1]
    for polynomial in polynomials:
        next_coefficients = [0]*(len(coefficients)+max(0, polynomial.bit_length()-1))
        for i, a in enumerate(coefficients):
            for j in range(polynomial.bit_length()):
                next_coefficients[i+j] ^= a*((polynomial >> j) & 1)
        coefficients = next_coefficients
    return sum(a << i for i, a in enumerate(coefficients))


def literal_entry(length, labels, mask_b, mask_c):
    answer = 0
    for internal in product((0, 1), repeat=length-1):
        path = (0,)+internal+(0,)
        terms = []
        for i in range(length):
            x, y, u, v = ((a >> i) & 1 for a in labels)
            matrix = local_entry(x, y, u, v, mask_b, mask_c)
            terms.append(matrix[path[i]][path[i+1]])
        answer ^= literal_product(terms)
    return answer


def run():
    counts = Counter(); records = []
    cases = ((2, 0, 0, None), (2, 0x6996, 0x1248, None),
             (3, 0x318c, 0xc33c, None), (3, 0x6996, 0x1248, None),
             (4, 0x18e7, 0xa5a5, (0, 3, 5, 10, 12)))
    for length, b, c, code in cases:
        data = tensor(length, b, c, code)
        q = data['dimensions'][0]; labels = data['physical_words']
        for index, (value, high) in enumerate(data['values']):
            x, y, u, v = ((index//q**j) % q for j in range(4))
            source = tuple(labels[a] for a in (x, y, u, v))
            assert value == literal_entry(length, source, b, c) and high == 0
            assert value.bit_count() % 2 == int(x == v and y == u)
            counts['literal_path_entries'] += 1
            counts['literal_endpoint_entries'] += 1
        oracle = RectangularOracle(data['values'], q, q)
        references = []
        for cut in (0, 1):
            matrix = flatten(data['values'], q, cut)
            reverse = [[reverse_polynomial(p, length) for p in row] for row in matrix]
            zero, _ = smith_profile(matrix, length)
            infinity, _ = smith_profile(reverse, length)
            references.append((zero, tuple(v-length for v in infinity)))
            for polynomial, orders in ((matrix, zero), (reverse, infinity)):
                for precision in sorted({1, 2, max(orders)+1}):
                    assert toeplitz_rank(polynomial, precision) == sum(max(0, precision-v) for v in orders)
                    counts['independent_toeplitz_ranks'] += 1
        for cut in (0, 1):
            _, _, orders = oracle.evaluate([0.]*(oracle.parameters-1)+[1.],
                                            [1/q]*q, [1/q]*q, cut, True)
            assert tuple(orders[:q*q]) == references[cut][1]
            assert tuple(orders[q*q:]) == references[1-cut][0]
            counts['native_complete_profile_pairs'] += 1
        records.append({'length': length, 'physical_words': labels, 'mask_b': b,
                        'mask_c': c, 'profiles_zero_then_infinity': references})
    # Transposition followed by reversal of the word exchanges the two masks.
    for b, c in ((0x6996, 0x1248), (0x18e7, 0xa5a5)):
        first = tensor(3, b, c); second = tensor(3, c, b)
        reverse_word = lambda x: ((x & 1) << 2) | (x & 2) | ((x >> 2) & 1)
        for i, value in enumerate(first['values']):
            coordinates = [reverse_word((i//8**j) % 8) for j in range(4)]
            other = sum(a*8**j for j, a in enumerate(coordinates))
            assert value == second['values'][other]
            counts['mask_exchange_reversal_entries'] += 1
    root = Path(__file__).resolve().parent
    return {'scope': 'Exact closed auxiliary endpoints, path compiler, and finite complete-profile controls; no rate improvement.',
            'checks': dict(counts), 'cases': records, 'mask_catalogue_size': len(product_masks()),
            'source_sha256': {name: sha256((root/name).read_bytes()).hexdigest()
                              for name in ('closed_auxiliary_blocks.py', Path(__file__).name)}}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
