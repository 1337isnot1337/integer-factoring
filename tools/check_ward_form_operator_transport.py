#!/usr/bin/env python3
"""Form-variable insertion transitions in the eight-dimensional Ward algebra."""

import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path
import random

from check_modular_crossing import rank_two
from check_ward_stable_operators import (
    apply, bits, closed_operator_interface, conjugate, inverse, local_operators,
    stable_columns, tensor_product, unflatten,
)


def insertion_transitions():
    interface = closed_operator_interface()
    basis = interface['homogeneous_basis_matrices']
    columns = stable_columns()
    pivots = {}
    for j, value in enumerate(basis):
        coordinate = 1 << j
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value, coordinate
                break
            old, old_coordinate = pivots[pivot]
            value ^= old
            coordinate ^= old_coordinate
    def coordinate(value):
        result = 0
        while value:
            old, old_coordinate = pivots[value.bit_length() - 1]
            value ^= old
            result ^= old_coordinate
        return result
    def decompose(operator):
        result = [0] * len(basis)
        for i in range(16):
            for j in range(16):
                coefficient = sum(((operator[16 * b + j] >> (16 * a + i)) & 1) << (4 * a + b)
                                  for a in range(4) for b in range(4))
                for beta in bits(coordinate(coefficient)):
                    result[beta] ^= 1 << (16 * i + j)
        rebuilt = [0] * 64
        for a, b in zip(basis, result):
            for j, value in enumerate(tensor_product(unflatten(a, 4), unflatten(b, 16))):
                rebuilt[j] ^= value
        assert rebuilt == operator
        return result
    rs, os, _, _, _, _ = local_operators()
    zero_checks = 0
    for first, second in ((rs, rs), (os, os), (rs, os), (os, rs)):
        assert all(apply(first, column) == 0 for column in second)
        zero_checks += len(second)
    transition_r = [decompose(conjugate(tensor_product(rs, unflatten(a, 4)), columns)) for a in basis]
    transition_o = [decompose(conjugate(tensor_product(os, unflatten(a, 4)), columns)) for a in basis]
    doubled = []
    for original in range(1024):
        old, new = original // 16, original % 16
        result = 0
        for row in bits(columns[old]):
            physical, ancillary = row // 4, row % 4
            for last in bits(columns[16 * ancillary + new]):
                result ^= 1 << (64 * physical + last)
        doubled.append(result)
    reverse = inverse(doubled)
    count = 0
    for first in (rs, os):
        earlier = [conjugate(tensor_product(first, unflatten(a, 4)), columns) for a in basis]
        for second, emission in ((rs, interface['incidence_emission']), (os, interface['creation_emission'])):
            expected = [0] * 1024
            for operator, coefficient in zip(earlier, emission):
                for j, value in enumerate(tensor_product(operator, unflatten(coefficient, 16))):
                    expected[j] ^= value
            source = tensor_product(tensor_product(first, second), [1 << i for i in range(4)])
            actual = [apply(reverse, apply(source, column)) for column in doubled]
            assert actual == expected
            count += len(actual)
    return {'ancillary_algebra_dimension': len(basis),
            'physical_insertion_transitions_incidence': transition_r,
            'physical_insertion_transitions_creation': transition_o,
            'local_square_zero_and_mixed_zero_columns': zero_checks,
            'two_site_ordered_insertion_columns': count,
            'transported_insertion_operator_matrices': 2 * len(basis)}


def form_minor_checks():
    rng = random.Random(320711)
    rs, os, _, _, _, _ = local_operators()
    def local_action(vector, operator, site):
        result = 0
        for source in bits(vector):
            digit = source >> (4 * site) & 15
            base = source ^ (digit << (4 * site))
            for target in bits(operator[digit]):
                result ^= 1 << (base | (target << (4 * site)))
        return result
    checked = 0
    for columns in ((1, 1), (1, 2), (1, 2, 3), (1, 2, 4), (3, 5, 6)):
        pairs, forms = len(columns), max(columns).bit_length()
        for _ in range(256):
            source = rng.randrange(1 << (4 * pairs))
            actual = [1 << source]
            for a in range(forms):
                following = [0] * (len(actual) + 1)
                for degree, vector in enumerate(actual):
                    for site, column in enumerate(columns):
                        if column >> a & 1:
                            following[degree] ^= local_action(vector, rs, site)
                            following[degree + 1] ^= local_action(vector, os, site)
                actual = following
            expected = [0] * (forms + 1)
            for selected in combinations(range(pairs), forms):
                if rank_two([columns[i] for i in selected]) != forms:
                    continue
                value = [1 << source]
                for site in selected:
                    following = [0] * (len(value) + 1)
                    for degree, vector in enumerate(value):
                        following[degree] ^= local_action(vector, rs, site)
                        following[degree + 1] ^= local_action(vector, os, site)
                    value = following
                expected = [a ^ b for a, b in zip(expected, value)]
            assert actual == expected
            checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in (Path(__file__).name, 'check_ward_stable_operators.py')},
              'insertion_transport': insertion_transitions(),
              'form_product_vs_site_minor_controls': form_minor_checks()}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
