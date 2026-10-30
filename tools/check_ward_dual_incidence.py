#!/usr/bin/env python3
"""Exact controls for both-sided incidence and maximal Ward charge descent."""

import argparse
from functools import lru_cache
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import random

from check_modular_crossing import rank_two
from check_noncommutative_crossing_radical_ward import local_clifford
from check_pfaffian_charge_gauges import E, N, exterior, kronecker
from complement_shell_tensor import multiply
from explore_complement_edge_pencils import evaluate_local
from pfaffian_particle_hole import particle_hole_mask


def bits(value):
    while value:
        bit = value & -value
        yield bit.bit_length() - 1
        value ^= bit


def determinant(rows):
    return int(rank_two(rows) == len(rows))


def generic_local_checks():
    stride = 8
    table = evaluate_local((2, 1 << stride, 1 << (stride * stride), 1, 1, 1))
    clifford = local_clifford()
    derivatives = []
    for value in table:
        derivative = 0
        for exponent in bits(value):
            for place in (1, stride, stride * stride):
                if exponent // place % stride & 1:
                    derivative ^= 1 << (exponent - place)
        derivatives.append(derivative)
    count = 0
    for x, y, u, v in product(range(4), repeat=4):
        index = ((x * 4 + y) * 4 + u) * 4 + v
        left = right = 0
        for mode in range(2):
            first_order = 0
            for state in bits(clifford[mode][u]):
                first_order ^= table[((x * 4 + y) * 4 + state) * 4 + v]
            inserted = 0
            if u >> mode & 1:
                inserted ^= table[((x * 4 + y) * 4 + (u ^ (1 << mode))) * 4 + v]
            if x >> (mode ^ 1) & 1:
                inserted ^= multiply(2, table[((((x ^ (1 << (mode ^ 1))) * 4 + y) * 4 + u) * 4 + v)])
            if y >> (mode ^ 1) & 1:
                inserted ^= table[((x * 4 + (y ^ (1 << (mode ^ 1)))) * 4 + u) * 4 + v]
            if v >> mode & 1:
                inserted ^= table[((x * 4 + y) * 4 + u) * 4 + (v ^ (1 << mode))]
            assert first_order == inserted
            if x >> (mode ^ 1) & 1:
                source = x ^ (1 << (mode ^ 1))
                for state in bits(clifford[mode][y]):
                    left ^= table[((source * 4 + state) * 4 + u) * 4 + v]
            if x >> mode & 1:
                source = x ^ (1 << mode)
                for state in bits(clifford[mode][u]):
                    right ^= table[((source * 4 + y) * 4 + state) * 4 + v]
        assert left == right == derivatives[index]
        count += 1
    return {'both_incidence_operators_equal_generic_derivative': count,
            'dual_first_order_identities': 2 * count}


def product_checks():
    rng = random.Random(302611)
    clifford = local_clifford()
    checked = 0
    for columns in ((1, 2), (1, 3), (1, 2, 3), (1, 2, 4), (3, 5, 6)):
        pairs, forms = len(columns), max(columns).bit_length()
        stride = 2 * pairs + 3
        table = evaluate_local((2, 1 << stride, 1 << (stride * stride), 1, 1, 1))
        derivative = []
        for value in table:
            result = 0
            for exponent in bits(value):
                for place in (1, stride, stride * stride):
                    if exponent // place % stride & 1:
                        result ^= 1 << (exponent - place)
            derivative.append(result)

        @lru_cache(None)
        def entry(x, y, u, v):
            result = 1
            for a, b, c, d in zip(x, y, u, v):
                result = multiply(result, table[((a * 4 + b) * 4 + c) * 4 + d])
            return result

        def acted(side, x, y, u, v):
            @lru_cache(None)
            def recurse(form, external, state):
                if form == forms:
                    return entry(external, state if side == 'Y' else y,
                                 state if side == 'U' else u, v)
                result = 0
                for site, column in enumerate(columns):
                    if not column >> form & 1:
                        continue
                    for mode in range(2):
                        target = mode ^ int(side == 'Y')
                        if external[site] >> target & 1:
                            new_external = external[:site] + (external[site] ^ (1 << target),) + external[site + 1:]
                            for value in bits(clifford[mode][state[site]]):
                                new_state = state[:site] + (value,) + state[site + 1:]
                                result ^= recurse(form + 1, new_external, new_state)
                return result
            return recurse(0, x, y if side == 'Y' else u)

        for _ in range(512):
            x, y, u, v = (tuple(rng.randrange(4) for _ in range(pairs)) for _ in range(4))
            expected = 0
            for selected in combinations(range(pairs), forms):
                rows = [sum(((columns[site] >> a) & 1) << j for j, site in enumerate(selected))
                        for a in range(forms)]
                if not determinant(rows):
                    continue
                value = 1
                for site, (a, b, c, d) in enumerate(zip(x, y, u, v)):
                    index = ((a * 4 + b) * 4 + c) * 4 + d
                    value = multiply(value, derivative[index] if site in selected else table[index])
                expected ^= value
            assert acted('Y', x, y, u, v) == expected
            assert acted('U', x, y, u, v) == expected
            checked += 1
    return checked


def gauge_checks():
    gauge = exterior(kronecker([[1, 0], [N, E]], [[E, N], [0, 1]]), scalar_parity=True)
    diagonal = off_diagonal = 0
    for a, y, b, z in product(range(4), repeat=4):
        value = gauge[particle_hole_mask(a, y)][particle_hole_mask(b, z)]
        if (a, y) == (b, z):
            assert value == (E if a.bit_count() == 1 else 1)
            diagonal += 1
        elif value:
            assert value == N
            assert a.bit_count() == b.bit_count() + 1
            assert abs(y.bit_count() - z.bit_count()) == 1
            off_diagonal += 1
    return {'diagonal_entries': diagonal, 'positive_shift_entries': off_diagonal}


def finite_difference_checks():
    rng = random.Random(617251)
    checked = 0
    for size in range(6):
        for shifts in range(size + 3):
            for _ in range(40):
                columns = [rng.randrange(1 << size) for _ in range(shifts)]
                offset = rng.randrange(1 << shifts)
                value = 0
                for choices in range(1 << shifts):
                    selected = choices ^ offset
                    rows = [0] * size
                    for i, column in enumerate(columns):
                        if selected >> i & 1:
                            for a in range(size):
                                if column >> a & 1:
                                    rows[a] ^= column
                    current = determinant(rows)
                    if shifts < size:
                        assert current == 0
                    value ^= current
                expected = determinant(columns) if shifts == size else 0
                assert value == expected
                checked += 1
    return checked


def laplace_checks():
    rng = random.Random(115722)
    checked = 0
    for size in range(7):
        for _ in range(80):
            columns = [rng.randrange(1 << size) for _ in range(size)]
            full_rows = [sum(((column >> i) & 1) << j for j, column in enumerate(columns))
                         for i in range(size)]
            for insertions in range(size + 1):
                value = 0
                for chosen_rows in combinations(range(size), insertions):
                    rest = [i for i in range(size) if i not in chosen_rows]
                    first = [sum(((columns[j] >> i) & 1) << j for j in range(insertions))
                             for i in chosen_rows]
                    second = [sum(((columns[j] >> i) & 1) << (j - insertions)
                                  for j in range(insertions, size)) for i in rest]
                    value ^= determinant(first) * determinant(second)
                assert value == determinant(full_rows)
                checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = {'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'generic_local_identities': generic_local_checks(),
              'multisite_product_derivative_controls': product_checks(),
              'exact_gauge_support': gauge_checks(),
              'maximal_shift_finite_difference_controls': finite_difference_checks(),
              'combined_creation_shift_minor_controls': laplace_checks()}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
