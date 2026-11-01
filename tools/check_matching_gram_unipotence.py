#!/usr/bin/env python3
"""Finite exact checks for the all-even-particle matching Gram theorem.

The all-size theorem is proved using Jack integrality and an integer
annihilator. Finite matrix checks here are corroboration, not extrapolation.
"""

import argparse
import json

from check_matching_code_lift import (
    cycle_type, gram_cycle_value, matching_partner,
)
from check_modular_crossing import column_bits, determinant_two, inverse_two, rank_two
from check_partition_channels import partitions
from fermionic_piece_polynomial import monomial_piece_polynomial, polynomial_product
from matching_code_lift import decode, labels, modes, permutation_action, specialize
from matching_homogeneous_lift import homogeneous_polynomial


def product_two(left, right):
    answer = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            row ^= bit
            value ^= right[bit.bit_length() - 1]
        answer.append(value)
    return answer


def power_two(matrix, exponent):
    result = [1 << i for i in range(len(matrix))]
    while exponent:
        if exponent & 1:
            result = product_two(result, matrix)
        exponent >>= 1
        if exponent:
            matrix = product_two(matrix, matrix)
    return result


def transpose_square(rows, size):
    """Transpose a binary square in O(size log size) Python integer steps."""
    power = 1 << (size - 1).bit_length()
    values = rows[:] + [0] * (power - len(rows))
    step = power // 2
    mask = (1 << step) - 1
    while step:
        for first in range(power):
            if first & step:
                continue
            change = ((values[first] >> step) ^ values[first + step]) & mask
            values[first] ^= change << step
            values[first + step] ^= change
        step //= 2
        if step:
            mask ^= mask << step
    return values[:size]


def column_blocks(rows, dimension):
    """The exact (x,u) | (v,y) flattening, via binary block transposes."""
    answer = [0] * dimension ** 2
    mask = (1 << dimension) - 1
    for x in range(dimension):
        for v in range(dimension):
            block = [(rows[x * dimension + y] >> (v * dimension)) & mask
                     for y in range(dimension)]
            for u, value in enumerate(transpose_square(block, dimension)):
                answer[x * dimension + u] |= value << (v * dimension)
    return answer


def collapsed_profile(vertices, direct_check=False):
    """Full exact ranks, using the stabilizer of one matching for entries."""
    labs = labels(vertices)
    dimension = len(labs)
    particles = vertices // 2
    generators = []
    for i in range(particles):
        permutation = list(range(vertices))
        permutation[2 * i], permutation[2 * i + 1] = permutation[2 * i + 1], permutation[2 * i]
        generators.append(permutation_action(labs, permutation))
    for i in range(particles - 1):
        permutation = list(range(vertices))
        permutation[2 * i], permutation[2 * i + 2] = permutation[2 * i + 2], permutation[2 * i]
        permutation[2 * i + 1], permutation[2 * i + 3] = permutation[2 * i + 3], permutation[2 * i + 1]
        generators.append(permutation_action(labs, permutation))
    seen = bytearray(dimension ** 3)
    base = [[0] * dimension for _ in range(2 * particles + 1)]
    orbit_count = 0
    partners = [matching_partner(label, vertices) for label in labs]
    for code in range(dimension ** 3):
        if seen[code]:
            continue
        y, v, u = code % dimension, (code // dimension) % dimension, code // dimension ** 2
        polynomial = homogeneous_polynomial(vertices, labs[0], labs[v], labs[u], labs[y], collapsed=True)
        assert polynomial.bit_count() % 2 == gram_cycle_value(partners[0], partners[v]) * int(y == u)
        degrees = [k for k in range(2 * particles + 1) if (polynomial >> k) & 1]
        pending = [(y, v, u)]
        seen[code] = 1
        while pending:
            y, v, u = pending.pop()
            for degree in degrees:
                base[degree][y] |= 1 << (v * dimension + u)
            for action in generators:
                target = action[y], action[v], action[u]
                target_code = target[0] + dimension * target[1] + dimension ** 2 * target[2]
                if not seen[target_code]:
                    seen[target_code] = 1
                    pending.append(target)
        orbit_count += 1
        if orbit_count % 500 == 0:
            polynomial_product.cache_clear()

    actions = []
    for label in labs:
        permutation = []
        for edge in label:
            permutation.extend(i for i in range(vertices) if (edge >> i) & 1)
        actions.append(permutation_action(labs, permutation))
    mask = (1 << dimension) - 1
    profile = []
    direct_entries = 0
    for degree, rows in enumerate(base):
        if not any(rows):
            continue
        full = [0] * dimension ** 2
        blocks = [[(row >> (v * dimension)) & mask for v in range(dimension)] for row in rows]
        distinct = set(value for row in blocks for value in row)
        for x, action in enumerate(actions):
            transformed = {}
            for value in distinct:
                image = 0
                remaining = value
                while remaining:
                    bit = remaining & -remaining
                    remaining ^= bit
                    image |= 1 << action[bit.bit_length() - 1]
                transformed[value] = image
            for y, row in enumerate(blocks):
                full[x * dimension + action[y]] = sum(
                    transformed[value] << (action[v] * dimension) for v, value in enumerate(row))
        columns = column_blocks(full, dimension)
        if direct_check:
            assert columns == column_bits(full, dimension)
            for x in range(dimension):
                for y in range(dimension):
                    for v in range(dimension):
                        for u in range(dimension):
                            polynomial = homogeneous_polynomial(vertices, labs[x], labs[v], labs[u], labs[y], collapsed=True)
                            assert ((full[x * dimension + y] >> (v * dimension + u)) & 1) == ((polynomial >> degree) & 1)
                            direct_entries += 1
        row_rank, column_rank = rank_two(full), rank_two(columns)
        profile.append([degree, row_rank, column_rank])
        print(json.dumps({'collapsed_vertices': vertices, 'degree': degree,
                          'row_rank': row_rank, 'column_rank': column_rank}), flush=True)
    return {'vertices': vertices, 'dimension': dimension, 'stabilizer_orbits': orbit_count,
            'profile': profile, 'maximum_rank_product': max(r * c for _, r, c in profile),
            'direct_entries': direct_entries}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--collapsed-eight', action='store_true')
    parser.add_argument('--output')
    args = parser.parse_args()
    reports = []
    determinant_checks = 0
    for particles in range(1, 6):
        vertices = 2 * particles
        labs = labels(vertices)
        partners = [matching_partner(label, vertices) for label in labs]
        gram = [sum(gram_cycle_value(first, second) << j for j, second in enumerate(partners))
                for first in partners]
        nilpotent = [row ^ ((1 << i) if particles % 2 == 0 else 0) for i, row in enumerate(gram)]
        bound = len(partitions(particles))
        assert not any(power_two(nilpotent, bound))
        if particles % 2 == 0:
            inverse = inverse_two(gram, len(labs))
            power = 1 << (bound - 1).bit_length()
            assert power_two(gram, power - 1) == inverse
        if particles <= 4:
            for i, x in enumerate(labs):
                for j, y in enumerate(labs):
                    matrix = [sum(int(not (e & f)) << column for column, f in enumerate(y)) for e in x]
                    assert determinant_two(matrix) == ((gram[i] >> j) & 1)
                    determinant_checks += 1
        reports.append({'particles': particles, 'dimension': len(labs),
                        'gram_rank': rank_two(gram), 'nilpotence_bound': bound})

    # k=2 is deliberately excluded from the even-eigenvalue statement.
    triangle = [6, 5, 3]
    assert product_two(triangle, triangle) == triangle and any(triangle)

    # G_6^2 is not I: compare a pair with cycle type (1,5).
    labs = labels(12)
    partners = [matching_partner(label, 12) for label in labs]
    first = partners[0]
    second = next(partner for partner in partners if cycle_type(first, partner) == (1, 5))
    off_diagonal_square = sum(gram_cycle_value(first, middle) * gram_cycle_value(middle, second)
                              for middle in partners) % 2
    assert off_diagonal_square == 1

    wick_checks = 0
    labs = labels(4)
    for x, v, u, y in __import__('itertools').product(labs, repeat=4):
        polynomial = homogeneous_polynomial(4, x, v, u, y)
        one_variable, offset = monomial_piece_polynomial(modes(x), modes(v), modes(u), modes(y))
        assert offset == -2 * sum(modes(x))
        assert specialize(4, polynomial) == one_variable
        collapsed = 0
        while polynomial:
            bit = polynomial & -polynomial
            polynomial ^= bit
            low, high = decode(4, bit.bit_length() - 1)
            exponent = 2 * low + sum(high)
            collapsed ^= 1 << exponent
        expected = homogeneous_polynomial(4, x, v, u, y, collapsed=True)
        assert all(not ((collapsed >> degree) & 1) for degree in range(collapsed.bit_length()) if degree % 8)
        assert sum(((collapsed >> (8 * degree)) & 1) << degree
                   for degree in range((collapsed.bit_length() + 7) // 8)) == expected
        wick_checks += 1

    collapsed = [collapsed_profile(4, direct_check=True)]
    assert collapsed[0]['maximum_rank_product'] == 16
    if args.collapsed_eight:
        collapsed.append(collapsed_profile(8))
        assert collapsed[-1]['maximum_rank_product'] == 7328 ** 2
    evidence = {
        'scope': 'finite exact checks; all-size unipotence is proved in the companion note',
        'gram_profiles': reports, 'determinant_checks': determinant_checks,
        'homogeneous_wick_checks': wick_checks,
        'negative_controls': {'A2_is_nonzero_idempotent': True,
                              'G6_squared_off_diagonal_cycle_type_1_5': off_diagonal_square},
        'collapsed_profiles': collapsed,
    }
    if args.output:
        with open(args.output, 'w') as handle:
            json.dump(evidence, handle, indent=2)
            handle.write('\n')
    print(json.dumps({key: value for key, value in evidence.items() if key != 'collapsed_profiles'}), flush=True)


if __name__ == '__main__':
    main()
