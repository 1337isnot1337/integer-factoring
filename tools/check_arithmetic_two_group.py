#!/usr/bin/env python3
"""Coordinate-independent 2-group algorithm, checked against known product groups."""
import json
from math import lcm, prod
import random

from arithmetic_two_group import power, project_to_two_part, two_group_basis


def rank_binary(rows):
    pivots = {}
    for row in rows:
        while row:
            bit = row.bit_length()-1
            if bit not in pivots:
                pivots[bit] = row
                break
            row ^= pivots[bit]
    return len(pivots)


def check_group(moduli, generators, annihilator=None):
    identity = (0,) * len(moduli)
    operations = [0]
    def multiply(a, b):
        operations[0] += 1
        return tuple((x+y) % modulus for x, y, modulus in zip(a, b, moduli))
    def inverse(a):
        return tuple(-x % modulus for x, modulus in zip(a, moduli))
    if annihilator is None:
        annihilator = lcm(*moduli) if moduli else 1
    projected = project_to_two_part(generators, annihilator, identity, multiply)
    exponents = [(modulus & -modulus).bit_length()-1 for modulus in moduli]
    basis = two_group_basis(projected, identity, multiply, inverse, max(exponents, default=0))
    assert sum(order for _, order in basis) == sum(exponents)
    tops = []
    for element, order in basis:
        assert power(element, 1 << order, identity, multiply) == identity
        top = power(element, 1 << (order-1), identity, multiply)
        assert top != identity
        mask = 0
        for i, (coordinate, modulus) in enumerate(zip(top, moduli)):
            assert coordinate in (0, modulus//2) if modulus % 2 == 0 else coordinate == 0
            if coordinate:
                mask |= 1 << i
        tops.append(mask)
    assert rank_binary(tops) == sum(exponent > 0 for exponent in exponents)
    return {'moduli': moduli, 'basis_order_logs': sorted(order for _, order in basis),
            'group_operations': operations[0], 'annihilator': annihilator}


def main():
    rng = random.Random(19921001)
    records = []
    # A single top power of each input misses a direction here; extraction
    # must form a combination of the two order-eight generators.
    mixed = check_group([8, 2], [(1, 0), (1, 1)], 64)
    assert mixed['basis_order_logs'] == [1, 3]
    swap = check_group([8, 2], [(0, 1), (1, 1)], 64)
    assert swap['basis_order_logs'] == [1, 3]
    for dimension in range(7):
        for trial in range(100):
            exponents = [rng.randrange(1, 7) for _ in range(dimension)]
            odd = [rng.choice((1, 3, 5, 9)) for _ in range(dimension)]
            moduli = [2**a*b for a, b in zip(exponents, odd)]
            generators = [tuple(int(i == j) for i in range(dimension)) for j in range(dimension)]
            # Integer elementary generator operations preserve the generated
            # group, including when the cyclic factor orders differ.
            if dimension > 1:
                for _ in range(12):
                    i, j = rng.sample(range(dimension), 2)
                    generators[i] = tuple((a+b) % modulus for a, b, modulus in zip(generators[i], generators[j], moduli))
            generators += [(0,) * dimension] * 2
            rng.shuffle(generators)
            exponent = lcm(*moduli) if moduli else 1
            # Deliberately oversaturate the 2-part and add odd factors: only
            # a multiple of the exponent is known, not the exact group order.
            annihilator = exponent * 2**rng.randrange(0, 9) * rng.choice((1, 3, 15))
            records.append(check_group(moduli, generators, annihilator))
    assert check_group([3, 5], [(1, 0), (0, 1)], 15)['basis_order_logs'] == []
    print(json.dumps({'scope': 'Exact generic 2-part basis extraction; no factoring runtime claim',
                      'product_groups_checked': len(records)+3,
                      'maximum_product_group_size': max(prod(row['moduli']) for row in records),
                      'maximum_basis_rank': max(len(row['basis_order_logs']) for row in records),
                      'maximum_group_operations': max(row['group_operations'] for row in records),
                      'mixed_generator_control': mixed, 'order_swap_control': swap,
                      'odd_group_control': 'passed'}, indent=2))


if __name__ == '__main__':
    main()
