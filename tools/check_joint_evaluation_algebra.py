#!/usr/bin/env python3
"""One good evaluation generates the complete five-letter energy algebra."""
from __future__ import annotations

import argparse
from collections import deque
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from f4_recentered_lattices import add_equation
from joint_energy_conjugations import (IDENTITY, MULTIPLICATION, N, Space, binary_product,
                                      field_product, kernel, model, unit_characters)

ROOT = Path(__file__).resolve().parents[1]


@lru_cache(maxsize=200000)
def polynomial_product(a, b):
    if a.bit_count() > b.bit_count():
        a, b = b, a
    result = 0
    while a:
        bit = a & -a
        result ^= b << (bit.bit_length()-1)
        a ^= bit
    return result


def polynomial_divmod(a, b):
    quotient, degree = 0, b.bit_length()
    assert b
    while a and a.bit_length() >= degree:
        shift = a.bit_length()-degree
        quotient ^= 1 << shift
        a ^= b << shift
    return quotient, a


def add_scaled(first, second, scalar):
    result = first.copy()
    for position, value in second.items():
        result[position] = result.get(position, 0) ^ polynomial_product(scalar, value)
        if not result[position]:
            del result[position]
    return result


def matrix_product(first, second):
    rows, result = {}, {}
    for position, value in second.items():
        rows.setdefault(position//N, []).append((position % N, value))
    for position, a in first.items():
        r, k = divmod(position, N)
        for c, b in rows.get(k, ()):
            position = N*r+c
            result[position] = result.get(position, 0) ^ polynomial_product(a, b)
            if not result[position]:
                del result[position]
    return result


def constant_algebra(data):
    masks = [sum(((1 << N)-1) << (N*r) for r in block) for block in data['blocks']]
    coefficients = [sum(((data['operator'][r][c][0] >> k) & 1) << (N*r+c)
                        for r in range(N) for c in range(N)) for k in range(19)]
    space, pending = Space([]), deque()
    def add(value):
        nonlocal space
        try:
            space.binary_coordinates(value)
        except KeyError:
            space = Space(space.basis+[value])
            pending.append(value)
    for block in data['blocks']:
        add(sum(1 << (N*r+r) for r in block))
    while pending:
        value = pending.popleft()
        for mask in masks:
            add(value & mask)
        for coefficient in coefficients:
            add(binary_product(coefficient, value))
    assert len(space.basis) == 169
    return space


def point_algebra(data, point):
    value = [0, 0]
    for r, row in enumerate(data['operator']):
        for c, (polynomial, high) in enumerate(row):
            assert not high
            entry, power = 0, 1
            for k in range(19):
                if polynomial >> k & 1:
                    entry ^= power
                power = MULTIPLICATION[power][point]
            for plane in range(2):
                value[plane] |= ((entry >> plane) & 1) << (N*r+c)
    pivots, pending = {}, deque()
    def add(matrix):
        count = len(pivots)
        assert add_equation(pivots, *matrix, 0)
        if len(pivots) != count:
            pending.append(matrix)
    masks = [sum(((1 << N)-1) << (N*r) for r in block) for block in data['blocks']]
    for block in data['blocks']:
        add((sum(1 << (N*r+r) for r in block), 0))
    while pending:
        matrix = pending.popleft()
        for mask in masks:
            add(tuple(v & mask for v in matrix))
        add(field_product(value, matrix))
    return len(pivots)


def polynomial_algebra(data):
    groups = data['blocks']
    generators = {}
    for i, rows in enumerate(groups):
        for j, columns in enumerate(groups):
            matrix = {N*r+c: data['operator'][r][c][0] for r in rows for c in columns
                      if data['operator'][r][c][0]}
            if matrix:
                generators[i, j] = matrix
    pivots, pending = {}, deque()
    insertions = 0
    def add(block, value):
        nonlocal insertions
        local = pivots.setdefault(block, {})
        while value:
            pivot = max(value)
            if pivot not in local:
                local[pivot] = value
                pending.append((block, value))
                insertions += 1
                return
            previous = local[pivot]
            quotient, remainder = polynomial_divmod(value[pivot], previous[pivot])
            value = add_scaled(value, previous, quotient)
            if remainder:
                assert value[pivot] == remainder
                local[pivot] = value
                pending.append((block, value))
                insertions += 1
                value = previous
    for i, block in enumerate(groups):
        add((i, i), {N*r+r: 1 for r in block})
    steps = 0
    while pending:
        (i, j), value = pending.popleft()
        for h in range(len(groups)):
            if (h, i) in generators:
                add((h, j), matrix_product(generators[h, i], value))
        steps += 1
    # Independent final module-membership replay of every left-generator product.
    def reduce(block, value):
        local = pivots.get(block, {})
        while value:
            pivot = max(value)
            assert pivot in local
            quotient, remainder = polynomial_divmod(value[pivot], local[pivot][pivot])
            assert not remainder
            value = add_scaled(value, local[pivot], quotient)
    closure_controls = 0
    for (i, j), local in pivots.items():
        for value in local.values():
            for h in range(len(groups)):
                if (h, i) in generators:
                    reduce((h, j), matrix_product(generators[h, i], value))
                    closure_controls += 1
    records = []
    x_total = y_total = 0
    for block, local in sorted(pivots.items()):
        for pivot, value in sorted(local.items()):
            polynomial, x, y = value[pivot], 0, 0
            while not polynomial & 1:
                polynomial >>= 1
                x += 1
            while polynomial.bit_count() % 2 == 0:
                polynomial, remainder = polynomial_divmod(polynomial, 3)
                assert not remainder
                y += 1
            assert polynomial == 1
            x_total += x
            y_total += y
            records.append(dict(block=list(block), pivot=pivot,
                                pivot_exponents=[x, y],
                                vector=[list(pair) for pair in sorted(value.items())]))
    assert len(records) == 169 and (x_total, y_total) == (1452, 350)
    return records, dict(insertions=insertions, processed_members=steps,
                         final_closure_controls=closure_controls,
                         triangular_minor_exponents=[x_total, y_total])


def certificate():
    data = model()
    full = constant_algebra(data)
    basis, statistics = polynomial_algebra(data)
    coefficient_controls = 0
    for record in basis:
        coefficients = {}
        for position, polynomial in record['vector']:
            while polynomial:
                bit = polynomial & -polynomial
                power = bit.bit_length()-1
                coefficients[power] = coefficients.get(power, 0) ^ (1 << position)
                polynomial ^= bit
        for value in coefficients.values():
            full.binary_coordinates(value)
            coefficient_controls += 1
    dimensions = [point_algebra(data, a) for a in range(4)]
    assert dimensions == [19, 11, 169, 169]
    # The full graded commutant has exactly four idempotents over every field.
    commutant = data['full'].basis
    assert all(binary_product(a, b) == binary_product(b, a) for a, b in product(commutant, repeat=2))
    _, characters = unit_characters(data)
    radical = [data['full'].matrix((v, 0))[0] for v in kernel(characters, 4)]
    assert len(radical) == 2 and all(binary_product(a, b) == 0 for a, b in product(radical, repeat=2))
    support = [0, 4, 6, 7, 8, 11, 12, 13, 16, 17, 18, 20, 24]
    first = sum(1 << (N*i+i) for i in support)
    second = IDENTITY ^ first
    assert binary_product(first, second) == binary_product(second, first) == 0
    assert len(Space(radical+[first, second]+commutant).basis) == 4
    assert len(Space(radical+[first, second]).basis) == 4
    return dict(coefficient_algebra_dimension=169, polynomial_basis=basis, statistics=statistics,
                coefficient_membership_controls=coefficient_controls, direct_F4_point_dimensions=dimensions,
                full_commutant_dimension=4, commutant_idempotents=[0, first, second, IDENTITY],
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('check_joint_evaluation_algebra.py', 'joint_energy_conjugations.py',
                                            'f4_recentered_lattices.py', 'search_f4_recentered_crossings.py')},
                scope='Every characteristic-two field and every finite point a other than zero and one; the energy projections and A(a) generate the complete constant coefficient algebra.',
                rational_gauge_conclusion='Every nonconstant scalar-rational projector gauge regular and invertible at zero, one and infinity that yields a polynomial transformed seed uses one of four neutral projectors and leaves the tensor unchanged.',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-evaluation-algebra.json')
    args = parser.parse_args()
    result = certificate()
    if args.generate:
        args.evidence.write_text(json.dumps(result, separators=(',', ':'))+'\n')
    else:
        assert json.loads(args.evidence.read_text()) == result
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('polynomial_basis', 'source_sha256', 'commutant_idempotents')}))


if __name__ == '__main__':
    main()
