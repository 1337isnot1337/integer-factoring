#!/usr/bin/env python3
"""Exact joint-energy conjugation algebra of the retained five-letter crossing."""
from __future__ import annotations

from functools import lru_cache
from itertools import product

from check_controlled_permutation_smith import binary_rank
from search_f4_recentered_crossings import five_letter_source

N = 25
MASK = (1 << N)-1
MULTIPLICATION = ((0, 0, 0, 0), (0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2))
INVERSE = (0, 1, 3, 2)
IDENTITY = sum(1 << (N*i+i) for i in range(N))


def binary_product(a, b):
    rows = [(b >> (N*i)) & MASK for i in range(N)]
    result = 0
    for i in range(N):
        row, value = (a >> (N*i)) & MASK, 0
        while row:
            bit = row & -row
            value ^= rows[bit.bit_length()-1]
            row ^= bit
        result |= value << (N*i)
    return result


def scale(value, scalar):
    a, b = value
    return value if scalar == 1 else (b, a ^ b) if scalar == 2 else (a ^ b, a) if scalar == 3 else (0, 0)


def field_product(a, b):
    x, y = a
    z, w = b
    yw = binary_product(y, w)
    return (binary_product(x, z) ^ yw,
            binary_product(x, w) ^ binary_product(y, z) ^ yw)


class Space:
    def __init__(self, vectors):
        self.basis, self.pivots = [], {}
        for vector in vectors:
            value, coordinates = vector, 0
            while value:
                pivot = value.bit_length()-1
                if pivot not in self.pivots:
                    self.pivots[pivot] = value, coordinates ^ (1 << len(self.basis))
                    self.basis.append(vector)
                    break
                a, b = self.pivots[pivot]
                value, coordinates = value ^ a, coordinates ^ b

    def binary_coordinates(self, vector):
        coordinates = 0
        while vector:
            pivot = vector.bit_length()-1
            a, b = self.pivots[pivot]
            vector, coordinates = vector ^ a, coordinates ^ b
        return coordinates

    def coordinates(self, matrix):
        return tuple(self.binary_coordinates(v) for v in matrix)

    def matrix(self, coordinates):
        result = [0, 0]
        for plane, bits in enumerate(coordinates):
            while bits:
                bit = bits & -bits
                result[plane] ^= self.basis[bit.bit_length()-1]
                bits ^= bit
        return tuple(result)


def kernel(equations, dimension):
    pivots = {}
    for value in equations:
        while value:
            pivot = value.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    basis = []
    for free in range(dimension):
        if free in pivots:
            continue
        vector = 1 << free
        for pivot in sorted(pivots):
            if (pivots[pivot] & vector).bit_count() % 2:
                vector ^= 1 << pivot
        assert all(not ((row & vector).bit_count() % 2) for row in pivots.values())
        basis.append(vector)
    return basis


def commutator_equations(matrix, positions):
    equations = [0]*(N*N)
    for variable, (i, j) in enumerate(positions):
        elementary = 1 << (N*i+j)
        value = binary_product(elementary, matrix) ^ binary_product(matrix, elementary)
        while value:
            bit = value & -value
            equations[bit.bit_length()-1] ^= 1 << variable
            value ^= bit
    return equations


@lru_cache(None)
def model():
    values, heights, energy_x, energy_y = five_letter_source()
    assert energy_x == energy_y
    half = tuple(v//2 for v in energy_x)
    energies = [half[i % 5]+half[i//5] for i in range(N)]
    blocks = {}
    for i, energy in enumerate(energies):
        blocks.setdefault(energy, []).append(i)
    positions = [(i, j) for i in range(N) for j in range(N) if energies[i] == energies[j]]
    weighted = {}
    operator = [[(0, 0)]*N for _ in range(N)]
    for index, (p, high) in enumerate(values):
        assert high == 0
        x, y, u, v = ((index//5**i) % 5 for i in range(4))
        row, column = x+5*y, v+5*u
        operator[row][column] = (p, 0)
        for k in range(19):
            if p >> k & 1:
                power = energies[column]-k
                weighted[power] = weighted.get(power, 0) ^ (1 << (N*row+column))
    powers = sorted(weighted)
    assert powers[:2] == [-8, -7]
    equations = commutator_equations(weighted[-8], positions)
    first = kernel(equations, len(positions))
    second = kernel(equations+commutator_equations(weighted[-7], positions), len(positions))
    matrix = lambda vector: sum(1 << (N*i+j) for k, (i, j) in enumerate(positions) if vector >> k & 1)
    leading, full = list(map(matrix, first)), list(map(matrix, second))
    assert len(positions) == 69 and len(leading) == 10 and len(full) == 4
    assert all(binary_product(a, b) == binary_product(b, a) for a in full for b in weighted.values())
    return dict(values=values, heights=heights, half=half, energies=energies,
                blocks=list(blocks.values()), positions=positions, weighted=weighted,
                operator=operator, leading=Space(leading), full=Space(full))


def inverse(matrix, blocks):
    result = [0, 0]
    for block in blocks:
        d = len(block)
        a = [[sum(((matrix[b] >> (N*r+c)) & 1) << b for b in range(2)) for c in block]
             + [int(r == c) for c in block] for r in block]
        for j in range(d):
            pivot = next((i for i in range(j, d) if a[i][j]), None)
            if pivot is None:
                return None
            a[j], a[pivot] = a[pivot], a[j]
            coefficient = INVERSE[a[j][j]]
            a[j] = [MULTIPLICATION[coefficient][v] for v in a[j]]
            for i in range(d):
                if i != j:
                    coefficient = a[i][j]
                    a[i] = [v ^ MULTIPLICATION[coefficient][w] for v, w in zip(a[i], a[j])]
        for i, row in enumerate(block):
            for j, column in enumerate(block):
                for b in range(2):
                    result[b] |= ((a[i][d+j] >> b) & 1) << (N*row+column)
    return tuple(result)


def apply(columns, vector):
    low, high = vector
    result, bits = [0, 0], low | high
    while bits:
        bit = bits & -bits
        i = bit.bit_length()-1
        coefficient = ((low >> i) & 1) | (((high >> i) & 1) << 1)
        a, b = scale(columns[i], coefficient)
        result[0], result[1] = result[0] ^ a, result[1] ^ b
        bits ^= bit
    return tuple(result)


def orbit(data):
    algebra = data['leading']
    coefficient = data['weighted'][-7]
    module = Space(binary_product(binary_product(a, coefficient), b)
                   for a, b in product(algebra.basis, repeat=2))
    assert len(module.basis) == 16
    generators = []
    for basis in algebra.basis:
        for scalar in (1, 2, 3):
            a, b = scale((basis, 0), scalar)
            unit = (a ^ IDENTITY, b)
            inv = inverse(unit, data['blocks'])
            if inv is None:
                continue
            assert field_product(unit, inv) == (IDENTITY, 0)
            action = [module.coordinates(field_product(field_product(unit, (m, 0)), inv))
                      for m in module.basis]
            left = [algebra.coordinates(field_product(unit, (m, 0))) for m in algebra.basis]
            generators.append((action, left))
    initial = module.coordinates((coefficient, 0))
    seen = {initial: algebra.coordinates((IDENTITY, 0))}
    pending = [initial]
    while pending:
        state = pending.pop()
        for action, left in generators:
            target = apply(action, state)
            if target not in seen:
                seen[target] = apply(left, seen[state])
                pending.append(target)
    return list(seen.values()), len(generators)


def binary_inverse(rows):
    n = len(rows)
    a = [r | (1 << (n+i)) for i, r in enumerate(rows)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i] >> j & 1)
        a[j], a[pivot] = a[pivot], a[j]
        for i in range(n):
            if i != j and a[i] >> j & 1:
                a[i] ^= a[j]
    return [row >> n for row in a]


def unit_characters(data):
    characters, full_characters = set(), set()
    for block in data['blocks']:
        d = len(block)
        submatrix = lambda m: [sum(((m >> (N*r+c)) & 1) << j for j, c in enumerate(block)) for r in block]
        matrices = list(map(submatrix, data['leading'].basis))
        def flag(prefix):
            if len(prefix) == d:
                return prefix
            for vector in range(1, 1 << d):
                if binary_rank(prefix+[vector]) == len(prefix):
                    continue
                next_prefix = prefix+[vector]
                images = [sum(((row & vector).bit_count() % 2) << i for i, row in enumerate(a)) for a in matrices]
                if all(binary_rank(next_prefix+[v]) == len(next_prefix) for v in images):
                    result = flag(next_prefix)
                    if result:
                        return result
            return None
        columns = flag([])
        assert columns is not None
        change = [sum(((v >> i) & 1) << j for j, v in enumerate(columns)) for i in range(d)]
        inv = binary_inverse(change)
        def small_product(a, b):
            return [sum(((row & sum(((b[k] >> j) & 1) << k for k in range(d))).bit_count() % 2) << j
                        for j in range(d)) for row in a]
        def transform(a):
            return small_product(small_product(inv, a), change)
        transformed = list(map(transform, matrices))
        full = [transform(submatrix(a)) for a in data['full'].basis]
        assert all(not (a[i] & ((1 << i)-1)) for a in transformed for i in range(d))
        for i in range(d):
            characters.add(sum(((a[i] >> i) & 1) << j for j, a in enumerate(transformed)))
            full_characters.add(sum(((a[i] >> i) & 1) << j for j, a in enumerate(full)))
    assert len(characters) == binary_rank(list(characters)) == 6
    assert len(full_characters) == binary_rank(list(full_characters)) == 2
    return sorted(characters), sorted(full_characters)


def conjugated_values(data, coordinates):
    unit = data['leading'].matrix(coordinates)
    inv = inverse(unit, data['blocks'])
    assert inv is not None
    coefficient = lambda a, r, c: sum(((a[b] >> (N*r+c)) & 1) << b for b in range(2))
    temp, output = [[(0, 0)]*N for _ in range(N)], [[(0, 0)]*N for _ in range(N)]
    for block in data['blocks']:
        for r, k in product(block, repeat=2):
            scalar = coefficient(unit, r, k)
            if scalar:
                for c in range(N):
                    a, b = scale(data['operator'][k][c], scalar)
                    temp[r][c] = temp[r][c][0] ^ a, temp[r][c][1] ^ b
    for block in data['blocks']:
        for k, c in product(block, repeat=2):
            scalar = coefficient(inv, k, c)
            if scalar:
                for r in range(N):
                    a, b = scale(temp[r][k], scalar)
                    output[r][c] = output[r][c][0] ^ a, output[r][c][1] ^ b
    return tuple(output[x+5*y][v+5*u] for v in range(5) for u in range(5) for y in range(5) for x in range(5))
