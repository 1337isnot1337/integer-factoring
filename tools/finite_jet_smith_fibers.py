#!/usr/bin/env python3
"""F4 finite Hasse equations inside a fixed two-cut regular-frame fiber.

The binary frame constructor is unchanged. Its endpoint-preserving affine
coordinates are extended to F4, then finite-place equations are reduced into
the same free coordinates. All polynomial entries remain present.
"""
from __future__ import annotations

from functools import lru_cache

from entrywise_reflection_constraints import field_dot, hasse_row
from full_alphabet_two_cut_lattices import CompressedFullAlphabetFiber
from physical_split_crossings import F4_PRODUCT
from solve_companion_affine_degree import coefficient, scaled


def insert(pivots, vector, rhs):
    a, b = vector
    while a or b:
        pivot = (a | b).bit_length() - 1
        value = coefficient((a, b), pivot)
        if pivot not in pivots:
            inverse = (0, 1, 3, 2)[value]
            pivots[pivot] = scaled((a, b), inverse), F4_PRODUCT[inverse][rhs]
            return True
        stored, answer = pivots[pivot]
        first, second = scaled(stored, value)
        a, b = a ^ first, b ^ second
        rhs ^= F4_PRODUCT[value][answer]
    return not rhs


class FiniteJetFiber:
    def __init__(self, lattice, floors):
        self.lattice = lattice
        self.binary = CompressedFullAlphabetFiber(lattice, floors)
        self.variables = self.binary.variables
        self.free = self.binary.free
        self.positions = {original: index for index, original in enumerate(self.free)}
        self.floors = self.binary.floors
        self.caps = lattice.caps
        self.source = tuple((value, 0) for value in lattice.base)

    def reduce_equation(self, vector, rhs):
        a, b = vector
        first = second = 0
        while a or b:
            pivot = (a | b).bit_length() - 1
            value = coefficient((a, b), pivot)
            if pivot in self.binary.base_pivots:
                row, constant = self.binary.base_pivots[pivot]
                low, high = scaled((row, 0), value)
                a, b = a ^ low, b ^ high
                rhs ^= value if constant else 0
            else:
                position = self.positions[pivot]
                first ^= (value & 1) << position
                second ^= (value >> 1) << position
                a ^= (value & 1) << pivot
                b ^= (value >> 1) << pivot
        return (first, second), rhs

    @lru_cache(maxsize=None)
    def jet_column(self, entry, point, order):
        degree = self.lattice.degree
        jet = hasse_row(point, order, degree)
        mask = (1 << degree) - 1
        # Delta p=(w+1)f with deg f<degree: coefficient k of f sees
        # Hasse(w^k)+Hasse(w^(k+1)).
        vector = tuple(((plane ^ (plane >> 1)) & mask) << (entry * degree) for plane in jet)
        source_value = field_dot(jet, self.source[entry])
        return self.reduce_equation(vector, source_value)

    def jet_equation(self, entry, point, order, target=0):
        vector, rhs = self.jet_column(entry, point, order)
        return vector, rhs ^ target

    def solve(self, equations, targets=None):
        targets = self.floors if targets is None else tuple(targets)
        if len(targets) != len(self.floors) or any(t < low or t > high for t, low, high in zip(targets, self.floors, self.caps)):
            raise ValueError('Every row target must lie in the constructed frame window.')
        pivots = {}
        for target, rows in zip(targets, self.binary.equations):
            for power, vector, rhs in rows:
                if power >= target:
                    break
                if not insert(pivots, (vector, 0), rhs):
                    return None
        for vector, rhs in equations:
            if not insert(pivots, vector, rhs):
                return None
        return pivots

    def increment(self, pivots, row, target):
        if target >= self.caps[row]:
            return None
        result = dict(pivots)
        for power, vector, rhs in self.binary.equations[row]:
            if power > target:
                break
            if power == target and not insert(result, (vector, 0), rhs):
                return None
        return result

    def sample(self, pivots, rng=None):
        result = (0, 0)
        if rng is not None:
            for index in range(self.variables):
                if index not in pivots:
                    scalar = rng.randrange(4)
                    result = (result[0] ^ ((scalar & 1) << index),
                              result[1] ^ ((scalar >> 1) << index))
        for pivot, (vector, rhs) in sorted(pivots.items()):
            value = field_dot(vector, result) ^ rhs
            result = (result[0] ^ ((value & 1) << pivot),
                      result[1] ^ ((value >> 1) << pivot))
        assert all(field_dot(vector, result) == rhs for vector, rhs in pivots.values())
        original = []
        for plane in range(2):
            remaining, vector = result[plane], 0
            while remaining:
                bit = remaining & -remaining
                remaining ^= bit
                vector ^= 1 << self.free[bit.bit_length() - 1]
            original.append(self.lattice.assignment(self.binary.base_pivots, vector, homogeneous=bool(plane)))
        return self.lattice.seed(*original)

    def verify_jets(self, values, requests):
        for entry, point, order, target in requests:
            assert field_dot(hasse_row(point, order, self.lattice.degree), values[entry]) == target
        for entry, value in enumerate(values):
            x, y, u, v = self.lattice.labels[entry]
            assert ((value[0].bit_count() & 1) | ((value[1].bit_count() & 1) << 1)) == int(x == v and y == u)
