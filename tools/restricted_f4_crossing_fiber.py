#!/usr/bin/env python3
"""Project exact crossing-frame equations onto prescribed tensor directions."""
from __future__ import annotations

from f4_recentered_lattices import add_equation, assignment, dot, scale_row


def divide_one_plus_w(polynomial):
    quotient = 0
    while polynomial.bit_length() > 1:
        bit = 1 << (polynomial.bit_length()-2)
        quotient ^= bit
        polynomial ^= bit ^ (bit << 1)
    if polynomial:
        raise ValueError('A tensor direction must vanish at the strict endpoint.')
    return quotient


class RestrictedCrossingFiber:
    """The full linear span of supplied endpoint-zero polynomial directions.

    The parameter list is reduced to independent original directions over
    F4.  No new tensor directions are inferred from an approximate rank.
    This object can itself be passed to RecenteredFiber to compress common
    floors in the much smaller prescribed parameter space.
    """
    def __init__(self, lattice, directions):
        self.full_lattice = lattice
        self.q, self.size, self.degree = lattice.q, lattice.size, lattice.degree
        self.orders, self.caps = lattice.orders, lattice.caps
        self.floors = (0,)*len(self.orders)
        self.base = lattice.base
        self.source_direction_count = len(directions)
        self.selected_direction_indices = []
        self.coordinate_directions = []
        image_pivots = {}
        for index, direction in enumerate(directions):
            if len(direction) != self.q**4:
                raise ValueError('Each direction needs every literal tensor entry.')
            low = high = 0
            for entry, pair in enumerate(direction):
                if len(pair) != 2 or any(not isinstance(value, int) or value < 0 or value.bit_length() > self.degree+1
                                         for value in pair):
                    raise ValueError('Every direction must fit the same complete polynomial degree.')
                low |= divide_one_plus_w(pair[0]) << (entry*self.degree)
                high |= divide_one_plus_w(pair[1]) << (entry*self.degree)
            before = len(image_pivots)
            assert add_equation(image_pivots, low, high, 0)
            if len(image_pivots) > before:
                self.selected_direction_indices.append(index)
                self.coordinate_directions.append((low, high))
        self.variables = len(self.coordinate_directions)
        self.image_pivots = image_pivots
        self.equations = []
        for row in lattice.equations:
            projected = []
            for power, a, b, rhs in row:
                low = high = 0
                for column, direction in enumerate(self.coordinate_directions):
                    coefficient = dot(a, b, *direction)
                    low |= (coefficient & 1) << column
                    high |= (coefficient >> 1) << column
                if low or high or rhs:
                    projected.append((power, low, high, rhs))
            self.equations.append(projected)

    def solve(self, targets):
        if len(targets) != len(self.orders) or any(t < 0 or t > cap for t, cap in zip(targets, self.caps)):
            raise ValueError('Use targets inside the complete coefficient window.')
        pivots = {}
        for target, row in zip(targets, self.equations):
            for power, a, b, rhs in row:
                if power >= target:
                    break
                if not add_equation(pivots, a, b, rhs):
                    return None
        return pivots

    def seed(self, low, high):
        if low < 0 or high < 0 or max(low, high).bit_length() > self.variables:
            raise ValueError('Parameter values must fit the independent direction space.')
        original_low = original_high = 0
        for column, direction in enumerate(self.coordinate_directions):
            coefficient = ((low >> column) & 1) | (((high >> column) & 1) << 1)
            a, b = scale_row(*direction, coefficient)
            original_low ^= a
            original_high ^= b
        return self.full_lattice.seed(original_low, original_high)

    def sample(self, pivots, rng):
        free = ((1 << self.variables)-1) ^ sum(1 << pivot for pivot in pivots)
        parameters = assignment(pivots, rng.getrandbits(self.variables) & free,
                                rng.getrandbits(self.variables) & free)
        return self.seed(*parameters)

    def profiles(self, values):
        return self.full_lattice.profiles(values)

    def contains(self, values):
        if len(values) != len(self.base):
            return False
        low = high = 0
        try:
            for entry, (pair, base) in enumerate(zip(values, self.base)):
                if len(pair) != 2 or any(not isinstance(value, int) or value < 0 or value.bit_length() > self.degree+1
                                         for value in pair):
                    return False
                low |= divide_one_plus_w(pair[0] ^ base[0]) << (entry*self.degree)
                high |= divide_one_plus_w(pair[1] ^ base[1]) << (entry*self.degree)
        except ValueError:
            return False
        pivots = dict(self.image_pivots)
        assert add_equation(pivots, low, high, 0)
        return len(pivots) == self.variables

    def check_witness(self, values, targets):
        assert self.contains(values)
        return self.full_lattice.check_witness(values, targets)
