#!/usr/bin/env python3
"""Exact F4 row-flag changes inside equal-order Smith blocks.

The reference tensor and its common-floor coefficient fiber are fixed.
Only the regular row frame changes.  Constant invertible changes within an
equal-order block preserve its common floor and therefore reuse the same
free affine coordinates exactly.
"""
from __future__ import annotations

from collections import defaultdict

from f4_recentered_lattices import (
    AffineCoordinateMap, add_equation, f4_rank, multiply_matrix,
    polynomial_multiply, scalar_multiply, scale_row, shifted, unpack,
)


def add_scaled(first, second, scalar):
    a, b = scale_row(second[0], second[1], scalar)
    return first[0] ^ a, first[1] ^ b, first[2] ^ scalar_multiply(second[2], scalar)


class SmithFlagChannels:
    """Column-resolved coefficient rows in the existing floor coordinates."""
    def __init__(self, lattice, fiber):
        if fiber.lattice is not lattice:
            raise ValueError('The coordinate fiber must belong to this lattice.')
        self.lattice, self.fiber = lattice, fiber
        self.variables = fiber.variables
        self.floors, self.caps, self.orders = fiber.floors, lattice.caps, lattice.orders
        self.blocks = []
        for side in range(2):
            blocks = defaultdict(list)
            for local, order in enumerate(self.orders[side*lattice.size:(side+1)*lattice.size]):
                blocks[order].append(side*lattice.size+local)
            for block in blocks.values():
                assert len({self.floors[i] for i in block}) == 1
                assert len({self.caps[i] for i in block}) == 1
                self.blocks.append(tuple(block))
        coordinates = AffineCoordinateMap(lattice.variables, fiber.base_pivots)
        assert coordinates.free == fiber.free
        self.channels = []
        variable_mask = (1 << lattice.degree)-1
        for side, (left, right) in enumerate(lattice.frames):
            scalar = lattice.scalar_shifts[side]
            h = lattice.heights[4*side:4*side+4]
            source = lattice.local_matrix(lattice.base, side)
            base = multiply_matrix(multiply_matrix(left, source, lattice.precisions[side]),
                                   right, lattice.precisions[side])
            start = side*lattice.size
            tables = [[[[0, 0] for _ in range(lattice.size)]
                       for power in range(self.floors[start+a], self.caps[start+a])]
                      for a in range(lattice.size)]
            for entry, bits in enumerate(lattice.labels):
                i, j = lattice.location(bits, side)
                offset = entry*lattice.degree
                weight = sum(row[label] for row, label in zip(h, bits))
                for a in range(lattice.size):
                    if not left[a][i]:
                        continue
                    floor, cap = self.floors[start+a], self.caps[start+a]
                    mask = (1 << (2*cap))-1
                    for b in range(lattice.size):
                        if not right[j][b]:
                            continue
                        outer = polynomial_multiply(left[a][i], right[j][b])
                        outer = (outer ^ (outer << 2)) & mask
                        for plane, coefficients in enumerate(unpack(outer)):
                            if not coefficients:
                                continue
                            if side:
                                width = coefficients.bit_length()
                                reverse = int(format(coefficients, '0'+str(width)+'b')[::-1], 2)
                                for power in range(floor, cap):
                                    at = power-scalar-weight
                                    section = shifted(reverse, at-width+1) & variable_mask
                                    tables[a][power-floor][b][plane] ^= section << offset
                            else:
                                for power in range(floor, cap):
                                    at = power-scalar-weight+1
                                    section = shifted(coefficients, -at) & variable_mask
                                    tables[a][power-floor][b][plane] ^= section << offset
            for a in range(lattice.size):
                floor, cap = self.floors[start+a], self.caps[start+a]
                row = []
                for power in range(floor, cap):
                    entries = []
                    for b in range(lattice.size):
                        low, high = tables[a][power-floor][b]
                        rhs = (base[a][b] >> (2*power)) & 3
                        entries.append(coordinates.reduce(low, high, rhs))
                    row.append(tuple(entries))
                self.channels.append(tuple(row))
        # Column information was independently reconstructed.  Removing it
        # must reproduce every existing reduced equation, in the same order.
        for row, (floor, channels) in enumerate(zip(self.floors, self.channels)):
            flattened = [(floor+power, a, b, c)
                         for power, entries in enumerate(channels) for a, b, c in entries if a or b or c]
            assert flattened == fiber.equations[row]

    def identity_mixes(self):
        n = self.lattice.size
        return tuple(tuple(tuple(int(i == j) for j in range(n)) for i in range(n)) for _ in range(2))

    def random_mixes(self, rng, sweeps=3):
        n = self.lattice.size
        result = [[[int(i == j) for j in range(n)] for i in range(n)] for _ in range(2)]
        for block in self.blocks:
            if len(block) < 2:
                continue
            side = block[0]//n
            indices = [index % n for index in block]
            for _ in range(sweeps*len(indices)):
                i, j = rng.sample(indices, 2)
                scalar = rng.randrange(1, 4)
                result[side][i] = [a ^ scalar_multiply(b, scalar)
                                   for a, b in zip(result[side][i], result[side][j])]
        return tuple(tuple(map(tuple, matrix)) for matrix in result)

    def mix(self, matrices):
        return MixedSmithFlagFiber(self, matrices)


class MixedSmithFlagFiber:
    """The same complete floor fiber, with different selected row conditions."""
    def __init__(self, channels, matrices):
        self.channels = channels
        self.lattice, self.fiber = channels.lattice, channels.fiber
        self.variables = channels.variables
        self.floors, self.caps = channels.floors, channels.caps
        self.matrices = tuple(tuple(map(tuple, matrix)) for matrix in matrices)
        n = self.lattice.size
        if len(matrices) != 2:
            raise ValueError('Use one row-frame change at each place.')
        self.equations = []
        self.frames = []
        for side, matrix in enumerate(matrices):
            if len(matrix) != n or any(len(row) != n or any(not isinstance(c, int) or not 0 <= c < 4 for c in row)
                                       for row in matrix):
                raise ValueError('A row-frame change must be an explicit square F4 matrix.')
            if f4_rank(matrix) != n:
                raise ValueError('The row-frame change must be invertible.')
            for i, row in enumerate(matrix):
                if any(c and channels.orders[side*n+i] != channels.orders[side*n+j] for j, c in enumerate(row)):
                    raise ValueError('Only equal-order blocks preserve the stated common floor.')
                floor, cap = self.floors[side*n+i], self.caps[side*n+i]
                equations = []
                for power in range(cap-floor):
                    for b in range(n):
                        equation = (0, 0, 0)
                        for j, coefficient in enumerate(row):
                            if coefficient:
                                equation = add_scaled(equation, channels.channels[side*n+j][power][b], coefficient)
                        if any(equation):
                            equations.append((floor+power, *equation))
                self.equations.append(equations)
            left, right = self.lattice.frames[side]
            self.frames.append((multiply_matrix(matrix, left), right))

    def solve(self, targets):
        if len(targets) != len(self.floors) or any(t < f or t > c for t, f, c in zip(targets, self.floors, self.caps)):
            raise ValueError('Every row target must remain inside its complete coefficient window.')
        pivots = {}
        for target, equations in zip(targets, self.equations):
            for power, low, high, rhs in equations:
                if power >= target:
                    break
                if not add_equation(pivots, low, high, rhs):
                    return None
        return pivots

    def sample(self, pivots, rng):
        return self.fiber.sample(pivots, rng)

    def check_witness(self, values, targets):
        self.lattice.validate_endpoint(values)
        if any(max(pair).bit_length() > self.lattice.degree+1 for pair in values):
            raise AssertionError('The witness exceeded its complete polynomial degree.')
        for side, (left, right) in enumerate(self.frames):
            source = self.lattice.local_matrix(values, side)
            precision = self.lattice.precisions[side]
            transformed = multiply_matrix(multiply_matrix(left, source, precision), right, precision)
            for i, row in enumerate(transformed):
                target = targets[side*self.lattice.size+i]
                assert all(not (entry & ((1 << (2*target))-1)) for entry in row)
        return self.lattice.profiles(values)
