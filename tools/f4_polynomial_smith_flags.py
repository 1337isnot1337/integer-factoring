#!/usr/bin/env python3
"""Polynomial Smith row transvections that reuse complete floor coordinates."""
from __future__ import annotations

from f4_recentered_lattices import multiply_matrix, scalar_multiply
from f4_smith_flag_fibers import MixedSmithFlagFiber, add_scaled


def random_transvections(channels, rng, sweeps=2, extra_degree=2):
    if sweeps < 1 or extra_degree < 0:
        raise ValueError('Use positive sweeps and a nonnegative extra degree.')
    n = channels.lattice.size
    operations = []
    for side in range(2):
        for _ in range(sweeps*n):
            i, j = rng.sample(range(n), 2)
            a, b = side*n+i, side*n+j
            shift = max(0, channels.orders[a]-channels.orders[b],
                        channels.floors[a]-channels.floors[b], channels.caps[a]-channels.caps[b])
            shift += rng.randrange(extra_degree+1)
            operations.append((side, i, j, shift, rng.randrange(1, 4)))
    return tuple(operations)


class PolynomialSmithFlagFiber(MixedSmithFlagFiber):
    """Apply exact row_i += scalar*z**shift*row_j operations to both jets.

    The shift preserves the reference row orders, the common-floor module,
    and availability of every coefficient below the recipient row cap.
    An omitted coefficient below a common floor is known to be zero for the
    entire coefficient fiber.  No unknown coefficient above a cap is used.
    """
    def __init__(self, channels, operations):
        self.channels = channels
        self.lattice, self.fiber = channels.lattice, channels.fiber
        self.variables = channels.variables
        self.floors, self.caps = channels.floors, channels.caps
        self.operations = tuple(tuple(operation) for operation in operations)
        n = self.lattice.size
        rows = [[list(entries) for entries in row] for row in channels.channels]
        matrices = [[[int(i == j) for j in range(n)] for i in range(n)] for _ in range(2)]
        for operation in self.operations:
            if len(operation) != 5 or any(not isinstance(value, int) for value in operation):
                raise ValueError('A transvection is (place, recipient, source, degree, F4 coefficient).')
            side, i, j, shift, coefficient = operation
            if side not in (0, 1) or not 0 <= i < n or not 0 <= j < n or i == j or coefficient not in (1, 2, 3):
                raise ValueError('Use distinct valid rows and a nonzero F4 coefficient.')
            a, b = side*n+i, side*n+j
            required = max(0, channels.orders[a]-channels.orders[b],
                           self.floors[a]-self.floors[b], self.caps[a]-self.caps[b])
            if shift < required:
                raise ValueError('The transvection does not preserve all stated row modules and windows.')
            for power in range(self.floors[a], self.caps[a]):
                at = power-shift
                if at < self.floors[b]:
                    continue
                assert at < self.caps[b]
                rows[a][power-self.floors[a]] = [add_scaled(first, second, coefficient)
                    for first, second in zip(rows[a][power-self.floors[a]], rows[b][at-self.floors[b]])]
            precision = self.lattice.precisions[side]
            mask = (1 << (2*precision))-1
            matrices[side][i] = [(first ^ (scalar_multiply(second, coefficient) << (2*shift))) & mask
                                 for first, second in zip(matrices[side][i], matrices[side][j])]
        self.matrices = tuple(tuple(map(tuple, matrix)) for matrix in matrices)
        self.equations = []
        self.coefficient_channels = tuple(tuple(map(tuple, row)) for row in rows)
        for floor, row in zip(self.floors, rows):
            self.equations.append([(floor+power, a, b, rhs)
                                   for power, entries in enumerate(row) for a, b, rhs in entries if a or b or rhs])
        self.frames = []
        for side, matrix in enumerate(matrices):
            left, right = self.lattice.frames[side]
            self.frames.append((multiply_matrix(matrix, left, self.lattice.precisions[side]), right))
