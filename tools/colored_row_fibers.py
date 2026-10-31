#!/usr/bin/env python3
"""Exact F4 row-frame cells inside a common binary coefficient fiber.

The common floor is unchanged by every displayed elementary row operation.
Binary elimination therefore removes most variables before genuine F4
coefficient equations are assembled. No uncomputed coefficient is set to zero.
"""
from __future__ import annotations

from f4_lattice_linear import MUL, INV, scale, coefficient, dot
from full_alphabet_two_cut_lattices import cut_location
from weyl_lattice_filter_synthesis import polynomial_degrees, matrix_multiply


def insert_equation(pivots, low, high, rhs):
    """Insert an affine equation in place; return false exactly if inconsistent."""
    while low | high:
        pivot = (low | high).bit_length()-1
        value = ((low >> pivot) & 1) | (((high >> pivot) & 1) << 1)
        if pivot not in pivots:
            low, high = scale((low, high), INV[value])
            pivots[pivot] = (low, high, MUL[INV[value]][rhs])
            return True
        old_low, old_high, old_rhs = pivots[pivot]
        a, b = scale((old_low, old_high), value)
        low ^= a; high ^= b; rhs ^= MUL[value][old_rhs]
    return rhs == 0


def sample_assignment(pivots, variables, rng):
    free_mask = ((1 << variables)-1) ^ sum(1 << p for p in pivots)
    low = rng.getrandbits(variables) & free_mask
    high = rng.getrandbits(variables) & free_mask
    for pivot, (a, b, rhs) in sorted(pivots.items()):
        value = rhs ^ dot((a, b), (low, high))
        low |= (value & 1) << pivot
        high |= (value >> 1) << pivot
    return low, high


class ColoredRowFiber:
    def __init__(self, lattice, floors, operations, *, prepared=None):
        self.lattice = lattice
        self.floors = tuple(floors)
        self.operations = tuple(tuple(op) for op in operations)
        if len(floors) != len(lattice.orders):
            raise ValueError('Use one common floor for each row at both places.')
        if prepared is None:
            pivots = lattice.solve(floors)
            if pivots is None:
                raise ValueError('The common coefficient fiber is empty.')
            free = tuple(i for i in range(lattice.variables) if i not in pivots)
            positions = {old: new for new, old in enumerate(free)}
            rows = []
            for equations, columns in zip(lattice.equations, lattice.equation_columns):
                entries = {}
                for (power, value, rhs), column in zip(equations, columns):
                    reduced = 0
                    while value:
                        pivot = value.bit_length()-1
                        if pivot in pivots:
                            old, constant = pivots[pivot]
                            value ^= old; rhs ^= constant
                        else:
                            reduced ^= 1 << positions[pivot]; value ^= 1 << pivot
                    if reduced or rhs:
                        entries[power, column] = (reduced, 0, rhs)
                rows.append(entries)
            prepared = pivots, free, tuple(rows)
        self.prepared = prepared
        self.base_pivots, self.free, original_rows = prepared
        self.variables = len(self.free)
        rows = [dict(row) for row in original_rows]
        for side, i, j, power, scalar in self.operations:
            if side not in (0, 1) or i == j or not (0 <= i < lattice.size and 0 <= j < lattice.size):
                raise ValueError('Each elementary operation stays within one cut and has distinct rows.')
            if scalar not in (1, 2, 3) or power < 0:
                raise ValueError('Use a nonzero F4 coefficient and nonnegative polynomial power.')
            target, source = side*lattice.size+i, side*lattice.size+j
            if power < max(0, lattice.orders[target]-lattice.orders[source]):
                raise ValueError('The row operation must preserve the baseline row lattice.')
            if power < floors[target]-floors[source]:
                raise ValueError('The row operation must preserve the common floor lattice.')
            if lattice.caps[target]-power > lattice.caps[source]:
                raise ValueError('The coefficient windows do not cover the transformed target row.')
            for (exponent, column), (a, b, rhs) in rows[source].items():
                exponent += power
                if exponent >= lattice.caps[target]:
                    continue
                a, b = scale((a, b), scalar)
                rhs = MUL[scalar][rhs]
                old = rows[target].get((exponent, column), (0, 0, 0))
                new = old[0] ^ a, old[1] ^ b, old[2] ^ rhs
                if any(new):
                    rows[target][exponent, column] = new
                else:
                    rows[target].pop((exponent, column), None)
        self.equations = []
        for floor, row in zip(floors, rows):
            assert all(power >= floor for power, column in row)
            self.equations.append(tuple((power, *value) for (power, column), value in sorted(row.items())))

    def solve(self, targets):
        if len(targets) != len(self.floors) or any(t < f or t > c for t, f, c in
                zip(targets, self.floors, self.lattice.caps)):
            raise ValueError('The targets must stay in the completely constructed coefficient window.')
        pivots = {}
        for target, equations in zip(targets, self.equations):
            for power, low, high, rhs in equations:
                if power >= target:
                    break
                if not insert_equation(pivots, low, high, rhs):
                    return None
        return pivots

    def increment(self, pivots, row, target):
        result = dict(pivots)
        for power, low, high, rhs in self.equations[row]:
            if power < target:
                continue
            if power > target:
                break
            if not insert_equation(result, low, high, rhs):
                return None
        return result

    def sample(self, pivots, rng):
        planes = sample_assignment(pivots, self.variables, rng)
        lifted = []
        for plane, value in enumerate(planes):
            original = 0
            while value:
                bit = value & -value; value ^= bit
                original |= 1 << self.free[bit.bit_length()-1]
            lifted.append(self.lattice.assignment(self.base_pivots, original, homogeneous=bool(plane)))
        return self.lattice.seed(*lifted)

    def check_witness(self, values, targets):
        """Check divisibility and full, untruncated transformed tensor matrices."""
        from check_rectangular_defect_support import pdivmod
        model = self.lattice
        for entry, (low, high) in enumerate(values):
            assert not pdivmod(low ^ model.base[entry], model.factors[entry])[1]
            assert not pdivmod(high, model.factors[entry])[1]
            assert max(low, high).bit_length() <= model.degree+1
        for side, (left, right) in enumerate(model.bases):
            heights = model.heights[4*side:4*side+4]
            transformed = []
            for plane in range(2):
                matrix = [[0]*model.size for _ in range(model.size)]
                for entry, pair in enumerate(values):
                    i, j = cut_location(model.labels[entry], model.q, side)
                    weight = sum(h[label] for h, label in zip(heights, model.labels[entry]))
                    for power in polynomial_degrees(pair[plane]):
                        exponent = model.scalar_shifts[side]+weight+(-power if side == 0 else power)
                        assert exponent >= 0
                        matrix[i][j] ^= 1 << exponent
                transformed.append(matrix_multiply(matrix_multiply(left, matrix), right))
            for place, i, j, power, scalar in self.operations:
                if place != side:
                    continue
                for col in range(model.size):
                    a, b = scale((transformed[0][j][col], transformed[1][j][col]), scalar)
                    transformed[0][i][col] ^= a << power
                    transformed[1][i][col] ^= b << power
            for plane in transformed:
                for i, row in enumerate(plane):
                    mask = (1 << targets[side*model.size+i])-1
                    assert all(value & mask == 0 for value in row)
        profiles = model.profiles(values)
        for side in range(2):
            correction = model.clearing+model.pole if side == 0 else -model.clearing
            bounds = sorted(targets[side*model.size+i]-model.scalar_shifts[side]+correction
                            for i in range(model.size))
            for k in range(1, model.size+1):
                assert sum(profiles[side][:k]) >= sum(bounds[:k])
        return profiles
