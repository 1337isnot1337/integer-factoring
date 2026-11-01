#!/usr/bin/env python3
"""Complete-entry crossing fibers in regular frames over F4.

The reference tensor may have arbitrary F4 coefficients.  Every entry is
varied by (w+1)f(w), with deg(f)<degree, so the literal endpoint is fixed.
Row divisibility is imposed after exact regular frame changes on the two
cuts.  A fiber is a sufficient, explicitly computable family; it does not
parameterize every tensor with the same determinantal lower bounds.
"""
from __future__ import annotations

from itertools import product

from entrywise_reflection_constraints import field_dot, hasse_row
from f4_regular_smith_frames import bounded_frame, matrix_product, polynomial_product
from finite_jet_smith_fibers import insert
from full_alphabet_two_cut_lattices import cut_location, shifted
from native_matrix_smith_f4 import MatrixOracle
from physical_split_crossings import F4_PRODUCT
from solve_companion_affine_degree import coefficient, scaled
from weyl_lattice_filter_synthesis import polynomial_degrees


def reduce_row(pivots, vector, rhs):
    """Reduce an F4 equation, retaining its free original coordinates."""
    a, b = vector
    remaining = (0, 0)
    while a or b:
        pivot = (a | b).bit_length() - 1
        scalar = coefficient((a, b), pivot)
        if pivot in pivots:
            row, answer = pivots[pivot]
            low, high = scaled(row, scalar)
            a, b = a ^ low, b ^ high
            rhs ^= F4_PRODUCT[scalar][answer]
        else:
            low, high = (scalar & 1) << pivot, (scalar >> 1) << pivot
            remaining = remaining[0] ^ low, remaining[1] ^ high
            a, b = a ^ low, b ^ high
    return remaining, rhs


def back_substitute(pivots, value=(0, 0)):
    for pivot, (row, rhs) in sorted(pivots.items()):
        assert not ((value[0] | value[1]) >> pivot) & 1
        scalar = field_dot(row, value) ^ rhs
        value = value[0] ^ ((scalar & 1) << pivot), value[1] ^ ((scalar >> 1) << pivot)
    return value


class F4FullAlphabetLattice:
    def __init__(self, values, q, degree, heights, raise_cap=8):
        if q < 2 or degree < 1 or raise_cap < 1 or len(values) != q**4 or len(heights) != 8:
            raise ValueError('Use a full tensor, eight physical heights, and positive degree/window.')
        if any(len(row) != q or any(not isinstance(value, int) for value in row) for row in heights):
            raise ValueError('Physical heights must be integer lists of length q.')
        if any(len(pair) != 2 or any(not isinstance(value, int) or value < 0 or value.bit_length() > degree + 1
                                    for value in pair) for pair in values):
            raise ValueError('F4 entries use two nonnegative coefficient planes within the declared degree.')
        self.q, self.size, self.degree = q, q*q, degree
        self.variables = q**4 * degree
        self.base = tuple(map(tuple, values))
        self.heights = tuple(map(tuple, heights))
        self.labels = tuple(tuple((entry // q**i) % q for i in range(4)) for entry in range(q**4))
        self.check_endpoint(self.base)
        self.bases, self.scalar_shifts, self.orders, self.caps, self.equations = [], [], [], [], []
        mask = (1 << degree) - 1
        for side in range(2):
            h = self.heights[4*side:4*side+4]
            scalar = -sum(min(row) for row in h) + (degree if side == 0 else 0)
            weights = [sum(row[label] for row, label in zip(h, bits)) for bits in self.labels]
            source = self.weighted_matrix(self.base, side, scalar)
            orders, left, right, diagonal = bounded_frame(source, raise_cap)
            self.bases.append((left, right))
            self.scalar_shifts.append(scalar)
            self.orders.extend(orders)
            caps = [order + raise_cap for order in orders]
            self.caps.extend(caps)
            # Each table entry is one F4 linear equation in all polynomial
            # coefficients.  Coefficient bitplanes use the same variable order.
            low_tables = [[[0] * cap for _ in range(self.size)] for cap in caps]
            high_tables = [[[0] * cap for _ in range(self.size)] for cap in caps]
            for entry, bits in enumerate(self.labels):
                i, j = cut_location(bits, q, side)
                offset = entry * degree
                for a in range(self.size):
                    if left[a][i] == (0, 0):
                        continue
                    cap = caps[a]
                    truncation = (1 << cap) - 1
                    for b in range(self.size):
                        if right[j][b] == (0, 0):
                            continue
                        outer = polynomial_product(left[a][i], right[j][b], cap)
                        outer = tuple((plane ^ (plane << 1)) & truncation for plane in outer)
                        if outer == (0, 0):
                            continue
                        for plane, tables in zip(outer, (low_tables, high_tables)):
                            if not plane:
                                continue
                            if side == 1:
                                width = plane.bit_length()
                                reverse = int(format(plane, '0' + str(width) + 'b')[::-1], 2)
                                for power in range(cap):
                                    at = power - scalar - weights[entry]
                                    section = shifted(reverse, at - width + 1) & mask
                                    tables[a][b][power] ^= section << offset
                            else:
                                for power in range(cap):
                                    at = power - scalar - weights[entry] + 1
                                    section = shifted(plane, -at) & mask
                                    tables[a][b][power] ^= section << offset
            for a, cap in enumerate(caps):
                rows = []
                for power in range(cap):
                    for b in range(self.size):
                        vector = low_tables[a][b][power], high_tables[a][b][power]
                        rhs = coefficient(diagonal[a][b], power)
                        if vector != (0, 0) or rhs:
                            rows.append((power, vector, rhs))
                self.equations.append(tuple(rows))
        self.orders, self.caps = tuple(self.orders), tuple(self.caps)

    def check_endpoint(self, values):
        for bits, pair in zip(self.labels, values):
            x, y, u, v = bits
            value = (pair[0].bit_count() & 1) | ((pair[1].bit_count() & 1) << 1)
            if value != int(x == v and y == u):
                raise ValueError('The full tensor must have the literal crossing endpoint.')

    def weighted_matrix(self, values, side, scalar=None):
        h = self.heights[4*side:4*side+4]
        if scalar is None:
            scalar = self.scalar_shifts[side]
        matrix = [[(0, 0)] * self.size for _ in range(self.size)]
        for bits, pair in zip(self.labels, values):
            i, j = cut_location(bits, self.q, side)
            weight = sum(row[label] for row, label in zip(h, bits))
            result = []
            for plane in pair:
                value = 0
                for power in polynomial_degrees(plane):
                    exponent = scalar + weight + (-power if side == 0 else power)
                    assert exponent >= 0
                    value ^= 1 << exponent
                result.append(value)
            matrix[i][j] = tuple(result)
        return matrix

    def solve(self, targets):
        if len(targets) != 2*self.size or any(t < 0 or t > cap for t, cap in zip(targets, self.caps)):
            raise ValueError('Targets must stay within the completely constructed coefficient window.')
        pivots = {}
        for target, rows in zip(targets, self.equations):
            for power, vector, rhs in rows:
                if power >= target:
                    break
                if not insert(pivots, vector, rhs):
                    return None
        return pivots

    def seed(self, value=(0, 0)):
        mask = (1 << self.degree) - 1
        result = []
        for entry, base in enumerate(self.base):
            first = (value[0] >> (entry*self.degree)) & mask
            second = (value[1] >> (entry*self.degree)) & mask
            result.append((base[0] ^ first ^ (first << 1), base[1] ^ second ^ (second << 1)))
        return tuple(result)

    def profiles(self, values):
        result = []
        for side in range(2):
            matrix = [[(0, 0)] * self.size for _ in range(self.size)]
            for bits, pair in zip(self.labels, values):
                i, j = cut_location(bits, self.q, side)
                matrix[i][j] = tuple(pair)
            h = self.heights[4*side:4*side+4]
            first, second = ((0, 2), (1, 3)) if side == 0 else ((0, 1), (2, 3))
            rows = [h[first[0]][a] + h[first[1]][b] for a, b in product(range(self.q), repeat=2)]
            columns = [h[second[0]][a] + h[second[1]][b] for a, b in product(range(self.q), repeat=2)]
            observed = MatrixOracle(matrix, self.size).profile(rows, columns, infinity=side == 0)['orders']
            assert all(abs(value - round(value)) < 1e-8 for value in observed)
            result.append(tuple(map(round, observed)))
        return tuple(result)

    def check_witness(self, values, targets):
        self.check_endpoint(values)
        for side, (left, right) in enumerate(self.bases):
            precision = max(targets[side*self.size:(side+1)*self.size])
            matrix = self.weighted_matrix(values, side)
            transformed = matrix_product(matrix_product(left, matrix, precision), right, precision)
            for i, row in enumerate(transformed):
                mask = (1 << targets[side*self.size+i]) - 1
                assert all(not ((value[0] | value[1]) & mask) for value in row)
        profiles = self.profiles(values)
        for side in range(2):
            lower = sorted(targets[side*self.size+i] - self.scalar_shifts[side] for i in range(self.size))
            for k in range(1, self.size+1):
                assert sum(profiles[side][:k]) >= sum(lower[:k])
        return profiles


class CompressedF4Fiber:
    def __init__(self, lattice, floors):
        self.lattice, self.floors = lattice, tuple(floors)
        self.base_pivots = lattice.solve(self.floors)
        if self.base_pivots is None:
            raise ValueError('The proposed common row fiber is empty.')
        self.free = tuple(i for i in range(lattice.variables) if i not in self.base_pivots)
        self.positions = {value: i for i, value in enumerate(self.free)}
        self.variables = len(self.free)
        self.equations = []
        for rows in lattice.equations:
            result = []
            for power, vector, rhs in rows:
                reduced, rhs = self.reduce_equation(vector, rhs)
                if reduced != (0, 0) or rhs:
                    result.append((power, reduced, rhs))
            self.equations.append(tuple(result))

    def reduce_equation(self, vector, rhs):
        reduced, rhs = reduce_row(self.base_pivots, vector, rhs)
        result = []
        for plane in reduced:
            value = 0
            while plane:
                bit = plane & -plane
                plane ^= bit
                value ^= 1 << self.positions[bit.bit_length()-1]
            result.append(value)
        return tuple(result), rhs

    def solve(self, targets=None, extra_equations=()):
        targets = self.floors if targets is None else tuple(targets)
        if len(targets) != len(self.floors) or any(t < low or t > cap for t, low, cap in zip(targets, self.floors, self.lattice.caps)):
            raise ValueError('Targets must stay inside the completely reduced fiber window.')
        pivots = {}
        for target, rows in zip(targets, self.equations):
            for power, vector, rhs in rows:
                if power >= target:
                    break
                if not insert(pivots, vector, rhs):
                    return None
        for vector, rhs in extra_equations:
            if not insert(pivots, vector, rhs):
                return None
        return pivots

    def increment(self, pivots, row, target):
        if target >= self.lattice.caps[row]:
            return None
        result = dict(pivots)
        for power, vector, rhs in self.equations[row]:
            if power > target:
                break
            if power == target and not insert(result, vector, rhs):
                return None
        return result

    def sample(self, pivots, rng=None):
        value = (0, 0)
        if rng is not None:
            for index in range(self.variables):
                if index not in pivots:
                    scalar = rng.randrange(4)
                    value = value[0] ^ ((scalar & 1) << index), value[1] ^ ((scalar >> 1) << index)
        value = back_substitute(pivots, value)
        assert all(field_dot(row, value) == rhs for row, rhs in pivots.values())
        original = []
        for plane in value:
            result = 0
            while plane:
                bit = plane & -plane
                plane ^= bit
                result ^= 1 << self.free[bit.bit_length()-1]
            original.append(result)
        return self.lattice.seed(back_substitute(self.base_pivots, tuple(original)))

    def jet_equation(self, entry, point, order, target=0):
        degree = self.lattice.degree
        row = hasse_row(point, order, degree)
        mask = (1 << degree) - 1
        vector = tuple(((plane ^ (plane >> 1)) & mask) << (entry*degree) for plane in row)
        rhs = field_dot(row, self.lattice.base[entry]) ^ target
        return self.reduce_equation(vector, rhs)
