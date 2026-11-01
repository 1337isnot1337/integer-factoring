#!/usr/bin/env python3
"""Complete affine coefficient equations with genuinely F4-valued frames.

Each F4 equation is represented by two coupled binary equations. The right
regular frame can be omitted from row divisibility, reducing construction
work. Historical binary-frame sources and receipts remain unchanged.
"""
from __future__ import annotations

from itertools import product

from full_alphabet_two_cut_lattices import (FullAlphabetTwoCutLattice,
    CompressedFullAlphabetFiber, cut_location, shifted)
from native_matrix_smith_f4 import MatrixOracle
from windowed_f4_smith import matrix_product, pack, planes, regular_smith_frames


class F4AffineLattice(FullAlphabetTwoCutLattice):
    def __init__(self, values, q, degree, heights, raise_cap=8):
        if q < 2 or degree < 1 or raise_cap < 1 or len(values) != q**4 or len(heights) != 8:
            raise ValueError('Supply a complete positive-degree crossing and eight height lists.')
        if any(len(row) != q or any(type(v) is not int for v in row) for row in heights):
            raise ValueError('Use one integer height per physical letter.')
        if any(len(pair) != 2 or any(type(v) is not int or v < 0 or v.bit_length() > degree+1
                                    for v in pair) for pair in values):
            raise ValueError('Supply F4 bit-plane polynomial coefficients within the degree bound.')
        self.q, self.size, self.degree = q, q*q, degree
        self.coefficient_variables = q**4*degree
        self.variables = 2*self.coefficient_variables
        self.base = tuple(pack(pair) for pair in values)
        self.heights = tuple(map(tuple, heights))
        self.labels = tuple(tuple((i//q**j) % q for j in range(4)) for i in range(q**4))
        for pair, (x, y, u, v) in zip(values, self.labels):
            endpoint = (pair[0].bit_count() % 2) | ((pair[1].bit_count() % 2) << 1)
            if endpoint != int(x == v and y == u):
                raise ValueError('The source must have the literal crossing endpoint.')
        self.orders, self.caps, self.equations = [], [], []
        self.scalar_shifts, self.left_frames, self.frame_windows = [], [], []
        self.profile_retries = []
        variable_mask, count = (1 << degree)-1, self.coefficient_variables
        for side in range(2):
            h = self.heights[4*side:4*side+4]
            scalar = -sum(min(row) for row in h)+(degree if side == 0 else 0)
            self.scalar_shifts.append(scalar)
            weights = [sum(row[label] for row, label in zip(h, bits)) for bits in self.labels]
            source = [[0]*self.size for _ in range(self.size)]
            for entry, packed in enumerate(self.base):
                i, j = cut_location(self.labels[entry], q, side)
                power = 0
                while packed:
                    coefficient = packed & 3
                    if coefficient:
                        exponent = scalar+weights[entry]+(-power if side == 0 else power)
                        assert exponent >= 0
                        source[i][j] ^= coefficient << (2*exponent)
                    packed >>= 2
                    power += 1
            orders, left, _, _, metadata = regular_smith_frames(source, raise_cap)
            self.orders.extend(orders)
            caps = [v+raise_cap for v in orders]
            self.caps.extend(caps)
            self.left_frames.append(left)
            self.frame_windows.append(metadata)
            baseline = matrix_product(left, source, metadata['output_precision'])
            tables = [[[(0, 0)]*cap for _ in range(self.size)] for cap in caps]
            for entry, bits in enumerate(self.labels):
                i, j = cut_location(bits, q, side)
                offset = entry*degree
                for a in range(self.size):
                    if not left[a][i]:
                        continue
                    outer = left[a][i] ^ (left[a][i] << 2)
                    low, high = planes(outer)
                    low_width, high_width = low.bit_length(), high.bit_length()
                    low_reverse = int(format(low, 'b')[::-1], 2) if low else 0
                    high_reverse = int(format(high, 'b')[::-1], 2) if high else 0
                    for power in range(caps[a]):
                        if side == 0:
                            at = power-scalar-weights[entry]+1
                            lo = shifted(low, -at) & variable_mask
                            hi = shifted(high, -at) & variable_mask
                        else:
                            at = power-scalar-weights[entry]
                            lo = shifted(low_reverse, at-low_width+1) & variable_mask
                            hi = shifted(high_reverse, at-high_width+1) & variable_mask
                        old_low, old_high = tables[a][j][power]
                        tables[a][j][power] = old_low ^ (lo << offset), old_high ^ (hi << offset)
            for a, cap in enumerate(caps):
                equations = []
                for power in range(cap):
                    for b in range(self.size):
                        low, high = tables[a][b][power]
                        rhs = (baseline[a][b] >> (2*power)) & 3
                        binary_low = low | (high << count)
                        binary_high = high | ((low ^ high) << count)
                        if binary_low or rhs & 1:
                            equations.append((power, binary_low, rhs & 1))
                        if binary_high or rhs >> 1:
                            equations.append((power, binary_high, rhs >> 1))
                self.equations.append(equations)
        self.orders, self.caps = tuple(self.orders), tuple(self.caps)

    def seed(self, assignment):
        count, mask = self.coefficient_variables, (1 << self.degree)-1
        if type(assignment) is not int or assignment < 0 or assignment >> self.variables:
            raise ValueError('The assignment must fit both binary coefficient planes.')
        values = []
        for entry, base in enumerate(self.base):
            offset = entry*self.degree
            low = (assignment >> offset) & mask
            high = (assignment >> (count+offset)) & mask
            q = pack((low, high))
            values.append(planes(base ^ q ^ (q << 2)))
        return tuple(values)

    def sample(self, pivots, rng):
        forbidden = sum(1 << p for p in pivots)
        free = rng.getrandbits(self.variables) & ~forbidden
        return self.seed(self.assignment(pivots, free))

    def profiles(self, values):
        profiles = []
        for side in range(2):
            matrix = [[(0, 0)]*self.size for _ in range(self.size)]
            for entry, pair in enumerate(values):
                i, j = cut_location(self.labels[entry], self.q, side)
                matrix[i][j] = tuple(pair)
            h = self.heights[4*side:4*side+4]
            rw = [h[0][a]+h[2 if side == 0 else 1][b]
                  for a, b in product(range(self.q), repeat=2)]
            cw = [h[1 if side == 0 else 2][a]+h[3][b]
                  for a, b in product(range(self.q), repeat=2)]
            for precision in (512, 1024, 2048, 4096, 8192):
                try:
                    orders = MatrixOracle(matrix, self.size, precision=precision).profile(
                        rw, cw, 1, infinity=side == 0)['orders']
                    break
                except ArithmeticError as error:
                    self.profile_retries.append({'side': side, 'precision': precision, 'error': str(error)})
            else:
                raise ArithmeticError('No complete native F4 profile resolved within the supported windows.')
            assert all(value == round(value) for value in orders)
            profiles.append(tuple(round(value) for value in orders))
        return tuple(profiles)

    def check_witness(self, values, targets):
        if len(values) != self.q**4 or len(targets) != 2*self.size:
            raise ValueError('Supply the complete tensor and both row-target vectors.')
        if any(t < 0 or t > c for t, c in zip(targets, self.caps)):
            raise ValueError('A target lies outside the complete coefficient window.')
        for pair, (x, y, u, v) in zip(values, self.labels):
            assert ((pair[0].bit_count() % 2) | ((pair[1].bit_count() % 2) << 1)) == int(x == v and y == u)
        for side, left in enumerate(self.left_frames):
            source = [[0]*self.size for _ in range(self.size)]
            h, scalar = self.heights[4*side:4*side+4], self.scalar_shifts[side]
            for entry, pair in enumerate(values):
                i, j = cut_location(self.labels[entry], self.q, side)
                weight = sum(row[label] for row, label in zip(h, self.labels[entry]))
                packed, power = pack(pair), 0
                while packed:
                    coefficient = packed & 3
                    if coefficient:
                        exponent = scalar+weight+(-power if side == 0 else power)
                        assert exponent >= 0
                        source[i][j] ^= coefficient << (2*exponent)
                    packed >>= 2
                    power += 1
            transformed = matrix_product(left, source, max(targets[side*self.size:(side+1)*self.size])+1)
            for i, row in enumerate(transformed):
                mask = (1 << (2*targets[side*self.size+i]))-1
                assert all(not (v & mask) for v in row)
        profiles = self.profiles(values)
        for side in range(2):
            bounds = sorted(targets[side*self.size+i]-self.scalar_shifts[side] for i in range(self.size))
            assert all(sum(profiles[side][:k]) >= sum(bounds[:k]) for k in range(1, self.size+1))
        return profiles


class CoupledF4Fiber(CompressedFullAlphabetFiber):
    def sample(self, pivots, rng):
        forbidden = sum(1 << p for p in pivots)
        value = rng.getrandbits(self.variables) & ~forbidden
        for pivot, (row, rhs) in sorted(pivots.items()):
            value |= (((row & value).bit_count() % 2) ^ rhs) << pivot
        original = 0
        while value:
            bit = value & -value
            value ^= bit
            original |= 1 << self.free[bit.bit_length()-1]
        return self.lattice.seed(self.lattice.assignment(self.base_pivots, original))
