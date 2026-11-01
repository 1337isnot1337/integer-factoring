"""Exact crossing-frame fibers with polynomial-module shift reuse.

Only primitive images are multiplied by the regular frames. Monomial
shifts are then applied exactly at each local place, retaining the
additional infinity precision needed when a shift lowers an order.
The reference and module primitives are binary; coefficients of the
affine fiber are arbitrary F4 elements.
"""

from native_matrix_smith_f4 import MatrixOracle
from f4_recentered_lattices import (
    RecenteredCrossingLattice, add_equation, assignment, regular_frame, unpack,
)
from neutral_reflected_crossing import polynomial_multiply
from source_slab_module_saturation import insert_binary


class FrameOnlyCrossingLattice:
    """The literal reference frames, without the full entry-coordinate system."""
    location = RecenteredCrossingLattice.location
    validate_endpoint = RecenteredCrossingLattice.validate_endpoint
    local_matrix = RecenteredCrossingLattice.local_matrix
    profiles = RecenteredCrossingLattice.profiles
    check_witness = RecenteredCrossingLattice.check_witness

    def __init__(self, values, q, degree, heights, raise_cap=5):
        if len(values) != q ** 4 or len(heights) != 8:
            raise ValueError('Supply the complete tensor and eight physical height lists.')
        if any(len(row) != q or any(not isinstance(value, int) for value in row) for row in heights):
            raise ValueError('Every physical height must be an explicit integer.')
        if any(len(pair) != 2 or any(not isinstance(value, int) or value < 0 or value.bit_length() > degree + 1
                                    for value in pair) for pair in values):
            raise ValueError('Every coefficient must fit the complete degree bound.')
        self.q, self.size, self.degree = q, q * q, degree
        self.base = tuple(map(tuple, values))
        self.heights = tuple(map(tuple, heights))
        self.labels = tuple(tuple((index // q ** leg) % q for leg in range(4)) for index in range(q ** 4))
        self.validate_endpoint(values)
        self.frames, self.diagonals = [], []
        self.orders, self.caps = [], []
        self.scalar_shifts, self.precisions = [], []
        for side in range(2):
            heights_side = self.heights[4 * side:4 * side + 4]
            scalar = -sum(min(row) for row in heights_side) + (degree if side == 0 else 0)
            self.scalar_shifts.append(scalar)
            matrix = self.local_matrix(values, side, scalar)
            native = MatrixOracle([[unpack(value) for value in row] for row in matrix], self.size).profile()['orders']
            assert all(abs(value - round(value)) < 1e-8 for value in native)
            precision = max(map(round, native)) + raise_cap
            orders, left, right, diagonal = regular_frame(matrix, precision)
            assert tuple(map(round, native)) == orders
            self.frames.append((left, right))
            self.diagonals.append(diagonal)
            self.orders.extend(orders)
            self.caps.extend(value + raise_cap for value in orders)
            self.precisions.append(precision)
        self.orders, self.caps = tuple(self.orders), tuple(self.caps)


def binary_matrix(matrix):
    result = []
    for row in matrix:
        values = []
        for value in row:
            low, high = unpack(value)
            if high:
                raise ValueError('This optimized constructor requires binary reference frames.')
            values.append(low)
        result.append(values)
    return result


def binary_matrix_product(left, right, precision):
    mask = (1 << precision) - 1
    result = [[0] * len(right[0]) for _ in left]
    sparse_right = [[(j, value) for j, value in enumerate(row) if value] for row in right]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if value:
                for j, other in sparse_right[k]:
                    result[i][j] ^= polynomial_multiply(value, other)
        result[i] = [value & mask for value in result[i]]
    return result


class SourceModuleCrossingFiber:
    """Exact full F4 span of the displayed binary slab/placement/shift list."""
    def __init__(self, values, source, heights, raise_cap=5):
        q = source['dimensions'][0]
        if source['dimensions'] != [q] * 4:
            raise ValueError('Use four equal physical alphabets.')
        degree = source['degree']
        self.full_lattice = FrameOnlyCrossingLattice(values, q, degree, heights, raise_cap)
        self.q, self.size, self.degree = q, q * q, degree
        self.base = self.full_lattice.base
        self.orders, self.caps = self.full_lattice.orders, self.full_lattice.caps
        self.floors = (0,) * len(self.orders)
        self.primitives = tuple(tuple(item['normalized_slab']) for item in source['primitives'])
        self.descriptors = tuple(map(tuple, source['direction_descriptors']))
        self.variables = len(self.descriptors)
        if len(set(self.descriptors)) != self.variables:
            raise ValueError('Direction descriptors must be distinct.')
        self.source_direction_count = self.variables
        self.selected_direction_indices = list(range(self.variables))
        self.image_pivots = [{} for _ in range(q)]
        self.groups = {}
        for column, (primitive, first, shift) in enumerate(self.descriptors):
            if not 0 <= primitive < len(self.primitives) or not 0 <= first < q or shift < 0:
                raise ValueError('Malformed polynomial-module descriptor.')
            row = self.primitives[primitive]
            if len(row) != q ** 3 or any(value < 0 or (value << shift).bit_length() > degree + 1
                                       or value.bit_count() % 2 for value in row):
                raise ValueError('Primitive must be endpoint zero and fit the complete degree.')
            image = sum((value << shift) << ((degree + 1) * j) for j, value in enumerate(row))
            if not insert_binary(self.image_pivots[first], image):
                raise ValueError('Supply an independent polynomial-module direction list.')
            self.groups.setdefault((primitive, first), []).append((column, shift))
        self.equations = []
        self.equation_positions = []
        self.primitive_frame_images = 0
        for side in range(2):
            caps = self.caps[side * self.size:(side + 1) * self.size]
            tables = [[[0] * cap for _ in range(self.size)] for cap in caps]
            left, right = map(binary_matrix, self.full_lattice.frames[side])
            for (primitive, first), columns in self.groups.items():
                literal = [(0, 0)] * q ** 4
                for j, value in enumerate(self.primitives[primitive]):
                    literal[first + q * j] = value, 0
                matrix = binary_matrix(self.full_lattice.local_matrix(literal, side))
                # Infinity shifts lower orders, so compute the extra coefficients.
                maximum_shift = max(shift for _, shift in columns)
                precision = max(caps) + (maximum_shift if side == 0 else 0)
                transformed = binary_matrix_product(binary_matrix_product(left, matrix, precision), right, precision)
                self.primitive_frame_images += 1
                for column, shift in columns:
                    coefficient = 1 << column
                    for a, row in enumerate(transformed):
                        mask = (1 << caps[a]) - 1
                        for b, value in enumerate(row):
                            value = ((value >> shift) if side == 0 else (value << shift)) & mask
                            while value:
                                bit = value & -value
                                value ^= bit
                                tables[a][b][bit.bit_length() - 1] ^= coefficient
            diagonal = self.full_lattice.diagonals[side]
            for a, cap in enumerate(caps):
                equations = []
                positions = []
                for power in range(cap):
                    for b in range(self.size):
                        low = tables[a][b][power]
                        rhs = (diagonal[a][b] >> (2 * power)) & 3
                        if low or rhs:
                            equations.append((power, low, 0, rhs))
                            positions.append((power, b))
                self.equations.append(equations)
                self.equation_positions.append(positions)

    def solve(self, targets):
        if len(targets) != len(self.orders) or any(target < 0 or target > cap for target, cap in zip(targets, self.caps)):
            raise ValueError('Targets must fit the complete row coefficient windows.')
        pivots = {}
        for target, row in zip(targets, self.equations):
            for power, low, high, rhs in row:
                if power >= target:
                    break
                if not add_equation(pivots, low, high, rhs):
                    return None
        return pivots

    def seed(self, low, high):
        if low < 0 or high < 0 or max(low, high).bit_length() > self.variables:
            raise ValueError('Parameter vector exceeds the prescribed module.')
        result = [list(pair) for pair in self.base]
        active = low | high
        while active:
            bit = active & -active
            active ^= bit
            column = bit.bit_length() - 1
            primitive, first, shift = self.descriptors[column]
            for j, value in enumerate(self.primitives[primitive]):
                if low & bit:
                    result[first + self.q * j][0] ^= value << shift
                if high & bit:
                    result[first + self.q * j][1] ^= value << shift
        return tuple(map(tuple, result))

    def sample(self, pivots, rng):
        free = ((1 << self.variables) - 1) ^ sum(1 << pivot for pivot in pivots)
        low, high = assignment(pivots, rng.getrandbits(self.variables) & free,
                              rng.getrandbits(self.variables) & free)
        return self.seed(low, high)

    def contains(self, values):
        if len(values) != self.q ** 4:
            return False
        for first in range(self.q):
            for plane in range(2):
                image = 0
                for j in range(self.q ** 3):
                    entry = first + self.q * j
                    if len(values[entry]) != 2 or any(value < 0 or value.bit_length() > self.degree + 1 for value in values[entry]):
                        return False
                    image |= (values[entry][plane] ^ self.base[entry][plane]) << ((self.degree + 1) * j)
                while image:
                    pivot = image.bit_length() - 1
                    if pivot not in self.image_pivots[first]:
                        return False
                    image ^= self.image_pivots[first][pivot]
        return True

    def profiles(self, values):
        return self.full_lattice.profiles(values)

    def check_witness(self, values, targets):
        assert self.contains(values)
        return self.full_lattice.check_witness(values, targets)
