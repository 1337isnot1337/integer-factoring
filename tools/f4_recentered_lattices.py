#!/usr/bin/env python3
"""Exact complete-entry crossing fibers with genuinely F4 reference frames.

Unlike the earlier binary-frame implementation, a retained nonbinary tensor
can be the reference for the next solve.  Every coefficient is still free,
and the strict crossing endpoint and complete degree bound remain explicit.
Only coefficient jets below a proved finite row cap enter the affine system.
"""
from __future__ import annotations

from itertools import product
import random

from native_matrix_smith_f4 import MatrixOracle
from weighted_smith_oracle_f4 import polynomial_multiply, scalar_multiply


INVERSE = (0, 1, 3, 2)


def pack(pair):
    """Change two binary coefficient planes to two bits per F4 coefficient."""
    low, high = pair
    result = 0
    while low:
        bit = low & -low
        low ^= bit
        result ^= 1 << (2 * (bit.bit_length() - 1))
    while high:
        bit = high & -high
        high ^= bit
        result ^= 2 << (2 * (bit.bit_length() - 1))
    return result


def unpack(value):
    low = high = 0
    power = 0
    while value:
        low |= (value & 1) << power
        high |= ((value >> 1) & 1) << power
        value >>= 2
        power += 1
    return low, high


def degree(value):
    return (value.bit_length() - 1) // 2 if value else -1


def order(value):
    if not value:
        raise ValueError('The zero polynomial has no finite order.')
    return ((value & -value).bit_length() - 1) // 2


def multiply_matrix(left, right, precision=None):
    mask = None if precision is None else (1 << (2 * precision)) - 1
    result = [[0] * len(right[0]) for _ in left]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if not value:
                continue
            for j, other in enumerate(right[k]):
                if other:
                    result[i][j] ^= polynomial_multiply(value, other)
        if mask is not None:
            result[i] = [value & mask for value in result[i]]
    return result


def f4_rank(rows):
    """Rank of literal scalar F4 rows, used only for independent frame checks."""
    work = [list(row) for row in rows]
    pivot_row = 0
    for column in range(len(work[0]) if work else 0):
        candidates = [i for i in range(pivot_row, len(work)) if work[i][column]]
        if not candidates:
            continue
        i = candidates[0]
        work[pivot_row], work[i] = work[i], work[pivot_row]
        inverse = INVERSE[work[pivot_row][column]]
        work[pivot_row] = [scalar_multiply(value, inverse) for value in work[pivot_row]]
        for i in range(pivot_row + 1, len(work)):
            factor = work[i][column]
            if factor:
                work[i] = [a ^ scalar_multiply(b, factor)
                           for a, b in zip(work[i], work[pivot_row])]
        pivot_row += 1
        if pivot_row == len(work):
            break
    return pivot_row


def regular_frame(source, precision):
    """Actual polynomial regular frames, with an exact finite-jet identity.

    Pivots are not inverted as rational functions.  A row/column elimination
    multiplies by its pivot unit, whose constant coefficient is nonzero.
    Thus each operation is invertible in F4[[z]], even after its polynomial
    representative is truncated modulo z**precision.
    """
    n = len(source)
    if not n or any(len(row) != n for row in source) or precision < 1:
        raise ValueError('Use a nonempty square matrix and positive precision.')
    mask = (1 << (2 * precision)) - 1
    work = [[value & mask for value in row] for row in source]
    left = [[int(i == j) for j in range(n)] for i in range(n)]
    right = [[int(i == j) for j in range(n)] for i in range(n)]
    orders = []
    for k in range(n):
        candidates = [(order(work[i][j]), i, j)
                      for i in range(k, n) for j in range(k, n) if work[i][j]]
        if not candidates:
            raise ValueError('The finite precision does not expose every promised pivot.')
        valuation, i, j = min(candidates)
        work[k], work[i] = work[i], work[k]
        left[k], left[i] = left[i], left[k]
        for values in (work, right):
            for row in values:
                row[k], row[j] = row[j], row[k]
        unit = work[k][k] >> (2 * valuation)
        assert unit & 3
        for i in range(k + 1, n):
            if not work[i][k]:
                continue
            assert order(work[i][k]) >= valuation
            quotient = work[i][k] >> (2 * valuation)
            for values in (work, left):
                values[i] = [(polynomial_multiply(unit, a) ^
                              polynomial_multiply(quotient, b)) & mask
                             for a, b in zip(values[i], values[k])]
            assert work[i][k] == 0
        for j in range(k + 1, n):
            if not work[k][j]:
                continue
            assert order(work[k][j]) >= valuation
            quotient = work[k][j] >> (2 * valuation)
            for values in (work, right):
                for row in values:
                    row[j] = (polynomial_multiply(unit, row[j]) ^
                              polynomial_multiply(quotient, row[k])) & mask
            assert work[k][j] == 0
        orders.append(valuation)
    assert orders == sorted(orders)
    assert all(not work[i][j] for i in range(n) for j in range(n) if i != j)
    assert multiply_matrix(multiply_matrix(left, source, precision), right, precision) == work
    for frame in (left, right):
        assert f4_rank([[value & 3 for value in row] for row in frame]) == n
    return tuple(orders), left, right, work


def scale_row(low, high, coefficient):
    if coefficient == 0:
        return 0, 0
    if coefficient == 1:
        return low, high
    if coefficient == 2:
        return high, low ^ high
    if coefficient == 3:
        return low ^ high, low
    raise ValueError('F4 coefficients are integers from zero through three.')


def add_equation(pivots, low, high, rhs):
    """Insert one F4 affine equation; normalized pivot coefficients are one."""
    while low | high:
        pivot = (low | high).bit_length() - 1
        coefficient = ((low >> pivot) & 1) | (((high >> pivot) & 1) << 1)
        if pivot not in pivots:
            inverse = INVERSE[coefficient]
            low, high = scale_row(low, high, inverse)
            pivots[pivot] = low, high, scalar_multiply(rhs, inverse)
            return True
        a, b, c = pivots[pivot]
        a, b = scale_row(a, b, coefficient)
        low ^= a
        high ^= b
        rhs ^= scalar_multiply(c, coefficient)
    return not rhs


def dot(low, high, a, b):
    return (((low & a).bit_count() ^ (high & b).bit_count()) & 1) | (
        (((low & b).bit_count() ^ (high & a).bit_count() ^ (high & b).bit_count()) & 1) << 1)


def assignment(pivots, low=0, high=0):
    for pivot, (a, b, rhs) in sorted(pivots.items()):
        if ((low | high) >> pivot) & 1:
            raise ValueError('A free assignment must leave pivot coordinates zero.')
        value = rhs ^ dot(a, b, low, high)
        low |= (value & 1) << pivot
        high |= ((value >> 1) & 1) << pivot
    return low, high


def shifted(value, shift):
    return value << shift if shift >= 0 else value >> (-shift)


class AffineCoordinateMap:
    """Literal free-coordinate substitution, with eight-column XOR tables.

    Each original variable is expressed in the free variables by triangular
    back substitution.  Tables only cache sums of those expressions.  They
    change neither the coefficient field nor the affine equation space.
    """
    def __init__(self, variables, pivots):
        self.free = tuple(i for i in range(variables) if i not in pivots)
        self.variables = len(self.free)
        self.width = self.variables+1  # Last coordinate is the affine constant.
        self.mask = (1 << self.width)-1
        self.expressions = []
        self.tables = []
        self.input_bytes = (variables+7)//8
        positions = {value: i for i, value in enumerate(self.free)}
        for variable in range(variables):
            if variable in positions:
                value = 1 << positions[variable]
            else:
                low, high, rhs = pivots[variable]
                assert (low >> variable) == 1 and (high >> variable) == 0
                value = self.combine(low ^ (1 << variable), high)
                value ^= ((rhs & 1) << self.variables) | ((rhs >> 1) << (self.width+self.variables))
            self.expressions.append(value)
            if len(self.expressions) % 8 == 0:
                self.append_table(self.expressions[-8:])
        if len(self.expressions) % 8:
            self.append_table(self.expressions[8*len(self.tables):])

    def append_table(self, expressions):
        table = [0]*(1 << len(expressions))
        for index in range(1, len(table)):
            bit = index & -index
            table[index] = table[index ^ bit] ^ expressions[bit.bit_length()-1]
        self.tables.append(table)

    def combine(self, low, high):
        """Substitute variables in a linear form, retaining its affine constant."""
        size = (max(low.bit_length(), high.bit_length())+7)//8
        lows, highs = low.to_bytes(size, 'little'), high.to_bytes(size, 'little')
        a = b = 0
        complete = min(size, len(self.tables))
        for i in range(complete):
            table = self.tables[i]
            a ^= table[lows[i]]
            b ^= table[highs[i]]
        if size > complete:
            # During triangular construction the last block is incomplete.
            assert size == complete+1
            first = 8*complete
            lo, hi = lows[-1], highs[-1]
            while lo:
                bit = lo & -lo
                lo ^= bit
                a ^= self.expressions[first+bit.bit_length()-1]
            while hi:
                bit = hi & -hi
                hi ^= bit
                b ^= self.expressions[first+bit.bit_length()-1]
        low_b, high_b = b & self.mask, b >> self.width
        return a ^ high_b ^ ((low_b ^ high_b) << self.width)

    def reduce(self, low, high, rhs):
        value = self.combine(low, high)
        a, b = value & self.mask, value >> self.width
        constant = (a >> self.variables) | ((b >> self.variables) << 1)
        variable_mask = (1 << self.variables)-1
        return a & variable_mask, b & variable_mask, rhs ^ constant


class RecenteredCrossingLattice:
    """All degree-bounded F4 coefficients around an arbitrary F4 crossing."""
    def __init__(self, values, q, degree_bound, heights, raise_cap=4):
        if q < 2 or len(values) != q**4 or len(heights) != 8:
            raise ValueError('Supply every q**4 entry and eight physical height lists.')
        if degree_bound < 1 or raise_cap < 1:
            raise ValueError('Degree and row allowance must be positive.')
        if any(len(row) != q or any(not isinstance(v, int) for v in row) for row in heights):
            raise ValueError('Every physical height must be an explicit integer.')
        if any(len(pair) != 2 or any(not isinstance(v, int) or v < 0 or
                                    v.bit_length() > degree_bound + 1 for v in pair)
               for pair in values):
            raise ValueError('Every F4 coefficient table must fit the complete degree bound.')
        self.q = q
        self.size = q*q
        self.degree = degree_bound
        self.variables = q**4 * degree_bound
        self.base = tuple(tuple(pair) for pair in values)
        self.heights = tuple(map(tuple, heights))
        self.labels = tuple(tuple((index // q**i) % q for i in range(4))
                            for index in range(q**4))
        self.validate_endpoint(values)
        self.frames = []
        self.orders = []
        self.caps = []
        self.equations = []
        self.scalar_shifts = []
        self.precisions = []
        variable_mask = (1 << self.degree) - 1
        for side in range(2):
            h = self.heights[4*side:4*side+4]
            scalar = -sum(min(row) for row in h) + (self.degree if side == 0 else 0)
            self.scalar_shifts.append(scalar)
            weights = [sum(row[label] for row, label in zip(h, bits)) for bits in self.labels]
            source = self.local_matrix(values, side, scalar)
            oracle = MatrixOracle([[unpack(value) for value in row] for row in source], self.size)
            native = oracle.profile()['orders']
            if any(abs(value - round(value)) > 1e-8 for value in native):
                raise ArithmeticError('Integer polynomial orders must be integers.')
            expected = tuple(map(round, native))
            precision = max(expected) + raise_cap
            orders, left, right, diagonal = regular_frame(source, precision)
            assert orders == expected
            self.frames.append((left, right))
            self.precisions.append(precision)
            self.orders.extend(orders)
            caps = [value + raise_cap for value in orders]
            self.caps.extend(caps)
            tables = [[[[0, 0] for _ in range(cap)] for _ in range(self.size)]
                      for cap in caps]
            for entry, bits in enumerate(self.labels):
                i, j = self.location(bits, side)
                offset = entry * self.degree
                for a in range(self.size):
                    if not left[a][i]:
                        continue
                    cap = caps[a]
                    truncation = (1 << (2*cap)) - 1
                    for b in range(self.size):
                        if not right[j][b]:
                            continue
                        outer = polynomial_multiply(left[a][i], right[j][b])
                        outer = (outer ^ (outer << 2)) & truncation
                        if not outer:
                            continue
                        for plane, coefficients in enumerate(unpack(outer)):
                            if not coefficients:
                                continue
                            if side == 1:
                                width = coefficients.bit_length()
                                reverse = int(format(coefficients, '0' + str(width) + 'b')[::-1], 2)
                                for power in range(cap):
                                    at = power - scalar - weights[entry]
                                    section = shifted(reverse, at-width+1) & variable_mask
                                    tables[a][b][power][plane] ^= section << offset
                            else:
                                for power in range(cap):
                                    at = power - scalar - weights[entry] + 1
                                    section = shifted(coefficients, -at) & variable_mask
                                    tables[a][b][power][plane] ^= section << offset
            for a, cap in enumerate(caps):
                rows = []
                for power in range(cap):
                    for b in range(self.size):
                        low, high = tables[a][b][power]
                        rhs = (diagonal[a][b] >> (2*power)) & 3
                        if low or high or rhs:
                            rows.append((power, low, high, rhs))
                self.equations.append(rows)
        self.orders = tuple(self.orders)
        self.caps = tuple(self.caps)

    def location(self, bits, side):
        x, y, u, v = bits
        return (self.q*x+u, self.q*y+v) if side == 0 else (self.q*x+y, self.q*u+v)

    def validate_endpoint(self, values):
        for bits, pair in zip(self.labels, values):
            x, y, u, v = bits
            endpoint = (pair[0].bit_count() & 1) | ((pair[1].bit_count() & 1) << 1)
            if endpoint != int(x == v and y == u):
                raise ValueError('The reference must have the literal strict crossing endpoint.')

    def local_matrix(self, values, side, scalar=None):
        h = self.heights[4*side:4*side+4]
        if scalar is None:
            scalar = self.scalar_shifts[side]
        result = [[0]*self.size for _ in range(self.size)]
        for bits, pair in zip(self.labels, values):
            i, j = self.location(bits, side)
            weight = sum(row[label] for row, label in zip(h, bits))
            for plane, polynomial in enumerate(pair):
                while polynomial:
                    bit = polynomial & -polynomial
                    polynomial ^= bit
                    power = bit.bit_length()-1
                    exponent = scalar + weight + (-power if side == 0 else power)
                    if exponent < 0:
                        raise ValueError('The explicit scalar clearing left a negative exponent.')
                    result[i][j] ^= (1 << plane) << (2*exponent)
        return result

    def solve(self, targets):
        if len(targets) != 2*self.size or any(t < 0 or t > cap for t, cap in zip(targets, self.caps)):
            raise ValueError('Targets must lie in the completely constructed coefficient window.')
        pivots = {}
        for target, equations in zip(targets, self.equations):
            for power, low, high, rhs in equations:
                if power >= target:
                    break
                if not add_equation(pivots, low, high, rhs):
                    return None
        return pivots

    def seed(self, low, high):
        mask = (1 << self.degree)-1
        values = []
        for entry, (a, b) in enumerate(self.base):
            x = (low >> (entry*self.degree)) & mask
            y = (high >> (entry*self.degree)) & mask
            values.append((a ^ x ^ (x << 1), b ^ y ^ (y << 1)))
        return tuple(values)

    def sample(self, pivots, rng):
        free_mask = ((1 << self.variables)-1) ^ sum(1 << pivot for pivot in pivots)
        low = rng.getrandbits(self.variables) & free_mask
        high = rng.getrandbits(self.variables) & free_mask
        return self.seed(*assignment(pivots, low, high))

    def profiles(self, values):
        output = []
        for side in range(2):
            matrix = self.local_matrix(values, side)
            raw = MatrixOracle([[unpack(value) for value in row] for row in matrix], self.size).profile()['orders']
            assert all(abs(value-round(value)) < 1e-8 for value in raw)
            output.append(tuple(round(value)-self.scalar_shifts[side] for value in raw))
        return tuple(output)

    def check_witness(self, values, targets):
        self.validate_endpoint(values)
        if any(max(pair).bit_length() > self.degree+1 for pair in values):
            raise AssertionError('A witness exceeded the source degree bound.')
        for side, (left, right) in enumerate(self.frames):
            matrix = self.local_matrix(values, side)
            precision = self.precisions[side]
            transformed = multiply_matrix(multiply_matrix(left, matrix, precision), right, precision)
            for i, row in enumerate(transformed):
                mask = (1 << (2*targets[side*self.size+i]))-1
                assert all(not (value & mask) for value in row)
        profiles = self.profiles(values)
        for side in range(2):
            bounds = sorted(targets[side*self.size+i]-self.scalar_shifts[side]
                            for i in range(self.size))
            for k in range(1, self.size+1):
                assert sum(profiles[side][:k]) >= sum(bounds[:k])
        return profiles


class RecenteredFiber:
    """Exact F4 coordinates on a shared affine floor fiber."""
    def __init__(self, lattice, floors, compression='auto'):
        if compression not in ('auto', 'direct', 'tables'):
            raise ValueError('Compression must be auto, direct, or tables.')
        self.lattice = lattice
        self.floors = tuple(floors)
        self.base_pivots = lattice.solve(floors)
        if self.base_pivots is None:
            raise ValueError('The common floor fiber is empty.')
        self.free = tuple(i for i in range(lattice.variables) if i not in self.base_pivots)
        self.variables = len(self.free)
        positions = {value: i for i, value in enumerate(self.free)}
        table_bytes = ((lattice.variables+7)//8)*256*(2*(self.variables+1)//8+48)
        use_tables = compression == 'tables' or (compression == 'auto' and lattice.variables >= 256 and table_bytes < 256_000_000)
        coordinate_map = AffineCoordinateMap(lattice.variables, self.base_pivots) if use_tables else None
        self.compression = 'tables' if use_tables else 'direct'
        self.equations = []
        for row in lattice.equations:
            reduced_rows = []
            for exponent, low, high, rhs in row:
                if coordinate_map is not None:
                    a, b, rhs = coordinate_map.reduce(low, high, rhs)
                    if a or b or rhs:
                        reduced_rows.append((exponent, a, b, rhs))
                    continue
                a = b = 0
                while low | high:
                    pivot = (low | high).bit_length()-1
                    coefficient = ((low >> pivot) & 1) | (((high >> pivot) & 1) << 1)
                    if pivot in self.base_pivots:
                        old_low, old_high, constant = self.base_pivots[pivot]
                        old_low, old_high = scale_row(old_low, old_high, coefficient)
                        low ^= old_low
                        high ^= old_high
                        rhs ^= scalar_multiply(constant, coefficient)
                    else:
                        a |= (coefficient & 1) << positions[pivot]
                        b |= ((coefficient >> 1) & 1) << positions[pivot]
                        low &= ~(1 << pivot)
                        high &= ~(1 << pivot)
                if a or b or rhs:
                    reduced_rows.append((exponent, a, b, rhs))
            self.equations.append(reduced_rows)

    def solve(self, targets):
        if len(targets) != len(self.floors) or any(t < f or t > c for t, f, c in
                                                 zip(targets, self.floors, self.lattice.caps)):
            raise ValueError('Targets must remain in the complete floor-fiber window.')
        pivots = {}
        for target, equations in zip(targets, self.equations):
            for power, low, high, rhs in equations:
                if power >= target:
                    break
                if not add_equation(pivots, low, high, rhs):
                    return None
        return pivots

    def sample(self, pivots, rng):
        free_mask = ((1 << self.variables)-1) ^ sum(1 << pivot for pivot in pivots)
        low, high = assignment(pivots, rng.getrandbits(self.variables) & free_mask,
                               rng.getrandbits(self.variables) & free_mask)
        original_low = original_high = 0
        for i, source in enumerate(self.free):
            original_low |= ((low >> i) & 1) << source
            original_high |= ((high >> i) & 1) << source
        return self.lattice.seed(*assignment(self.base_pivots, original_low, original_high))
