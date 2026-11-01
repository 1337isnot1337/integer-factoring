#!/usr/bin/env python3
"""Original top-table queries with a packed boundary and streamed settings.

No gate/edge metadata table, recovered-setting cache, or interpolation
matrix is retained.  A reference network is used only by the companion
checker, never by these query procedures.
"""

from math import gcd


class InvalidBoundary(Exception):
    pass


class BitTape:
    """Read-only little-endian input bytes, with constant-size bit reads.

    The caller supplies this input representation. In particular the
    sublinear-workspace query does not copy an n-bit integer by shifting
    it. from_int is only a convenient preparation method for tests.
    """

    def __init__(self, data):
        if not isinstance(data, bytes):
            raise TypeError("input tape must be immutable bytes")
        self.data = data
        last = len(data) - 1
        while last >= 0 and data[last] == 0:
            last -= 1
        self.length = 0 if last < 0 else 8 * last + data[last].bit_length()

    @classmethod
    def from_int(cls, value):
        if value < 0:
            raise ValueError("negative input integer")
        return cls(value.to_bytes((value.bit_length() + 7) // 8, "little"))

    def bit(self, index):
        if index < 0:
            raise ValueError("negative bit position")
        byte, shift = divmod(index, 8)
        return (self.data[byte] >> shift) & 1 if byte < len(self.data) else 0

    def word(self, start, length):
        result = 0
        for bit in range(length):
            result |= self.bit(start + bit) << bit
        return result


def edges(a):
    """Fixed shared edge order, generated with scalar loop counters."""
    for kind in ("p", "q"):
        for row in range(a):
            for column in range(a):
                yield (kind, row), ("mid", column)
    for row in range(a):
        for column in range(a):
            yield ("mid", row), ("ib", column)
    for q in range(a):
        for r in range(a):
            yield ("ib", r), ("carry", q * a + r)
    for index in range(a * a - 1):
        yield ("carry", index), ("carry", index + 1)


def group_of(gate, a, length):
    kind, index = gate
    return (index // a if kind == "carry" else index) // length


def lagrange_inplace(points, values, modulus):
    """Unit-denominator interpolation with two linear coefficient arrays."""
    count = len(points)
    if len(values) != count:
        raise ValueError("inconsistent interpolation input")
    result = [0] * count
    for chosen in range(count):
        current = [0] * count
        current[0] = 1
        degree = 0
        denominator = 1
        for other in range(count):
            if other == chosen:
                continue
            x = points[other]
            # Read the old coefficient at power-1 before overwriting it.
            for power in range(degree + 1, 0, -1):
                current[power] = (current[power - 1] - x * current[power]) % modulus
            current[0] = -x * current[0] % modulus
            degree += 1
            denominator = denominator * (points[chosen] - x) % modulus
        if gcd(denominator, modulus) != 1:
            raise ValueError("interpolation denominator is not a unit")
        scale = values[chosen] * pow(denominator, -1, modulus) % modulus
        for power in range(count):
            result[power] = (result[power] + scale * current[power]) % modulus
    return result


class StreamedTop:
    """One explicitly grouped original top-table query.

    The conceptual arithmetic storage bound is O(A*f), in addition to
    the packed boundary.  Python object overhead is not a claim about
    the finite-alphabet implementation proved in the manuscript.
    """

    def __init__(self, s, groups, group, boundary_tape, toy_base=None):
        self.h = 1 << s
        self.a = 2 * self.h
        self.d = 2 * self.h * self.h
        self.f = 4 * self.d + 1
        self.modulus = (1 << (4 * self.d)) + 1
        self.base = 1 << self.d if toy_base is None else toy_base
        if self.base <= 1 or self.base & (self.base - 1):
            raise ValueError("base is not a positive power of two")
        self.digit_bits = self.base.bit_length() - 1
        self.input_bits = self.a * self.h * self.digit_bits
        if not 2 <= groups <= self.a or self.a % groups:
            raise ValueError("invalid group count")
        if not 0 <= group < groups:
            raise ValueError("invalid group index")
        self.groups, self.group = groups, group
        self.length = self.a // groups
        self.edge_words = 8 * self.length * self.length + 1
        self.edge_bits = self.edge_words * self.f
        self.boundary_bits = (groups - 1) * self.edge_bits
        self.tape = boundary_tape
        self.omega = 16 % self.modulus
        self.inverse = pow(self.omega, -1, self.modulus)
        self.xi = pow(self.omega, self.a, self.modulus)
        self.normalization = pow(self.a * self.a, -1, self.modulus)
        self.cutoff = 0
        self.prefix = None
        self.max_arithmetic_depth = 0
        self.arithmetic_depth = 0
        self.interpolations = 0

    def inside(self, gate):
        return group_of(gate, self.a, self.length) == self.group

    def block_position(self, other):
        if other == self.group:
            raise ValueError("internal edge is not a boundary block")
        return other if other < self.group else other - 1

    def read_boundary(self, inside, outside):
        if not self.inside(inside) or self.inside(outside):
            raise ValueError("incorrect boundary orientation")
        other = group_of(outside, self.a, self.length)
        ordinal = 0
        for source, destination in edges(self.a):
            if self.inside(source) == self.inside(destination):
                continue
            local, remote = (source, destination) if self.inside(source) else (destination, source)
            if group_of(remote, self.a, self.length) != other:
                continue
            if local == inside and remote == outside:
                offset = self.block_position(other) * self.edge_bits + ordinal * self.f
                return (self.tape >> offset) & ((1 << self.f) - 1)
            ordinal += 1
        raise InvalidBoundary("requested edge was absent")

    def validate_padding(self):
        if self.tape < 0 or self.tape.bit_length() > self.boundary_bits:
            raise InvalidBoundary("invalid tape length")
        for other in range(self.groups):
            if other == self.group:
                continue
            words = 8 * self.length * self.length + (abs(other - self.group) == 1)
            used = words * self.f
            padding = self.edge_bits - used
            offset = self.block_position(other) * self.edge_bits + used
            if padding and ((self.tape >> offset) & ((1 << padding) - 1)):
                raise InvalidBoundary("nonzero padding")

    def residue(self, value):
        if not 0 <= value < self.modulus:
            raise InvalidBoundary("illegal arithmetic residue")
        return value

    def enter(self):
        self.arithmetic_depth += 1
        self.max_arithmetic_depth = max(self.max_arithmetic_depth, self.arithmetic_depth)

    def leave(self):
        self.arithmetic_depth -= 1

    def input_digits(self, kind, row):
        self.enter()
        try:
            gate = (kind, row)
            points, values = [], []
            for column in range(self.a):
                if self.inside(("mid", column)):
                    continue
                points.append(pow(self.xi, column, self.modulus))
                value = self.residue(self.read_boundary(gate, ("mid", column)))
                values.append(value * pow(self.inverse, row * column, self.modulus) % self.modulus)
                if len(points) == self.h:
                    break
            if len(points) != self.h:
                raise InvalidBoundary("too few prescribed outputs")
            self.interpolations += 1
            digits = lagrange_inplace(points, values, self.modulus)
            if any(not 0 <= digit < self.base for digit in digits):
                raise InvalidBoundary("input digit outside its range")
            if kind == "p":
                if row == 0 and self.prefix.length > max(0, self.input_bits - self.cutoff):
                    raise InvalidBoundary("prefix is too long")
                for position, value in enumerate(digits):
                    for bit in range(self.digit_bits):
                        index = self.digit_bits * (row + self.a * position) + bit
                        if index >= self.cutoff:
                            if ((value >> bit) & 1) != self.prefix.bit(index - self.cutoff):
                                raise InvalidBoundary("input prefix mismatch")
            return digits
        finally:
            self.leave()

    def input_output(self, kind, row, column):
        digits = self.input_digits(kind, row)
        return pow(self.omega, row * column, self.modulus) * sum(
            value * pow(self.xi, position * column, self.modulus)
            for position, value in enumerate(digits)) % self.modulus

    def middle_inputs(self, column):
        result = []
        for kind in ("p", "q"):
            vector = []
            for row in range(self.a):
                source = (kind, row)
                value = self.input_output(kind, row, column) if self.inside(source) else (
                    self.read_boundary(("mid", column), source))
                vector.append(self.residue(value))
            result.append(vector)
        return result

    def middle_output(self, column, row):
        self.enter()
        try:
            first, second = self.middle_inputs(column)
            total = 0
            for frequency in range(self.a):
                left = sum(value * pow(self.xi, index * frequency, self.modulus)
                           for index, value in enumerate(first)) % self.modulus
                right = sum(value * pow(self.xi, index * frequency, self.modulus)
                            for index, value in enumerate(second)) % self.modulus
                total = (total + left * right * pow(
                    self.inverse, self.a * frequency * row, self.modulus)) % self.modulus
            return pow(self.inverse, column * row, self.modulus) * total % self.modulus
        finally:
            self.leave()

    def inverse_inputs(self, row):
        incoming = []
        for column in range(self.a):
            source = ("mid", column)
            value = self.middle_output(column, row) if self.inside(source) else (
                self.read_boundary(("ib", row), source))
            incoming.append(self.residue(value))
        return incoming

    def inverse_output(self, row, q):
        self.enter()
        try:
            incoming = self.inverse_inputs(row)
            total = sum(value * pow(self.inverse, self.a * column * q, self.modulus)
                        for column, value in enumerate(incoming)) % self.modulus
            return self.normalization * total % self.modulus
        finally:
            self.leave()

    def query(self, number, cutoff, prefix, weight_bit=None, multiply=None, one=1):
        """Return None for an invalid boundary, otherwise its unique weight.

        weight_bit(kind, bit_index, bit_value) may weight either family.
        Arithmetic reconstruction does not invoke this callback.
        """
        if not isinstance(number, BitTape) or not isinstance(prefix, BitTape):
            raise TypeError("N and the prefix must be supplied as read-only bit tapes")
        self.cutoff, self.prefix = cutoff, prefix
        if not 0 <= cutoff <= self.input_bits:
            raise ValueError("invalid prefix arguments")
        if multiply is None:
            multiply = lambda left, right: left * right
        weight = one
        lo, hi = self.group * self.length, (self.group + 1) * self.length
        try:
            self.validate_padding()
            for kind in ("p", "q"):
                for row in range(lo, hi):
                    digits = self.input_digits(kind, row)
                    for column in range(self.a):
                        if self.inside(("mid", column)):
                            continue
                        value = pow(self.omega, row * column, self.modulus) * sum(
                            digit * pow(self.xi, position * column, self.modulus)
                            for position, digit in enumerate(digits)) % self.modulus
                        if value != self.read_boundary((kind, row), ("mid", column)):
                            raise InvalidBoundary("input output mismatch")
                    if weight_bit is not None:
                        for position, digit in enumerate(digits):
                            for bit in range(self.digit_bits):
                                index = self.digit_bits * (row + self.a * position) + bit
                                weight = multiply(weight, weight_bit(kind, index, (digit >> bit) & 1))
            for column in range(lo, hi):
                self.middle_inputs(column)
                for row in range(self.a):
                    if not self.inside(("ib", row)):
                        if self.middle_output(column, row) != self.read_boundary(
                                ("mid", column), ("ib", row)):
                            raise InvalidBoundary("middle output mismatch")
            for row in range(lo, hi):
                self.inverse_inputs(row)
                for q in range(self.a):
                    if not lo <= q < hi:
                        if self.inverse_output(row, q) != self.read_boundary(
                                ("ib", row), ("carry", q * self.a + row)):
                            raise InvalidBoundary("inverse output mismatch")
            first_carry, last_carry = lo * self.a, hi * self.a
            carry = 0 if first_carry == 0 else self.read_boundary(
                ("carry", first_carry), ("carry", first_carry - 1))
            for index in range(first_carry, last_carry):
                q, row = divmod(index, self.a)
                coefficient = self.inverse_output(row, q) if lo <= row < hi else (
                    self.read_boundary(("carry", index), ("ib", row)))
                self.residue(coefficient)
                if not 0 <= carry < 1 << self.f:
                    raise InvalidBoundary("carry outside its range")
                total = coefficient + carry
                expected = number.word(index * self.digit_bits, self.digit_bits)
                if total % self.base != expected:
                    raise InvalidBoundary("incorrect product digit")
                carry = total // self.base
            if last_carry == self.a * self.a:
                if carry:
                    raise InvalidBoundary("nonzero final carry")
            elif carry != self.read_boundary(("carry", last_carry - 1), ("carry", last_carry)):
                raise InvalidBoundary("outgoing carry mismatch")
            return weight
        except InvalidBoundary:
            return None


def pack_boundary(a, groups, group, f, boundary):
    """Checker/input codec; the query itself retains only its packed result."""
    length = a // groups
    block_bits = (8 * length * length + 1) * f
    result = 0
    for other in range(groups):
        if other == group:
            continue
        block = other if other < group else other - 1
        ordinal = 0
        for source, destination in edges(a):
            first = group_of(source, a, length)
            second = group_of(destination, a, length)
            if {first, second} != {group, other}:
                continue
            key = (source, destination) if first == group else (destination, source)
            value = boundary[key]
            if not 0 <= value < 1 << f:
                raise ValueError("boundary word does not fit")
            result |= value << (block * block_bits + ordinal * f)
            ordinal += 1
        assert ordinal == 8 * length * length + (abs(other - group) == 1)
    return result
