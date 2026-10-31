#!/usr/bin/env python3
"""Outward dyadic bounds for the periodic charged Gaussian continuation.

All arithmetic used by the certificate is integer arithmetic.  Exp and
log bounds use explicit Taylor remainders; no floating result proves a sign.
"""
from functools import lru_cache
from math import isqrt


class Dyadic:
    def __init__(self, precision=80):
        self.precision = precision
        self.scale = 1 << precision
        self.one = (self.scale, self.scale)
        self.zero = (0, 0)
        root = isqrt(2 * self.scale * self.scale)
        self.root_two = (root, root + 1)

    @staticmethod
    def ceil_div(a, b):
        return -((-a) // b)

    def rational(self, a, b=1):
        assert b > 0
        return a * self.scale // b, self.ceil_div(a * self.scale, b)

    @staticmethod
    def add(a, b):
        return a[0] + b[0], a[1] + b[1]

    @staticmethod
    def subtract(a, b):
        return a[0] - b[1], a[1] - b[0]

    def multiply(self, a, b):
        values = [x * y for x in a for y in b]
        return min(values) // self.scale, self.ceil_div(max(values), self.scale)

    def divide_integer(self, a, b):
        assert b > 0
        return a[0] // b, self.ceil_div(a[1], b)

    def inverse(self, a):
        assert 0 < a[0] <= a[1]
        square = self.scale * self.scale
        return square // a[1], self.ceil_div(square, a[0])

    def divide(self, a, b):
        return self.multiply(a, self.inverse(b))

    @lru_cache(maxsize=262144)
    def exp_negative_point(self, numerator):
        assert numerator >= 0
        shifts = 0
        while numerator > (self.scale << shifts) // 2:
            shifts += 1
        y = (numerator // (1 << shifts), self.ceil_div(numerator, 1 << shifts))
        term = total = self.one
        for index in range(1, 25):
            term = self.divide_integer(self.multiply(term, y), index)
            total = self.add(total, term)
        following = self.divide_integer(self.multiply(term, y), 25)
        # y<=1/2: the omitted positive exponential tail is <2 times
        # its first term.  Every preceding operation was outward rounded.
        result = total[0], total[1] + 2 * following[1]
        for _ in range(shifts):
            result = self.multiply(result, result)
        return self.inverse(result)

    def exp_negative(self, a):
        assert a[0] >= 0
        return self.exp_negative_point(a[1])[0], self.exp_negative_point(a[0])[1]

    def log_one_plus(self, a):
        assert a[0] >= 0
        z = self.divide(a, self.add(self.rational(2), a))
        square = self.multiply(z, z)
        assert square[1] < self.scale
        power = z
        total = self.zero
        for index in range(32):
            total = self.add(total, self.divide_integer(power, 2 * index + 1))
            power = self.multiply(power, square)
        remainder = self.divide_integer(
            self.divide(power, self.subtract(self.one, square)), 65)
        return 2 * total[0], 2 * (total[1] + remainder[1])

    @lru_cache(maxsize=262144)
    def softplus_point(self, numerator):
        thermal = self.log_one_plus(self.exp_negative_point(abs(numerator)))
        return (thermal[0] - min(numerator, 0),
                thermal[1] - min(numerator, 0))

    def softplus(self, a):
        # log(1+exp(-x)) is decreasing.
        return self.softplus_point(a[1])[0], self.softplus_point(a[0])[1]

    def paired_sum(self, scale, charge, terms):
        total = self.zero
        for index in range(1, terms + 1):
            odd = self.rational(2 * index - 1)
            for value in (self.add(odd, charge), self.subtract(odd, charge)):
                total = self.add(total, self.softplus(self.multiply(value, scale)))
        return total

    def continuation(self, scale, charge, terms=20):
        positive = self.add(self.paired_sum(scale, charge, terms),
                            self.paired_sum(scale, self.zero, terms))
        physical_scale = self.multiply(self.root_two, scale)
        half_charge = self.divide_integer(charge, 2)
        negative = self.multiply(self.multiply(self.rational(2), self.root_two),
                                 self.paired_sum(physical_scale, half_charge, terms))
        # The two omitted physical tails have first argument at least
        # sqrt(2)*t*(2N+1-|J|/2), then step 2sqrt(2)*t.
        assert 0 <= charge[0] <= charge[1] <= 2 * self.scale
        first = self.multiply(physical_scale,
                              self.subtract(self.rational(2 * terms + 1), half_charge))
        denominator = self.subtract(self.one, self.exp_negative(
            self.multiply(self.rational(2), physical_scale)))
        tail = self.multiply(self.multiply(self.rational(4), self.root_two),
                             self.divide(self.exp_negative(first), denominator))
        partial = self.subtract(positive, negative)
        # Positive omitted tails can be discarded for the required lower bound.
        return partial[0] - tail[1]

    def finite_pressure(self, n, scale, charge):
        """The complete binary critical pressure, with no tail truncation."""
        result = self.zero
        for index in range(n):
            odd = self.rational(2 * index + 1)
            for value in (self.add(odd, charge), self.subtract(odd, charge)):
                result = self.add(result, self.softplus(self.multiply(value, scale)))
        for index in range(2 * n):
            energy = self.rational(2 * index + 1 - 2 * n)
            result = self.add(result, self.softplus(self.multiply(energy, scale)))
        physical = self.zero
        physical_scale = self.divide(scale, self.root_two)
        for index in range(n):
            energy = self.rational(4 * index - 2 * (n - 1))
            for value in (self.add(energy, charge), self.subtract(energy, charge)):
                physical = self.add(physical, self.softplus(self.multiply(value, physical_scale)))
        return self.subtract(result, self.multiply(self.root_two, physical))


def split_rectangle(rectangle):
    lo_t, hi_t, lo_j, hi_j = rectangle
    # Compare relative scale width with charge width.  This choice affects
    # only the size of the finite certificate, not its validity.
    scale_axis = 2 * (hi_t - lo_t) >= (hi_j - lo_j) * (hi_t + lo_t)
    if scale_axis:
        middle = (lo_t + hi_t) / 2
        return ('0', (lo_t, middle, lo_j, hi_j)), ('1', (middle, hi_t, lo_j, hi_j))
    middle = (lo_j + hi_j) / 2
    return ('2', (lo_t, hi_t, lo_j, middle)), ('3', (lo_t, hi_t, middle, hi_j))
