#!/usr/bin/env python3
"""Long interval-product jets using parity and exact Newton power sums."""

from math import comb


class QuadraticJets:
    """(Z/2^precision)[w,Z]/(w^2+w+1,Z^(degree+1))."""

    def __init__(self, precision, degree):
        assert precision >= 1 and degree >= 0
        self.precision, self.size = precision, degree + 1
        self.modulus = 1 << precision

    def constant(self, a=0, b=0):
        return ((a % self.modulus, b % self.modulus),) + ((0, 0),) * (self.size - 1)

    def argument(self):
        values = list(self.constant(0, 1))
        if self.size > 1:
            values[1] = (1, 0)
        return tuple(values)

    def pair_product(self, x, y):
        a, b = x
        c, d = y
        return ((a * c - b * d) % self.modulus,
                (a * d + b * c - b * d) % self.modulus)

    def add(self, x, y):
        return tuple(((a + c) % self.modulus, (b + d) % self.modulus)
                     for (a, b), (c, d) in zip(x, y))

    def scale(self, x, scalar):
        return tuple((a * scalar % self.modulus, b * scalar % self.modulus)
                     for a, b in x)

    def multiply(self, x, y):
        result = [(0, 0)] * self.size
        for i, a in enumerate(x):
            for j, b in enumerate(y[:self.size - i]):
                c, d = self.pair_product(a, b)
                u, v = result[i + j]
                result[i + j] = ((u + c) % self.modulus, (v + d) % self.modulus)
        return tuple(result)

    def power(self, x, exponent):
        assert exponent >= 0
        result = self.constant(1)
        while exponent:
            if exponent % 2:
                result = self.multiply(result, x)
            exponent //= 2
            if exponent:
                x = self.multiply(x, x)
        return result

    def inverse(self, x):
        a, b = x[0]
        norm_inverse = pow(a * a - a * b + b * b, -1, self.modulus)
        inverse_constant = ((a - b) * norm_inverse % self.modulus,
                            -b * norm_inverse % self.modulus)
        result = [inverse_constant]
        for n in range(1, self.size):
            total = (0, 0)
            for j in range(1, n + 1):
                c, d = self.pair_product(x[j], result[n - j])
                total = (total[0] + c, total[1] + d)
            c, d = self.pair_product(inverse_constant, total)
            result.append((-c % self.modulus, -d % self.modulus))
        return tuple(result)


def prefix_power_sums(stop, degree):
    sums = []
    for j in range(degree + 1):
        numerator = stop ** (j + 1) - sum(comb(j + 1, i) * sums[i] for i in range(j))
        assert numerator % (j + 1) == 0
        sums.append(numerator // (j + 1))
    return sums


def interval_elementaries(lower, upper, degree):
    """First elementary symmetric functions of lower,...,upper-1, exactly."""
    assert 0 <= lower <= upper and degree >= 0
    high = prefix_power_sums(upper, degree)
    low = prefix_power_sums(lower, degree)
    sums = [a - b for a, b in zip(high, low)]
    elementary = [1]
    for k in range(1, degree + 1):
        numerator = sum((-1) ** (j - 1) * elementary[k - j] * sums[j]
                        for j in range(1, k + 1))
        assert numerator % k == 0
        elementary.append(numerator // k)
    return tuple(elementary)


def unit_progression_jet(ring, interval, unit, step):
    """Product (unit+t*step) for an even jet step and a unit base."""
    start, stop = interval
    assert 0 <= start <= stop and all(a % 2 == b % 2 == 0 for a, b in step)
    length = stop - start
    increment = ring.multiply(step, ring.inverse(unit))
    degree = min(ring.precision - 1, length)
    elementary = interval_elementaries(start, stop, degree)
    series, power = ring.constant(1), ring.constant(1)
    for k in range(1, degree + 1):
        power = ring.multiply(power, increment)
        series = ring.add(series, ring.scale(power, elementary[k]))
    return ring.multiply(ring.power(unit, length), series)


def interval_product_jet(ring, interval, argument=None, number=None):
    """Product (X-a) or (number-a*X), without visiting the interval."""
    lower, upper = interval
    assert 1 <= lower <= upper
    assert number is None or number % 2
    argument = ring.argument() if argument is None else argument
    # Its residue is outside F_2 inside F_4, making both parity bases units.
    assert argument[0][1] % 2
    result = ring.constant(1)
    for parity in (0, 1):
        start = -((parity - lower) // 2)
        stop = -((parity - upper) // 2)
        length = stop - start
        if not length:
            continue
        if number is None:
            unit = ring.add(argument, ring.constant(-parity))
            variable = ring.constant(1)
        else:
            unit = ring.add(ring.constant(number), ring.scale(argument, -parity))
            variable = argument
        result = ring.multiply(result, unit_progression_jet(
            ring, (start, stop), unit, ring.scale(variable, -2)))
    return result
