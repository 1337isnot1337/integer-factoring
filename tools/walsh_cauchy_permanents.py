#!/usr/bin/env python3
"""Exact principal permanents and finite scalar-decoder building blocks."""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import factorial

from cauchy_packet_polynomials import cauchy_polynomials, multiply, shift


def walsh_transform(values):
    result = list(values)
    count = len(result)
    if not count or count & (count - 1):
        raise ValueError('The Walsh transform length must be a power of two.')
    step = 1
    while step < count:
        for start in range(0, count, 2 * step):
            for index in range(start, start + step):
                left, right = result[index], result[index + step]
                result[index], result[index + step] = left + right, left - right
        step *= 2
    return result


def permanent(matrix, modulus):
    count = len(matrix)
    values = [0] * (1 << count)
    values[0] = 1
    for state in range(1, 1 << count):
        row = state.bit_count() - 1
        total = 0
        bits = state
        while bits:
            bit = bits & -bits
            bits ^= bit
            column = bit.bit_length() - 1
            total += matrix[row][column] * values[state ^ bit]
        values[state] = total % modulus
    return values[-1]


def principal_permanents(matrix, modulus):
    length = len(matrix)
    values = []
    for state in range(1 << length):
        selected = [x for x in range(length) if state >> x & 1]
        values.append(permanent([[matrix[x][y] for y in selected] for x in selected], modulus))
    return values


def principal_energies(length, center, subset):
    """Linear assignment cost and exact quadratic Cauchy determinant order."""
    a = length.bit_length() - 1
    selected = [x for x in range(length) if subset >> x & 1]
    absolute = square = 0
    for level in range(1, a + 2):
        period = 1 << level
        counts = [0] * period
        for x in selected:
            counts[x % period] += 1
        differences = [counts[r] - counts[(r - center) % period] for r in range(period)]
        absolute += sum(abs(value) for value in differences)
        square += sum(value * value for value in differences)
    assert not absolute % 2 and not square % 2
    return absolute // 2, square // 2


def minimum_assignment(weights):
    count = len(weights)
    costs = [None] * (1 << count)
    costs[0] = 0
    for state in range(1, 1 << count):
        row = state.bit_count() - 1
        bits = state
        best = None
        while bits:
            bit = bits & -bits
            bits ^= bit
            column = bit.bit_length() - 1
            value = costs[state ^ bit] + weights[row][column]
            best = value if best is None else min(best, value)
        costs[state] = best
    return costs[-1]


def rational_residue(value, modulus):
    value = Fraction(value)
    if not value.denominator & 1:
        raise ValueError('The rational constant is not two-adically integral.')
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


@lru_cache(maxsize=None)
def global_kernels(length):
    return cauchy_polynomials(length)


def local_kernel_polynomials(length, center, precision):
    """Independent global-polynomial construction, used for exact replay."""
    modulus = 1 << precision
    return {difference: [rational_residue(value, modulus)
                         for value in shift(polynomial, center, 2 * length)]
            for difference, polynomial in global_kernels(length).items()}


def local_kernel_jets(length, center, precision):
    """Compute the first precision coefficients from the local root product.

    Every coefficient of degree j is divisible by 2^j. Therefore these
    coefficients determine the complete polynomial modulo 2^precision.
    """
    modulus = 1 << precision
    dilation = 2 * length
    derivative_ratio = Fraction(
        (-1 if center % 2 else 1) * factorial(length + center) * factorial(length - 1 - center),
        factorial(length) * factorial(length - 1))
    unit = rational_residue(derivative_ratio, modulus)
    product = [1] + [0] * (precision - 1)
    factors = {}
    for root in range(-length, length):
        if root == center:
            continue
        factor = rational_residue(Fraction(dilation, center - root), modulus)
        factors[root] = factor
        for degree in range(precision - 1, 0, -1):
            product[degree] = (product[degree] + factor * product[degree - 1]) % modulus
    kernels = {}
    for difference in range(-length + 1, length):
        if difference == center:
            values = [unit * value % modulus for value in product]
        else:
            factor = factors[difference]
            quotient = [product[0]]
            for degree in range(1, precision):
                quotient.append((product[degree] - factor * quotient[-1]) % modulus)
            values = [0] + [unit * factor * value % modulus for value in quotient[:-1]]
        kernels[difference] = values
    return kernels


def principal_permanent_polynomials(length, center, precision, subsets=None):
    """Full local polynomials for selected principal subsets, by subset DP."""
    modulus = 1 << precision
    kernels = local_kernel_polynomials(length, center, precision)
    if subsets is None:
        subsets = range(1 << length)
    result = {}
    for subset in subsets:
        selected = [x for x in range(length) if subset >> x & 1]
        count = len(selected)
        values = [None] * (1 << count)
        values[0] = [1]
        for state in range(1, 1 << count):
            row = selected[state.bit_count() - 1]
            total = [0] * (state.bit_count() * (2 * length - 1) + 1)
            bits = state
            while bits:
                bit = bits & -bits
                bits ^= bit
                column = selected[bit.bit_length() - 1]
                term = multiply(values[state ^ bit], kernels[column - row], modulus)
                for degree, value in enumerate(term):
                    total[degree] += value
            values[state] = [value % modulus for value in total]
        result[subset] = values[-1]
    return result


def multiply_truncated(left, right, modulus, precision):
    result = [0] * precision
    for i, a in enumerate(left[:precision]):
        if a:
            for j, b in enumerate(right[:precision - i]):
                if b:
                    result[i + j] += a * b
    return [value % modulus for value in result]


def sparse_interval_jets(length, center, precision):
    """All initial-interval principal permanents by a precision-pruned DP.

    The routine has no asserted polynomial worst-case bound. It records
    the number of live subset states so finite search costs stay visible.
    """
    if center <= 0 or center % 2 == 0:
        raise ValueError('This routine expects a positive odd center.')
    count = length if center < precision else min(length, precision - 1)
    modulus = 1 << precision
    mask = modulus - 1
    kernels = {difference: [(degree, value) for degree, value in enumerate(polynomial) if value]
               for difference, polynomial in local_kernel_jets(length, center, precision).items()}
    states = {0: [1] + [0] * (precision - 1)}
    answers = [[1] + [0] * (precision - 1)]
    live_counts = []
    for row in range(count):
        following = {}
        for state, polynomial in states.items():
            unused = ((1 << count) - 1) ^ state
            nonzero = [(degree, value) for degree, value in enumerate(polynomial) if value]
            while unused:
                bit = unused & -unused
                unused ^= bit
                column = bit.bit_length() - 1
                value = [0] * precision
                for j, a in kernels[column - row]:
                    for i, b in nonzero:
                        if i + j < precision:
                            value[i + j] += a * b
                value = [entry & mask for entry in value]
                if any(value):
                    target = state | bit
                    if target in following:
                        following[target] = [(a + b) & mask for a, b in zip(following[target], value)]
                    else:
                        following[target] = value
        states = {state: value for state, value in following.items() if any(value)}
        answers.append(states.get((1 << (row + 1)) - 1, [0] * precision))
        live_counts.append(len(states))
    # If center>=precision, a principal interval of size>=precision
    # has at least precision non-shift entries in every permutation.
    answers.extend([[0] * precision for _ in range(length - count)])
    return answers, live_counts


def reflect_interval_jets(length, center, precision, positive):
    """Use the exact Cauchy reflection identity to obtain the negative disk."""
    modulus = 1 << precision
    ratio = [rational_residue(Fraction(length - center, length + center), modulus)]
    ratio.extend(rational_residue(Fraction(2 * length, length + center)**(degree + 1), modulus)
                 for degree in range(1, precision))
    power = [1] + [0] * (precision - 1)
    result = []
    for polynomial in positive:
        reflected = [value * (-1)**degree % modulus for degree, value in enumerate(polynomial)]
        result.append(multiply_truncated(reflected, power, modulus, precision))
        power = multiply_truncated(power, ratio, modulus, precision)
    return result


def interval_first_three_jets(length, precision):
    """Exact jets through order three at centers +/-1 and +/-3.

    The formula enumerates only permutations with at most three broken
    shift edges. It has polynomial cost in length, independent of 2^length.
    See the note for its three-break terms and the small direct controls.
    """
    if length < 4 or length & (length - 1) or precision < 4:
        raise ValueError('Use a power-of-two length at least four and precision at least four.')
    modulus = 1 << precision
    alpha = 2 * length
    reciprocal = [0] + [rational_residue(Fraction(alpha, d), modulus)
                        for d in range(1, length + 3)]
    valuations = [precision] + [(value & -value).bit_length() - 1 if value else precision
                                for value in reciprocal[1:]]
    square_convolution = [0] * (length + 1)
    cube_convolution = [0] * (length + 1)
    mixed = [0] * (length + 1)
    three_break = [0] * (length + 1)
    for k in range(2, length + 1):
        square_convolution[k] = sum(reciprocal[j] * reciprocal[k - j]
                                    for j in range(1, k)) % modulus
        mixed[k] = 2 * sum(reciprocal[j]**2 * reciprocal[k - j]
                           for j in range(1, k)) % modulus
    for k in range(3, length + 1):
        cube_convolution[k] = sum(reciprocal[j] * square_convolution[k - j]
                                  for j in range(1, k)) % modulus
    for middle in range(1, length - 1):
        for left in range(middle + 1, length):
            # The remaining denominator is below length, so its reciprocal
            # has valuation at least two. This only prunes zero terms.
            if valuations[middle] + valuations[left] + 2 >= precision:
                continue
            for right in range(left, min(length - 1, length + middle - left) + 1):
                if valuations[middle] + valuations[left] + valuations[right] < precision:
                    k = left + right - middle
                    three_break[k] -= (reciprocal[middle] * reciprocal[left] * reciprocal[right]
                                       * (1 if left == right else 2))
    three_break = [value % modulus for value in three_break]
    shift = local_kernel_jets(length, 1, precision)[1][:4]
    power = [1, 0, 0, 0]
    positive = [[1, 0, 0, 0]]
    for k in range(1, length + 1):
        power = multiply_truncated(power, shift, modulus, 4)
        tail = sum((k - j - 1) * reciprocal[j]**2 for j in range(1, k - 1)) % modulus
        reduced = [0, reciprocal[k],
                   (-reciprocal[k]**2 + square_convolution[k]) % modulus,
                   (reciprocal[k]**3 - mixed[k] + cube_convolution[k]
                    + three_break[k] - reciprocal[k] * tail) % modulus]
        positive.append(multiply_truncated(power, reduced, modulus, 4))

    def reflect(center, values):
        return [row[:4] for row in reflect_interval_jets(length, center, precision, values)]

    third_kernels = local_kernel_jets(length, 3, precision)
    unit = third_kernels[3][0]
    unit_power = 1
    positive_third = [[1, 0, 0, 0]]
    for k in range(1, length + 1):
        unit_power = unit_power * unit % modulus
        if k == 1:
            polynomial = third_kernels[0][:4]
        elif k == 2:
            left = multiply_truncated(third_kernels[0][:4], third_kernels[0][:4], modulus, 4)
            right = multiply_truncated(third_kernels[1][:4], third_kernels[-1][:4], modulus, 4)
            polynomial = [(x + y) % modulus for x, y in zip(left, right)]
        else:
            residual = [[reciprocal[k + i - j] for j in range(3)] for i in range(3)]
            polynomial = [0, 0, 0, unit_power * permanent(residual, modulus) % modulus]
        positive_third.append(polynomial)
    return {1: positive, -1: reflect(1, positive),
            3: positive_third, -3: reflect(3, positive_third)}


def restricted_interval_features(length, precision):
    jets = interval_first_three_jets(length, precision)
    features = [(-1, 0), (-1, 1), (-1, 2), (-3, 3), (-1, 3), (1, 1)]
    rows = [[int(k == 0) for k in range(length + 1)]]
    rows.extend([[jets[center][k][order] for k in range(length + 1)]
                 for center, order in features[1:]])
    return features, rows
