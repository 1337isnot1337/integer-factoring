"""Exact monotone retiming of scalar coefficient order/degree profiles.

The optimization applies to constant-basis diagonal cuts. Polynomial
Smith bases are not constant bases and require a fresh physical check.
All block comparisons and retained pressure bounds use rational integers.
"""

from fractions import Fraction
from itertools import product
from math import isqrt


def optimal_blocks(low, high):
    """Return an exact first-order certificate for the global infimum.

    Each tuple is (first index, last index, low total, high total).
    The infimum pressure is (sum sqrt(low total * high total))**2.
    """
    low, high = tuple(low), tuple(high)
    if len(low) != len(high) or not low or low[0] <= 0 or high[-1] <= 0:
        raise ValueError("normalize the common monomial before retiming")
    if any(not isinstance(n, int) or n < 0 for n in low + high):
        raise ValueError("rank multiplicities must be nonnegative integers")
    blocks = []
    for index, (a, b) in enumerate(zip(low, high)):
        if a + b == 0:
            continue
        blocks.append([index, index, a, b])
        while len(blocks) > 1 and blocks[-2][2] * blocks[-1][3] > blocks[-1][2] * blocks[-2][3]:
            right = blocks.pop()
            left = blocks[-1]
            left[1] = right[1]
            left[2] += right[2]
            left[3] += right[3]
    for index, block in enumerate(blocks):
        block[1] = blocks[index + 1][0] - 1 if index + 1 < len(blocks) else len(low) - 1
        assert block[2] > 0 and block[3] > 0
    return tuple(tuple(block) for block in blocks)


def root_interval(blocks, precision=40):
    """Rational outward interval for the optimal pressure square root."""
    denominator = 1 << precision
    lower = upper = 0
    for _, _, a, b in blocks:
        square = a * b * denominator ** 2
        floor = isqrt(square)
        lower += floor
        upper += floor + (floor * floor != square)
    return Fraction(lower, denominator), Fraction(upper, denominator)


def nearest_half_log(ratio, scale):
    """Nearest integer to log(ratio)/(2 log(scale)), using exact powers."""
    ratio, scale = Fraction(ratio), Fraction(scale)
    if ratio < 1 or scale <= 1:
        raise ValueError("requires ratio >= 1 and scale > 1")
    low, high = -1, 1
    while scale ** (2 * high + 1) < ratio:
        high *= 2
    while high - low > 1:
        middle = (low + high) // 2
        if scale ** (2 * middle + 1) < ratio:
            low = middle
        else:
            high = middle
    return high


def profile_pressure(low, high, degrees, scale):
    scale = Fraction(scale)
    if len(low) != len(high) or len(low) != len(degrees) or scale < 1:
        raise ValueError("incompatible profiles, schedule, or scale")
    a = sum(Fraction(n) * scale ** (-i) for n, i in zip(low, degrees))
    b = sum(Fraction(n) * scale ** i for n, i in zip(high, degrees))
    return a * b


def rational_schedule(low, high, accuracy=64):
    """Strict integer schedule within x**(D+1) of the optimal pressure.

    Here x=1+1/accuracy and D is the original degree. Reciprocity is not
    required. Higher accuracy can produce large exact rational receipts.
    """
    if not isinstance(accuracy, int) or accuracy < 1:
        raise ValueError("accuracy must be a positive integer")
    blocks = optimal_blocks(low, high)
    scale = Fraction(accuracy + 1, accuracy)
    first_ratio = Fraction(blocks[0][2], blocks[0][3])
    degrees = [None] * len(low)
    for start, stop, a, b in blocks:
        rounded = nearest_half_log(Fraction(a, b) / first_ratio, scale)
        for index in range(start, stop + 1):
            degrees[index] = rounded + index
    assert degrees[0] == 0
    assert all(a < b for a, b in zip(degrees, degrees[1:]))
    bound = profile_pressure(low, high, degrees, scale)
    factor = scale ** len(low)
    lower, upper = root_interval(blocks)
    return {
        "blocks": [list(block) for block in blocks],
        "retimed_degrees": degrees,
        "scale": [scale.numerator, scale.denominator],
        "pressure": [bound.numerator, bound.denominator],
        "approximation_factor": [factor.numerator, factor.denominator],
        "infimum_root_interval": [[lower.numerator, lower.denominator],
                                   [upper.numerator, upper.denominator]],
    }


def retime_coefficients(rows, degrees):
    if len(rows) != len(degrees) or degrees[0] != 0:
        raise ValueError("incompatible degree schedule")
    if any(a >= b for a, b in zip(degrees, degrees[1:])):
        raise ValueError("degree indices must be strictly increasing")
    result = [(0,) * len(rows[0]) for _ in range(degrees[-1] + 1)]
    for index, row in zip(degrees, rows):
        result[index] = tuple(row)
    return tuple(result)


def retime_binary_polynomial(value, degrees):
    """Map a polynomial bit mask over F2 by an injective degree schedule."""
    if value < 0 or value.bit_length() > len(degrees):
        raise ValueError("polynomial exceeds the degree schedule")
    if any(a >= b for a, b in zip(degrees, degrees[1:])):
        raise ValueError("degree indices must be strictly increasing")
    answer = 0
    while value:
        bit = value & -value
        answer |= 1 << degrees[bit.bit_length() - 1]
        value ^= bit
    return answer


def physical_binary_profiles(values, alphabet, degrees, zero_heights,
                             infinity_heights, scalar_shifts=(0, 0)):
    """Recompute both complete physical Smith profiles after retiming.

    Dense tensor entries use x+q*y+q**2*u+q**3*v. Each height argument
    contains four lists of integer physical heights. Opposite scalar shifts
    cancel and are permitted explicitly; no old Smith order is retimed.
    """
    from check_block_smith_crossings import reverse_polynomial, smith_profile

    q = alphabet
    if not isinstance(q, int) or q < 2 or len(values) != q ** 4:
        raise ValueError("incompatible physical alphabet")
    if len(scalar_shifts) != 2 or sum(scalar_shifts) != 0:
        raise ValueError("the displayed common scalar shifts must cancel")
    for heights in (zero_heights, infinity_heights):
        if len(heights) != 4 or any(len(row) != q for row in heights):
            raise ValueError("each physical leg needs one height per letter")
        if any(not isinstance(value, int) for row in heights for value in row):
            raise ValueError("use integer heights in the exact physical replay")
    changed = tuple(retime_binary_polynomial(value, degrees) for value in values)
    degree = max(value.bit_length() - 1 for value in changed)
    for x, y, u, v in product(range(q), repeat=4):
        value = changed[x + q * y + q * q * u + q ** 3 * v]
        if value.bit_count() % 2 != int(x == v and y == u):
            raise ValueError("retimed tensor does not have the strict crossing endpoint")
    profiles, metadata = [], []
    for infinity, heights, scalar in ((False, zero_heights, scalar_shifts[0]),
                                      (True, infinity_heights, scalar_shifts[1])):
        entries = []
        for labels in product(range(q), repeat=4):
            x, y, u, v = labels
            value = changed[x + q * y + q * q * u + q ** 3 * v]
            if not value:
                continue
            shift = scalar + sum(row[label] for row, label in zip(heights, labels))
            if infinity:
                value = reverse_polynomial(value, degree)
                shift -= degree
                row, column = q * x + u, q * y + v
            else:
                row, column = q * x + y, q * u + v
            order = (value & -value).bit_length() - 1
            entries.append((row, column, value >> order, shift + order))
        offset = min(entry[3] for entry in entries)
        matrix = [[0] * (q * q) for _ in range(q * q)]
        for row, column, unit, shift in entries:
            matrix[row][column] = unit << (shift - offset)
        maximum = max(value.bit_length() - 1 for row in matrix for value in row)
        orders, precision = smith_profile(matrix, maximum)
        assert len(orders) == q * q
        profiles.append([order + offset for order in orders])
        metadata.append({"clearing_offset": offset, "matrix_degree_bound": maximum,
                         "series_precision": precision})
    return {"degree_schedule": list(degrees), "actual_polynomial_degree": degree,
            "zero_orders": profiles[0], "infinity_orders": profiles[1],
            "profile_metadata": metadata, "literal_endpoint_entries": q ** 4,
            "total_physical_heights": [sum(heights[leg][letter]
                                           for heights in (zero_heights, infinity_heights)
                                           for leg in range(4)) for letter in range(q)]}
