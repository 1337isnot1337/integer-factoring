#!/usr/bin/env python3
"""Signed global Gaussian Pascal crossings over several prime fields.

This checker keeps ordinary determinants, Pfaffian signs, global mode order,
and both cut profiles explicit. All arithmetic is exact. It makes no claim
of an improved factoring exponent.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path

from check_pascal_complement_convolution import modular_complement


def trim(value, prime):
    value = [x % prime for x in value]
    while value and not value[-1]:
        value.pop()
    return tuple(value)


def add(first, second, prime, coefficient=1):
    result = list(first) + [0] * max(0, len(second) - len(first))
    for i, value in enumerate(second):
        result[i] = (result[i] + coefficient * value) % prime
    return trim(result, prime)


def multiply(first, second, prime, limit=None):
    if not first or not second:
        return ()
    size = len(first) + len(second) - 1
    if limit is not None:
        size = min(size, limit)
    result = [0] * size
    for i, a in enumerate(first):
        if not a or i >= size:
            continue
        for j, b in enumerate(second[:size - i]):
            result[i + j] = (result[i + j] + a * b) % prime
    return trim(result, prime)


def scale(value, coefficient, prime):
    return trim([coefficient * a for a in value], prime)


def exact_divide(numerator, denominator, prime):
    if not denominator:
        raise ZeroDivisionError("zero polynomial denominator")
    if not numerator:
        return ()
    result = [0] * max(0, len(numerator) - len(denominator) + 1)
    remaining = list(numerator)
    inverse = pow(denominator[-1], -1, prime)
    while len(remaining) >= len(denominator):
        shift = len(remaining) - len(denominator)
        coefficient = remaining[-1] * inverse % prime
        result[shift] = coefficient
        for i, value in enumerate(denominator):
            remaining[i + shift] = (remaining[i + shift] - coefficient * value) % prime
        while remaining and not remaining[-1]:
            remaining.pop()
    if remaining:
        raise ArithmeticError("Bareiss division was not exact")
    return trim(result, prime)


def determinant(matrix, prime):
    n = len(matrix)
    if not n:
        return (1,)
    work = [[tuple(value) for value in row] for row in matrix]
    previous = (1,)
    sign = 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if work[i][k]), None)
        if pivot is None:
            return ()
        if pivot != k:
            work[k], work[pivot] = work[pivot], work[k]
            sign = -sign
        current = work[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = add(multiply(current, work[i][j], prime),
                                multiply(work[i][k], work[k][j], prime), prime, -1)
                work[i][j] = exact_divide(numerator, previous, prime)
            work[i][k] = ()
        previous = current
    return scale(work[-1][-1], sign, prime)


def pfaffian(matrix, prime):
    n = len(matrix)
    if not n:
        return (1,)
    if n % 2:
        return ()
    result = ()
    for j in range(1, n):
        if not matrix[0][j]:
            continue
        remaining = [i for i in range(1, n) if i != j]
        minor = [[matrix[a][b] for b in remaining] for a in remaining]
        term = multiply(matrix[0][j], pfaffian(minor, prime), prime)
        result = add(result, term, prime, 1 if j % 2 else -1)
    return result


def binomial_power(exponent, prime):
    return trim([comb(exponent, j) for j in range(exponent + 1)], prime)


def delannoy_polynomial(a, b, prime):
    """(-1)^(a+b) [X^a Y^b] 1/(1-t(X+Y+XY))."""
    result = [0] * (a + b + 1)
    for k in range(max(a, b), a + b + 1):
        result[k] = (-1) ** (a + b) * comb(k, a) * comb(a, k - b)
    return trim(result, prime)


def pascal_kernel(indices, ambient, prime):
    indices = tuple(indices)
    n = len(indices)
    assert ambient % 2 == 0 and ambient > 2 * max(indices)
    matrix = [[()] * (2 * n) for _ in range(2 * n)]
    for i, a in enumerate(indices):
        for j, b in enumerate(indices):
            if a <= b:
                matrix[i][j] = scale(binomial_power(a + 1, prime), comb(b, a), prime)
            if b <= a:
                matrix[i][n + j] = (0,) * ambient + ((-1) ** (a - b) * comb(a, b) % prime,)
                matrix[i][n + j] = trim(matrix[i][n + j], prime)
            if i == j:
                matrix[n + i][j] = (1,)
            o = delannoy_polynomial(a, b, prime)
            coefficients = [0] * ambient
            for degree, value in enumerate(o):
                coefficients[ambient - 1 - degree] = value
            matrix[n + i][n + j] = trim(coefficients, prime)
    return matrix


def normalized_infinity(indices, prime):
    indices = tuple(indices)
    n = len(indices)
    matrix = [[()] * (2 * n) for _ in range(2 * n)]
    for i, a in enumerate(indices):
        for j, b in enumerate(indices):
            if a <= b:
                matrix[i][j] = (0,) * b + scale(binomial_power(a, prime), comb(b, a), prime)
            if b <= a:
                matrix[i][n + j] = trim((0,) * a + ((-1) ** (a - b) * comb(a, b),), prime)
            if i == j:
                matrix[n + i][j] = (0,) * a + (1,)
            matrix[n + i][n + j] = multiply((1, 1), delannoy_polynomial(a, b, prime), prime)
    return matrix


def series_inverse(unit, prime, precision):
    assert unit and unit[0]
    result = [0] * precision
    result[0] = pow(unit[0], -1, prime)
    for degree in range(1, precision):
        result[degree] = -result[0] * sum(unit[j] * result[degree - j]
                                         for j in range(1, min(degree + 1, len(unit)))) % prime
    return trim(result, prime)


def jet_orders(matrix, prime, precision):
    work = [[trim(value[:precision], prime) for value in row] for row in matrix]
    orders = []
    while work:
        candidates = [(next(k for k, a in enumerate(value) if a), i, j)
                      for i, row in enumerate(work) for j, value in enumerate(row) if value]
        if not candidates:
            return None
        valuation, row, column = min(candidates)
        work[0], work[row] = work[row], work[0]
        for values in work:
            values[0], values[column] = values[column], values[0]
        inverse = series_inverse(work[0][0][valuation:], prime, precision - valuation)
        updated = []
        for i in range(1, len(work)):
            quotient = multiply(work[i][0][valuation:], inverse, prime, precision - valuation)
            updated.append([add(work[i][j], multiply(quotient, work[0][j], prime, precision), prime, -1)
                            for j in range(1, len(work))])
        orders.append(valuation)
        work = updated
    assert orders == sorted(orders)
    return tuple(orders)


def smith_orders(matrix, prime, maximum_precision=512):
    precision = 4
    while precision <= maximum_precision:
        orders = jet_orders(matrix, prime, precision)
        if orders is not None:
            return orders, precision
        precision *= 2
    raise ArithmeticError("singular matrix or insufficient exact jet precision")


def weighted_orders(matrix, row_weights, column_weights, prime, infinity=False):
    entries = []
    minimum = 0
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if not value:
                continue
            coefficients = tuple(reversed(value)) if infinity else value
            shift = row_weights[i] + column_weights[j] - (len(value) - 1 if infinity else 0)
            while coefficients and not coefficients[0]:
                shift += 1
                coefficients = coefficients[1:]
            minimum = min(minimum, shift)
            entries.append((i, j, shift, coefficients))
    shifted = [[()] * len(matrix) for _ in matrix]
    for i, j, shift, value in entries:
        shifted[i][j] = (0,) * (shift - minimum) + value
    orders, precision = smith_orders(shifted, prime)
    return tuple(order + minimum for order in orders), precision


def rank_mod(matrix, prime):
    if not matrix:
        return 0
    work = [row[:] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((i for i in range(rank, len(work)) if work[i][column] % prime), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column] % prime, -1, prime)
        work[rank] = [value * inverse % prime for value in work[rank]]
        for i in range(rank + 1, len(work)):
            coefficient = work[i][column] % prime
            if coefficient:
                work[i] = [(a - coefficient * b) % prime for a, b in zip(work[i], work[rank])]
        rank += 1
    return rank


def toeplitz_rank(matrix, prime, length):
    n = len(matrix)
    expanded = [[0] * (n * length) for _ in range(n * length)]
    for degree in range(length):
        for i in range(n):
            for j in range(n):
                value = matrix[i][j][degree] if degree < len(matrix[i][j]) else 0
                if value:
                    for shift in range(length - degree):
                        expanded[(shift + degree) * n + i][shift * n + j] = value
    return rank_mod(expanded, prime)


def type_words(length, counts):
    result = []
    for occupied in product(*(tuple(combinations(range(length), k)) for k in counts)):
        word = [0] * length
        for mode, positions in enumerate(occupied):
            for position in positions:
                word[position] |= 1 << mode
        result.append(tuple(word))
    return result


def selected_positions(word, modes):
    return tuple(mode * len(word) + copy for mode in range(modes)
                 for copy, mask in enumerate(word) if mask >> mode & 1)


def global_kernel(local, copies):
    size = len(local)
    return [[local[i // copies][j // copies] if i % copies == j % copies else ()
             for j in range(size * copies)] for i in range(size * copies)]


def alternating_kernel(matrix, prime):
    n = len(matrix)
    result = [[()] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            result[i][n + j] = matrix[i][j]
            result[n + j][i] = scale(matrix[i][j], -1, prime)
    return result


def evaluate(value, point, prime):
    result = 0
    for coefficient in reversed(value):
        result = (result * point + coefficient) % prime
    return result


def check_type(indices, prime, copies, counts, counts_y=None):
    n = len(indices)
    ambient = 2 * (max(indices) + 1)
    local = pascal_kernel(indices, ambient, prime)
    matrix = global_kernel(local, copies)
    skew = alternating_kernel(matrix, prime)
    boundary_size = n * copies
    counts_y = counts if counts_y is None else counts_y
    labels = type_words(copies, counts)
    labels_y = type_words(copies, counts_y)
    positions = [selected_positions(word, n) for word in labels]
    positions_y = [selected_positions(word, n) for word in labels_y]
    total_x, total_y = sum(counts), sum(counts_y)
    selected_size = total_x + total_y
    endpoint_scalar = (-1) ** (total_x * (total_x - 1) // 2
                              + total_y * (total_y - 1) // 2)
    queries = 0
    nonzero = 0
    for x, y, u, v in product(range(len(labels)), range(len(labels_y)),
                               range(len(labels_y)), range(len(labels))):
        rows = positions[v] + tuple(boundary_size + i for i in positions_y[y])
        columns = positions_y[u] + tuple(boundary_size + i for i in positions[x])
        det = determinant([[matrix[i][j] for j in columns] for i in rows], prime)
        occupied = rows + tuple(2 * boundary_size + i for i in columns)
        pf = pfaffian([[skew[i][j] for j in occupied] for i in occupied], prime)
        assert pf == scale(det, (-1) ** (selected_size * (selected_size - 1) // 2), prime)
        assert evaluate(pf, -1, prime) == endpoint_scalar * int(x == v and y == u) % prime
        queries += 1
        nonzero += bool(pf)
    return {"indices": indices, "prime": prime, "copies": copies, "counts": counts,
            "counts_y": counts_y, "labels": len(labels), "labels_y": len(labels_y),
            "endpoint_scalar": endpoint_scalar % prime,
            "global_endpoint_and_pfaffian_checks": queries,
            "nonzero_polynomial_entries": nonzero}


def ramify_negative(value, prime):
    result = [0] * max(0, 2 * len(value) - 1)
    for degree, coefficient in enumerate(value):
        result[2 * degree] = coefficient * (-1) ** degree
    return trim(result, prime)


def subset_sums(values):
    result = [0]
    for value in values:
        result += [a + value for a in result]
    return tuple(sorted(result))


def check_cut_profiles(indices, prime):
    n = len(indices)
    ambient = 2 * (max(indices) + 1)
    matrix = [[ramify_negative(value, prime) for value in row]
              for row in pascal_kernel(indices, ambient, prime)]
    skew = alternating_kernel(matrix, prime)
    # Original variable order V,Y,B,X. Each displayed cut uses its own
    # contiguous Grassmann order; fixed physical charges make the change
    # back to the typed tensor a constant sign.
    zero_heights = ([0] * n + [2 * a + 1 for a in indices]
                    + [0] * n + [2 * a - 2 * ambient + 1 for a in indices])
    infinity_heights = ([2 * a - n + 2 for a in indices] + [-n] * n
                        + [2 * a - n + 1 for a in indices] + [2 * ambient - n - 1] * n)
    complement = modular_complement(indices, prime)
    infinity_one_particle = [2 * a + 1 - 2 * n for a in sorted(tuple(indices) + complement)]
    records = []
    for infinity, groups, heights, one_particle in (
            (False, (3, 1, 0, 2), zero_heights, [2 * a + 1 for a in indices] * 2),
            (True, (3, 2, 0, 1), infinity_heights, infinity_one_particle)):
        order = tuple(group * n + i for group in groups for i in range(n))
        half = 2 * n
        side_size = 1 << half
        flattening = [[()] * side_size for _ in range(side_size)]
        row_weights, column_weights = [], []
        for mask in range(side_size):
            row_weights.append(sum(heights[order[i]] for i in range(half) if mask >> i & 1))
            column_weights.append(sum(heights[order[half + i]] for i in range(half) if mask >> i & 1))
        for left in range(side_size):
            for right in range(side_size):
                occupied = tuple(order[i] for i in range(half) if left >> i & 1)
                occupied += tuple(order[half + i] for i in range(half) if right >> i & 1)
                flattening[left][right] = pfaffian([[skew[i][j] for j in occupied] for i in occupied], prime)
        orders, precision = weighted_orders(flattening, row_weights, column_weights, prime, infinity)
        assert orders == subset_sums(one_particle)
        records.append({"end": "infinity" if infinity else "zero", "one_particle_orders": one_particle,
                        "full_tensor_orders": orders, "exact_jet_precision": precision})
    return {"indices": indices, "prime": prime, "profiles": records}


def check_naive_power_sign():
    prime = 3
    local = pascal_kernel((0,), 2, prime)
    def endpoint(x, y, u, v):
        rows = ([0] if v else []) + ([1] if y else [])
        columns = ([0] if u else []) + ([1] if x else [])
        if len(rows) != len(columns):
            return 0
        value = determinant([[local[i][j] for j in columns] for i in rows], prime)
        return evaluate(value, -1, prime)
    aligned = endpoint(1, 1, 1, 1) * endpoint(0, 0, 0, 0) % prime
    separated = endpoint(1, 0, 0, 1) * endpoint(0, 1, 1, 0) % prime
    assert aligned == prime - 1 and separated == 1
    return {"prime": prime, "same_common_type": [1], "copies": 2,
            "naive_aligned_endpoint": aligned, "naive_separated_endpoint": separated,
            "ordinary_local_determinant_powers_have_no_common_endpoint_scalar": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    profiles = []
    independent_toeplitz_checks = 0
    for prime in (2, 3, 5, 7):
        for indices in ((0,), (1,), (0, 1), (0, 2), (1, 2), (0, 1, 2), (0, 1, 2, 4)):
            matrix = normalized_infinity(indices, prime)
            complement = modular_complement(indices, prime)
            expected = tuple(sorted(indices + complement))
            actual, precision = smith_orders(matrix, prime)
            assert actual == expected
            if len(indices) <= 2:
                for length in range(1, max(expected) + 2):
                    assert toeplitz_rank(matrix, prime, length) == sum(max(0, length - a) for a in expected)
                    independent_toeplitz_checks += 1
            profiles.append({"indices": indices, "prime": prime, "complement": complement,
                             "normalized_orders": actual, "exact_jet_precision": precision})
    typed = [check_type(indices, prime, copies, counts)
             for prime in (2, 3, 5, 7)
             for indices, copies, counts in (((0,), 3, (1,)), ((0, 2), 2, (1, 1)))]
    paired = [check_type(indices, prime, copies, counts_x, counts_y)
              for prime in (2, 3, 5, 7)
              for indices, copies, counts_x, counts_y in (
                  ((0,), 3, (1,), (2,)), ((0, 2), 2, (1, 0), (1, 1)),
                  ((0, 2), 2, (0, 1), (1, 0)))]
    cuts = [check_cut_profiles(indices, prime)
            for prime in (2, 3, 5, 7) for indices in ((0,), (0, 2))]
    report = {"status": "exact signed endpoints and both Gaussian profiles over fixed prime fields; no improved exponent",
              "normalized_pascal_profiles": profiles,
              "independent_toeplitz_checks": independent_toeplitz_checks,
              "global_typed_checks": typed, "paired_typed_checks": paired,
              "full_gaussian_cut_checks": cuts,
              "negative_control": check_naive_power_sign(),
              "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
