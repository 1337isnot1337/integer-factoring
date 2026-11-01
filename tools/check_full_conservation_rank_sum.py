#!/usr/bin/env python3
"""Exact controls for the characteristic-two full-conservation rank-sum proof.

The finite controls support the accompanying proof; they do not replace it.
No numerical optimizer or third-party package is used.
"""
import argparse
from collections import deque
from functools import lru_cache
import hashlib
import json
from math import comb, isqrt
from pathlib import Path
from random import Random


class Field:
    def __init__(self, order):
        assert order in (2, 3, 4)
        self.order = order

    def add(self, a, b):
        return a ^ b if self.order != 3 else (a + b) % 3

    def sub(self, a, b):
        return a ^ b if self.order != 3 else (a - b) % 3

    def mul(self, a, b):
        if self.order == 4:
            return ((0, 0, 0, 0), (0, 1, 2, 3),
                    (0, 2, 3, 1), (0, 3, 1, 2))[a][b]
        return (a * b) % self.order

    def inv(self, a):
        assert a
        return (0, 1, 3, 2)[a] if self.order == 4 else pow(a, -1, self.order)


def rank(matrix, field):
    rows = [row[:] for row in matrix]
    if not rows:
        return 0
    answer = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(answer, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[answer], rows[pivot] = rows[pivot], rows[answer]
        scalar = field.inv(rows[answer][col])
        rows[answer] = [field.mul(scalar, x) for x in rows[answer]]
        for i in range(answer + 1, len(rows)):
            scalar = rows[i][col]
            if scalar:
                rows[i] = [field.sub(a, field.mul(scalar, b))
                           for a, b in zip(rows[i], rows[answer])]
        answer += 1
        if answer == len(rows):
            break
    return answer


def transpose(matrix):
    return [list(col) for col in zip(*matrix)]


def identity(d):
    return [[int(i == j) for j in range(d)] for i in range(d)]


def matmul(a, b, field):
    assert len(a[0]) == len(b)
    out = [[0] * len(b[0]) for _ in a]
    for i, row in enumerate(a):
        for k, scalar in enumerate(row):
            if scalar:
                for j, value in enumerate(b[k]):
                    out[i][j] = field.add(out[i][j], field.mul(scalar, value))
    return out


def inverse(matrix, field):
    d = len(matrix)
    rows = [matrix[i][:] + identity(d)[i] for i in range(d)]
    for col in range(d):
        pivot = next(i for i in range(col, d) if rows[i][col])
        rows[col], rows[pivot] = rows[pivot], rows[col]
        scalar = field.inv(rows[col][col])
        rows[col] = [field.mul(scalar, x) for x in rows[col]]
        for i in range(d):
            if i != col:
                scalar = rows[i][col]
                rows[i] = [field.sub(a, field.mul(scalar, b))
                           for a, b in zip(rows[i], rows[col])]
    return [row[d:] for row in rows]


def independent_indices(vectors, field):
    chosen = []
    for i in range(len(vectors)):
        if rank([vectors[j] for j in chosen + [i]], field) > len(chosen):
            chosen.append(i)
    return chosen


def actual_basis_factorization(kernel, field):
    """Return K=C M R with actual columns C, actual rows R, invertible M."""
    row_indices = independent_indices(kernel, field)
    column_indices = independent_indices(transpose(kernel), field)
    k = len(row_indices)
    assert len(column_indices) == k
    if not k:
        return ([[] for _ in kernel], [], [], row_indices, column_indices)
    c = [[row[j] for j in column_indices] for row in kernel]
    r = [kernel[i][:] for i in row_indices]
    intersection = [[kernel[i][j] for j in column_indices] for i in row_indices]
    m = inverse(intersection, field)
    assert matmul(matmul(c, m, field), r, field) == kernel
    # Multiplication by the full-column/full-row rank outside factors is injective.
    assert rank(matmul(r, c, field), field) == rank(matmul(kernel, kernel, field), field)
    return c, m, r, row_indices, column_indices


def common_basis(left_vectors, right_vectors, d, field):
    """Exact linear-matroid intersection by shortest augmenting paths.

    This is used only to construct independently checkable Cauchy--Binet
    witnesses. Every returned selection is checked by fresh elimination.
    """
    count = len(left_vectors)
    assert count == len(right_vectors)
    chosen = set()
    caches = ({}, {})

    def independent(indices, side):
        key = tuple(sorted(indices))
        if key not in caches[side]:
            vectors = (left_vectors, right_vectors)[side]
            caches[side][key] = rank([vectors[i] for i in key], field) == len(key)
        return caches[side][key]

    while len(chosen) < d:
        outside = [i for i in range(count) if i not in chosen]
        parents = {i: None for i in outside if independent(chosen | {i}, 0)}
        queue = deque(parents)
        terminal = None
        while queue and terminal is None:
            node = queue.popleft()
            if node not in chosen:
                if independent(chosen | {node}, 1):
                    terminal = node
                    break
                neighbors = [old for old in sorted(chosen)
                             if independent((chosen - {old}) | {node}, 1)]
            else:
                neighbors = [new for new in outside
                             if independent((chosen - {node}) | {new}, 0)]
            for nxt in neighbors:
                if nxt not in parents:
                    parents[nxt] = node
                    queue.append(nxt)
        assert terminal is not None, "Cauchy--Binet guarantees a common basis"
        while terminal is not None:
            if terminal in chosen:
                chosen.remove(terminal)
            else:
                chosen.add(terminal)
            terminal = parents[terminal]
        assert independent(chosen, 0) and independent(chosen, 1)
    result = sorted(chosen)
    assert rank([left_vectors[i] for i in result], field) == d
    assert rank([right_vectors[i] for i in result], field) == d
    return result


def degree(values, width, field):
    coefficients = list(values)
    for bit in range(width):
        for mask in range(1 << width):
            if mask & (1 << bit):
                coefficients[mask] = field.sub(coefficients[mask], coefficients[mask ^ (1 << bit)])
    return max((mask.bit_count() for mask, value in enumerate(coefficients) if value), default=-1)


def convolution_rank(values, field):
    d = len(values)
    return rank([[values[a ^ b] for b in range(d)] for a in range(d)], field)


def cuts(kernel, field):
    return (sum(convolution_rank(row, field) for row in kernel),
            sum(convolution_rank(col, field) for col in transpose(kernel)))


def binary_rank(words):
    pivots = {}
    for word in words:
        while word:
            pivot = word.bit_length() - 1
            if pivot in pivots:
                word ^= pivots[pivot]
            else:
                pivots[pivot] = word
                break
    return len(pivots)


@lru_cache(maxsize=None)
def low_masks(width):
    d = 1 << width
    return tuple(sum(((1 << (1 << bit)) - 1) << start
                     for start in range(0, d, 2 << bit)) for bit in range(width))


@lru_cache(maxsize=None)
def binary_convolution_rank(word, width):
    rows = [word]
    for bit, low in enumerate(low_masks(width)):
        shift = 1 << bit
        rows += [((row & low) << shift) | ((row >> shift) & low) for row in rows]
    return binary_rank(rows)


def binary_degree(word, width):
    for bit, low in enumerate(low_masks(width)):
        word ^= (word & low) << (1 << bit)
    result = -1
    while word:
        low = word & -word
        result = max(result, (low.bit_length() - 1).bit_count())
        word ^= low
    return result


def dense_tensor_cuts(kernel, field):
    d = len(kernel)
    row = [[0] * (d * d) for _ in range(d * d)]
    column = [[0] * (d * d) for _ in range(d * d)]
    for x in range(d):
        for y in range(d):
            for u in range(d):
                v = x ^ y ^ u
                value = kernel[x ^ y][x ^ u]
                row[x * d + y][u * d + v] = value
                column[x * d + u][y * d + v] = value
    return rank(row, field), rank(column, field)


def verify_basis_degrees(vectors, width, field):
    d = 1 << width
    assert len(vectors) == d and rank(vectors, field) == d
    degrees = [degree(v, width, field) for v in vectors]
    for cutoff in range(width + 1):
        assert sum(x <= cutoff for x in degrees) <= sum(comb(width, j) for j in range(cutoff + 1))
    assert sum(degrees) >= width * d // 2
    return degrees


def cauchy_binet_witness(kernels, width, field):
    d = 1 << width
    factors = [actual_basis_factorization(k, field) for k in kernels]
    total = sum(len(item[1]) for item in factors)
    c = [[] for _ in range(d)]
    r = []
    m = [[0] * total for _ in range(total)]
    labels = []
    offset = 0
    for label, (ci, mi, ri, _, _) in enumerate(factors):
        k = len(mi)
        for s in range(d):
            c[s] += ci[s]
        r += ri
        for a in range(k):
            for b in range(k):
                m[offset + a][offset + b] = mi[a][b]
        labels += [label] * k
        offset += k
    mr = matmul(m, r, field)
    assert matmul(c, mr, field) == identity(d)
    a = common_basis(transpose(c), mr, d, field)
    ma = [m[i] for i in a]
    b = common_basis(transpose(ma), r, d, field)
    ca = [[row[j] for j in a] for row in c]
    rb = [r[i] for i in b]
    middle = [[m[i][j] for j in b] for i in a]
    assert rank(ca, field) == rank(rb, field) == rank(middle, field) == d
    column_counts = [sum(labels[j] == i for j in a) for i in range(len(kernels))]
    row_counts = [sum(labels[j] == i for j in b) for i in range(len(kernels))]
    assert column_counts == row_counts and sum(row_counts) == d
    row_degrees = verify_basis_degrees(rb, width, field)
    column_degrees = verify_basis_degrees(transpose(ca), width, field)
    profiles = [cuts(k, field) for k in kernels]
    for label, (rr, cc) in enumerate(profiles):
        assert rr >= sum(1 << row_degrees[t] for t, j in enumerate(b) if labels[j] == label)
        assert cc >= sum(1 << column_degrees[t] for t, j in enumerate(a) if labels[j] == label)
    # A separate, directed exact radical lower bound checks the final conclusion.
    scale = 1 << 40
    lower = sum(isqrt(rr * cc * scale * scale) for rr, cc in profiles)
    target_square = d ** 3 * scale * scale
    target_upper = isqrt(target_square)
    if target_upper * target_upper < target_square:
        target_upper += 1
    assert lower >= target_upper
    return {"counts": row_counts, "row_degree_sum": sum(row_degrees),
            "column_degree_sum": sum(column_degrees), "profiles": profiles}


def random_invertible(d, field, rng):
    out = identity(d)
    for _ in range(8 * d):
        i, j = rng.sample(range(d), 2)
        scalar = rng.randrange(1, field.order)
        out[i] = [field.add(a, field.mul(scalar, b)) for a, b in zip(out[i], out[j])]
    assert rank(out, field) == d
    return out


def square_zero_decomposition(d):
    answer = [[[0] * d for _ in range(d)] for _ in range(3)]
    for j in range(d // 2):
        a, b = 2 * j, 2 * j + 1
        answer[0][a][b] = answer[1][b][a] = 1
        for u in (a, b):
            for v in (a, b):
                answer[2][u][v] = 1
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = Random(20261030127)
    result = {'status': 'PASS', 'scope': 'Characteristic-two fixed additive conservation; no factoring exponent improvement'}

    count = 0
    for width in range(1, 5):
        for word in range(1 << (1 << width)):
            j = binary_degree(word, width)
            assert binary_convolution_rank(word, width) >= (0 if j < 0 else 1 << j)
            count += 1
    result['exhaustive_binary_function_degree_controls'] = count

    f4 = Field(4)
    count = 0
    for value in range(4 ** 4):
        function = [(value >> (2 * i)) & 3 for i in range(4)]
        j = degree(function, 2, f4)
        assert convolution_rank(function, f4) >= (0 if j < 0 else 1 << j)
        count += 1
    for width in (3, 4):
        for _ in range(40):
            function = [rng.randrange(4) for _ in range(1 << width)]
            j = degree(function, width, f4)
            assert convolution_rank(function, f4) >= 1 << j
            count += 1
    result['extension_field_function_degree_controls'] = count

    count = 0
    for value in range(1 << 16):
        rows = [(value >> (4 * i)) & 15 for i in range(4)]
        columns = [sum(((rows[i] >> j) & 1) << i for i in range(4)) for j in range(4)]
        square = []
        for row in rows:
            output = 0
            for j in range(4):
                if row >> j & 1:
                    output ^= rows[j]
            square.append(output)
        rr = sum(binary_convolution_rank(row, 2) for row in rows)
        cc = sum(binary_convolution_rank(col, 2) for col in columns)
        assert rr * cc >= 4 * binary_rank(square) ** 2
        count += 1
    result['exhaustive_binary_rank_square_controls'] = count

    witnesses = []
    arbitrary_count = nilpotent_count = dense_count = rank_square_count = 0
    for order in (2, 4):
        field = Field(order)
        for width in (1, 2, 3):
            d = 1 << width
            for trial in range(8):
                # At least three genuinely non-idempotent summands, without support restrictions.
                while True:
                    kernels = [[[rng.randrange(order) for _ in range(d)] for _ in range(d)] for _ in range(2)]
                    final = identity(d)
                    for kernel in kernels:
                        final = [[field.sub(a, b) for a, b in zip(row, other)]
                                 for row, other in zip(final, kernel)]
                    kernels.append(final)
                    if all(matmul(k, k, field) != k for k in kernels):
                        break
                witness = cauchy_binet_witness(kernels, width, field)
                arbitrary_count += 1
                if trial == 0:
                    witnesses.append({'field_order': order, 'width': width, 'family': 'arbitrary', **witness})
                for kernel in kernels:
                    rr, cc = cuts(kernel, field)
                    assert rr * cc >= d * rank(matmul(kernel, kernel, field), field) ** 2
                    rank_square_count += 1
                if trial < 2 and width <= 2:
                    for kernel in kernels:
                        assert dense_tensor_cuts(kernel, field) == cuts(kernel, field)
                        dense_count += 1
            for trial in range(4):
                p = random_invertible(d, field, rng)
                pinv = inverse(p, field)
                kernels = [matmul(matmul(p, k, field), pinv, field) for k in square_zero_decomposition(d)]
                assert all(rank(k, field) == d // 2 and rank(matmul(k, k, field), field) == 0 for k in kernels)
                witness = cauchy_binet_witness(kernels, width, field)
                nilpotent_count += 1
                if trial == 0:
                    witnesses.append({'field_order': order, 'width': width, 'family': 'three_square_zero', **witness})
    result['arbitrary_nonprojector_decompositions'] = arbitrary_count
    result['conjugated_three_square_zero_decompositions'] = nilpotent_count
    result['random_rank_square_controls'] = rank_square_count
    result['independent_dense_tensor_controls'] = dense_count
    result['exact_double_cauchy_binet_witnesses'] = arbitrary_count + nilpotent_count
    result['sample_witnesses'] = witnesses

    field = Field(2)
    d = 16
    dot = [[(s & w).bit_count() & 1 for w in range(d)] for s in range(d)]
    assert cuts(dot, field) == (30, 30) and rank(dot, field) == 4
    assert rank(matmul(dot, dot, field), field) == 0
    assert 30 * 30 < d * d * 4
    result['false_rank_product_negative_control'] = {'labels': d, 'cuts': [30, 30], 'rank_K': 4, 'rank_K_squared': 0}

    field = Field(3)
    d = 4
    walsh = [[1 if not ((s & w).bit_count() & 1) else 2 for w in range(d)] for s in range(d)]
    assert rank(walsh, field) == d and rank(matmul(walsh, walsh, field), field) == d
    assert cuts(walsh, field) == (d, d)
    assert d * d < d * d * d
    result['odd_characteristic_negative_control'] = {'field_order': 3, 'labels': d, 'cuts': [d, d], 'rank_K_squared': d,
                                                   'scope': 'Refutes characteristic-free degree/rank-square and invertible-target extensions'}
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
