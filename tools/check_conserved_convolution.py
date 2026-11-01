#!/usr/bin/env python3
"""Exact conserved-convolution blocks, rank bounds, and Smith controls."""

import argparse
from fractions import Fraction
from itertools import product
import json
import random

from crossing_rank_sum import binary_rank


class AbelianGroup:
    def __init__(self, moduli):
        self.moduli = tuple(moduli)
        self.elements = tuple(product(*(range(m) for m in moduli)))
        self.index = {a: i for i, a in enumerate(self.elements)}
        self.order = len(self.elements)
        self.add = tuple(tuple(self.index[tuple((x + y) % m for x, y, m in zip(a, b, moduli))]
                               for b in self.elements) for a in self.elements)
        self.neg = tuple(self.index[tuple((-x) % m for x, m in zip(a, moduli))] for a in self.elements)
        self.double_image = frozenset(self.add[i][i] for i in range(self.order))

    def sub(self, a, b):
        return self.add[a][self.neg[b]]


def matrices(group, kernels):
    d = group.order
    row = [0] * (d * d)
    column = [0] * (d * d)
    blocks = set()
    for w in range(d):
        for difference in range(d):
            if not ((kernels[w] >> difference) & 1):
                continue
            for x in range(d):
                v = group.sub(x, difference)
                y = group.sub(w, x)
                u = group.add[y][difference]
                row[d * x + y] ^= 1 << (d * v + u)
                column[d * x + u] ^= 1 << (d * v + y)
                j = group.sub(x, u)
                assert j == group.sub(v, y)
                blocks.add(j)
    return tuple(row), tuple(column), blocks


def column_blocks(group, kernels):
    d = group.order
    answer = []
    for j in range(d):
        block = []
        for x in range(d):
            value = 0
            for v in range(d):
                w = group.sub(group.add[x][v], j)
                difference = group.sub(x, v)
                if (kernels[w] >> difference) & 1:
                    value ^= 1 << v
            block.append(value)
        answer.append(tuple(block))
    return answer


def tensor_kernel(first, second, f, g):
    combined = AbelianGroup(first.moduli + second.moduli)
    d = combined.order
    answer = [0] * d
    for wi, w in enumerate(combined.elements):
        a = first.index[w[:len(first.moduli)]]
        b = second.index[w[len(first.moduli):]]
        for zi, z in enumerate(combined.elements):
            c = first.index[z[:len(first.moduli)]]
            e = second.index[z[len(first.moduli):]]
            answer[wi] |= (((f[a] >> c) & 1) & ((g[b] >> e) & 1)) << zi
    return combined, answer


def smith_valuations(coefficients, cap):
    size = len(coefficients[0])
    degree = len(coefficients) - 1
    previous_rank = previous_count = 0
    multiplicities = []
    for order in range(cap + 1):
        rows = []
        for block in range(order + 1):
            for r in range(size):
                value = 0
                for j in range(max(0, block - degree), block + 1):
                    value ^= coefficients[block - j][r] << (size * j)
                rows.append(value)
        current_rank = binary_rank(rows)
        count = current_rank - previous_rank
        assert previous_count <= count <= size
        multiplicities.append(count - previous_count)
        if count == size:
            return multiplicities
        previous_rank, previous_count = current_rank, count
    raise AssertionError('Insufficient truncation depth')


def cyclic_polynomial(d):
    coefficients = [[0] * d for _ in range(d)]
    for w in range(d):
        k = w
        while True:
            j = k
            while True:
                coefficients[w][d - 1 - w + j] ^= 1 << k
                if not j:
                    break
                j = (j - 1) & k
            if not k:
                break
            k = (k - 1) & w
    return coefficients


def verify():
    rng = random.Random(53729)
    cases = block_cases = decomposition_cases = tensor_cases = 0
    for moduli in ((2,), (3,), (4,), (2, 2), (8,), (4, 2), (3, 3)):
        group = AbelianGroup(moduli)
        d = group.order
        m = len(group.double_image)
        cross = [1] * d
        row, col, _ = matrices(group, cross)
        assert row == col == tuple(1 << i for i in range(d * d))
        if d <= 3:
            samples = range(1, 1 << (d * d))
        else:
            samples = [1 << i for i in range(d * d)] + [rng.getrandbits(d * d) or 1 for _ in range(128)]
        for packed in samples:
            kernels = [(packed >> (d * w)) & ((1 << d) - 1) for w in range(d)]
            row, col, nonzero = matrices(group, kernels)
            blocks = column_blocks(group, kernels)
            assert binary_rank(col) == sum(binary_rank(b) for b in blocks)
            assert len(nonzero) >= m
            assert binary_rank(col) >= len(nonzero) >= m
            cases += 1
            block_cases += d
        for _ in range(20):
            k = rng.randrange(2, 6)
            terms = [[rng.getrandbits(d) for _ in range(d)] for _ in range(k - 1)]
            final = cross[:]
            for term in terms:
                final = [a ^ b for a, b in zip(final, term)]
            terms.append(final)
            pairs = []
            for term in terms:
                row, col, _ = matrices(group, term)
                pairs.append((binary_rank(row), binary_rank(col)))
            assert max(r * c for r, c in pairs) * k >= d * d * m
            decomposition_cases += 1

    for first_moduli, second_moduli in (((2,), (4,)), ((3,), (2,)), ((4,), (4,))):
        first, second = AbelianGroup(first_moduli), AbelianGroup(second_moduli)
        for _ in range(8):
            f = [rng.getrandbits(first.order) for _ in range(first.order)]
            g = [rng.getrandbits(second.order) for _ in range(second.order)]
            combined, h = tensor_kernel(first, second, f, g)
            a, b, _ = matrices(first, f)
            c, e, _ = matrices(second, g)
            r, s, _ = matrices(combined, h)
            assert binary_rank(r) == binary_rank(a) * binary_rank(c)
            assert binary_rank(s) == binary_rank(b) * binary_rank(e)
            assert len(combined.double_image) == len(first.double_image) * len(second.double_image)
            tensor_cases += 1

    power_certificates = []
    for d in (4, 8, 16, 32):
        rational = (Fraction(2 * d * d, d + 2) ** (d + 2)
                    * Fraction(d * d * (d * d + 2 * d - 2), d - 1) ** (2 * d - 2)
                    * (2 * d * (2 * d - 1) * (d * d - 1)) ** d)
        assert rational.numerator ** 7 > rational.denominator ** 7 * d ** (80 * d)
        weights = (Fraction(d + 2, 4 * d), Fraction(2 * d - 2, 4 * d), Fraction(1, 4))
        t = Fraction(2 * d - 1, 2 * d)
        assert sum(weights) == 1
        assert sum(w * (i - t) for i, w in enumerate(weights)) == 0
        power_certificates.append({'labels': d, 'strict_power_test': 'pass'})
    assert 2 * 5 ** 2 > 7 ** 2
    assert 2 * 17 ** 2 > 24 ** 2

    d = 8
    group = AbelianGroup((d,))
    F = cyclic_polynomial(d)
    row_vals = [0] * (2 * d - 1)
    col_vals = [0] * (2 * d - 1)
    row_ranks = [0] * d
    col_ranks = [0] * d
    for w in range(d):
        blocks = [tuple(sum(((F[w][i] >> group.sub(x, v)) & 1) << v for v in range(d))
                        for x in range(d)) for i in range(d)]
        vals = smith_valuations(blocks, 2 * d)
        for i, count in enumerate(vals):
            row_vals[i] += count
        for i, block in enumerate(blocks):
            row_ranks[i] += binary_rank(block)
    for j in range(d):
        blocks = [column_blocks(group, [F[w][i] for w in range(d)])[j] for i in range(d)]
        vals = smith_valuations(blocks[::-1], 2 * d)
        for i, count in enumerate(vals):
            col_vals[i] += count
        for i, block in enumerate(blocks):
            col_ranks[i] += binary_rank(block)
    assert row_vals == list(range(1, 9)) + list(range(7, 0, -1))
    assert col_vals == [4] * 7 + [8] + [4] * 7
    full = []
    for i in range(d):
        r, c, _ = matrices(group, [F[w][i] for w in range(d)])
        assert binary_rank(r) == row_ranks[i]
        assert binary_rank(c) == col_ranks[i]
        full.append(r)
    assert all(__import__('functools').reduce(int.__xor__, (r[i] for r in full)) == 1 << i
               for i in range(d * d))

    return {
        'status': 'pass',
        'scope': 'full conserved-convolution label domains; no claim about arbitrary label restrictions',
        'nonzero_kernel_cases': cases,
        'column_blocks_checked': block_cases,
        'exact_crossing_decompositions': decomposition_cases,
        'tensor_product_closure_cases': tensor_cases,
        'rational_am_gm_certificates': power_certificates,
        'cyclic_eight': {'row_valuations': row_vals, 'reverse_column_valuations': col_vals,
                         'row_coefficient_ranks': row_ranks, 'column_coefficient_ranks': col_ranks},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output')
    args = parser.parse_args()
    encoded = json.dumps(verify(), indent=2, sort_keys=True) + '\n'
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end='')
