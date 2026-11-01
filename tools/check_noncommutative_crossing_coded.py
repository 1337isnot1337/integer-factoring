#!/usr/bin/env python3
"""Exact constant-cut and filtration checks for coded central quotients.

The bounds concern individual cuts. They do not establish a better complete
crossing decomposition or an improved factoring exponent.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path

from check_noncommutative_crossing_orbit import binary_rank, group_word_polynomial


def vector_rank(vectors):
    return binary_rank(list(vectors))


class CodedGroup:
    def __init__(self, columns, center):
        self.columns = tuple(columns)
        self.m, self.t = len(columns), center
        self.q = 1 << self.m
        self.n = 2 * self.m
        self.v = 1 << self.n
        self.width = self.n + self.t
        self.d = 1 << self.width
        assert all(0 < c < 1 << center for c in columns)
        assert vector_rank(columns) == center
        self.linear = tuple(self.apply(mask) for mask in range(self.q))
        self.generators = tuple(1 << i for i in range(self.width))
        self.coordinates = [0] * self.d
        for bits in range(self.d):
            a, b, c = self.decode(bits)
            self.coordinates[bits ^ (self.linear[a & b] << self.n)] = bits

    def apply(self, mask):
        value = 0
        for j, column in enumerate(self.columns):
            if mask >> j & 1:
                value ^= column
        return value

    def decode(self, value):
        return value & (self.q - 1), (value >> self.m) & (self.q - 1), value >> self.n

    def mul(self, left, right):
        a, b, c = self.decode(left)
        A, B, C = self.decode(right)
        return left ^ right ^ (self.linear[a & B] << self.n)


def independent_monomials(group):
    """The parallel extension of the binary column matroid, with two copies."""
    result = [0]
    for size in range(1, min(group.m, group.t) + 1):
        for selected in itertools.combinations(range(group.m), size):
            if vector_rank(group.columns[i] for i in selected) != size:
                continue
            for choices in range(1 << size):
                result.append(sum(1 << (i + group.m * ((choices >> j) & 1))
                                  for j, i in enumerate(selected)))
    return tuple(sorted(result))


def multiply_boolean(left, right):
    result = set()
    for a in left:
        for b in right:
            value = a | b
            if value in result:
                result.remove(value)
            else:
                result.add(value)
    return result


def insert_row(pivots, row):
    while row:
        lead = row.bit_length() - 1
        if lead not in pivots:
            pivots[lead] = row
            return
        row ^= pivots[lead]


def constant_rank_polynomial(group):
    monomials = independent_monomials(group)
    positions = {monomial: i for i, monomial in enumerate(monomials)}
    bases = [{} for _ in range(1 << group.t)]
    examined = 0
    for h in sorted(range(group.v), key=int.bit_count):
        a, b, _ = group.decode(h)
        z = group.linear[a & b]
        if len(bases[z]) == len(monomials):
            continue
        linear_forms = []
        for j in range(group.t):
            form = set()
            for i, column in enumerate(group.columns):
                if column >> j & 1:
                    if b >> i & 1:
                        form.add(1 << i)
                    if a >> i & 1:
                        form.add(1 << (group.m + i))
            linear_forms.append(form)
        polynomials = [{0}]
        for form in linear_forms:
            polynomials += [multiply_boolean(poly, form) for poly in polynomials]
        for poly in polynomials:
            assert all(monomial in positions for monomial in poly)
            row = sum(1 << positions[monomial] for monomial in poly)
            insert_row(bases[z], row)
        examined += 1
        if all(len(basis) == len(monomials) for basis in bases):
            break
    ranks = [len(basis) for basis in bases]
    return {"constant_column_rank": sum(ranks), "fiber_ranks": ranks,
            "independent_monomial_count": len(monomials),
            "matroid_upper_bound": (1 << group.t) * len(monomials),
            "degree_upper_bound": (1 << group.t) * sum(
                math.comb(group.n, j) for j in range(group.t + 1)),
            "horizontal_rows_examined": examined}


def constant_rank_direct(group):
    """Independent function-evaluation matrices, without polynomial expansion."""
    fibers = [[] for _ in range(1 << group.t)]
    for h in range(group.v):
        a, b, _ = group.decode(h)
        rows = [0] * (1 << group.t)
        for H in range(group.v):
            A, B, _ = group.decode(H)
            delta = group.linear[(a & B) ^ (A & b)]
            rows[delta] |= 1 << H
        fibers[group.linear[a & b]].extend(rows)
    return [binary_rank(rows) for rows in fibers]


def hilbert_coefficients(horizontal, center):
    result = [0] * (horizontal + 2 * center + 1)
    for a in range(horizontal + 1):
        for b in range(center + 1):
            result[a + 2 * b] += math.comb(horizontal, a) * math.comb(center, b)
    return result


def filtration_rank_bound(horizontal, center, shift):
    coefficients = hilbert_coefficients(horizontal, center)
    dimension = 1 << (horizontal + center)
    best = dimension
    for cutoff in range(len(coefficients) + 1):
        best = min(best, sum(coefficients[:cutoff]) + sum(coefficients[cutoff + shift:]))
    return best


def row_rank_upper_bounds(group):
    states = Counter()
    for a, b, c in itertools.product(range(group.q), range(group.q), range(1 << group.t)):
        vectors = [1 << j for j in range(group.t) if c >> j & 1]
        vectors += [column for i, column in enumerate(group.columns) if (a & b) >> i & 1]
        k, r = a.bit_count() + b.bit_count(), vector_rank(vectors)
        states[k, r, k + 2 * c.bit_count()] += 1
    maximum = group.n + 2 * group.t
    bounds = []
    for degree in range(maximum + 1):
        bound = sum(count * (group.d >> (k + r)) * filtration_rank_bound(k, r, degree)
                    for (k, r, weight), count in states.items() if weight >= degree)
        bounds.append(bound)
    assert bounds[0] == group.d * group.d and bounds[-1] == 1
    return bounds


def direct_row_ranks(group):
    maximum = group.n + 2 * group.t
    weights = (1,) * group.n + (2,) * group.t
    ranks = [0] * (maximum + 1)
    for S in range(group.d):
        kernel = group_word_polynomial(group, S, tuple(range(group.width)), weights)
        matrices = [[0] * group.d for _ in ranks]
        for y in range(group.d):
            for a, polynomial in kernel.items():
                u = group.mul(y, a)
                while polynomial:
                    bit = polynomial & -polynomial
                    matrices[bit.bit_length() - 1][y] ^= 1 << u
                    polynomial ^= bit
        for degree, rows in enumerate(matrices):
            ranks[degree] += binary_rank(rows)
    return ranks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--larger", action="store_true")
    args = parser.parse_args()
    cases = [(1, (1,)), (1, (1, 1)), (1, (1, 1, 1)),
             (2, (1, 2)), (2, (1, 2, 3)), (2, (1, 1, 2, 3)),
             (2, (1, 2, 1, 2)), (3, (1, 2, 4, 3))]
    if args.larger:
        cases += [(2, (1, 2, 3, 1, 2, 3, 1, 2)),
                  (3, (1, 2, 4, 3, 6, 7))]
    records = []
    for center, columns in cases:
        group = CodedGroup(columns, center)
        constant = constant_rank_polynomial(group)
        if group.m <= 4:
            assert constant["fiber_ranks"] == constant_rank_direct(group)
        record = {"horizontal_pairs": group.m, "center_bits": center,
                  "columns": list(columns), "label_bits": group.width,
                  **constant}
        if group.width <= 7:
            record["row_rank_upper_bounds_center_weight_two"] = row_rank_upper_bounds(group)
        if group.width <= 6:
            exact = direct_row_ranks(group)
            assert all(a <= b for a, b in zip(exact, record["row_rank_upper_bounds_center_weight_two"]))
            record["exact_row_ranks_center_weight_two"] = exact
        records.append(record)
    result = {"status": "exact constant-cut identities and proven row upper bounds; no full profile improvement",
              "records": records,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
