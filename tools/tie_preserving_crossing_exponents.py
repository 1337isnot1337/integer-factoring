#!/usr/bin/env python3
"""Exact linear constraints preserving tied leading two-by-two minors.

This constructs search directions. Preserving these ties does not imply
preservation of higher minors, weighted profiles, or a factoring rate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import lcm
from pathlib import Path

from relabelled_gaussian_sums import gaussian_tensor


def rref(rows, columns):
    basis = {}
    for source in rows:
        row = [Fraction(value) for value in source]
        for pivot, previous in sorted(basis.items()):
            if row[pivot]:
                factor = row[pivot]
                row = [a - factor * b for a, b in zip(row, previous)]
        pivot = next((i for i, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        factor = row[pivot]
        row = [value / factor for value in row]
        for old, previous in list(basis.items()):
            if previous[pivot]:
                factor = previous[pivot]
                basis[old] = [a - factor * b for a, b in zip(previous, row)]
        basis[pivot] = row
    assert all(len(row) == columns for row in basis.values())
    return basis


def nullspace(rows, columns):
    basis = rref(rows, columns)
    vectors = []
    for free in range(columns):
        if free in basis:
            continue
        vector = [Fraction(0)] * columns
        vector[free] = Fraction(1)
        for pivot, row in basis.items():
            vector[pivot] = -row[free]
        denominator = lcm(*(value.denominator for value in vector))
        vectors.append([int(value * denominator) for value in vector])
    return basis, vectors


def construct(count):
    q = 1 << count
    if count & (count - 1):
        raise ValueError('The displayed full-interval control uses a power-of-two mode count.')
    tensor = gaussian_tensor(tuple(range(count)), count)
    support = [mask for mask, value in enumerate(tensor) if value]
    position = {mask: i for i, mask in enumerate(support)}
    constraints = []
    counts = []
    for infinity in (False, True):
        cells = {}
        for mask in support:
            x, y, u, v = ((mask >> (count * leg)) & (q - 1) for leg in range(4))
            row, column = (q * x + u, q * y + v) if infinity else (q * x + y, q * u + v)
            value = tensor[mask]
            order = -(value.bit_length() - 1) if infinity else (value & -value).bit_length() - 1
            cells[row, column] = position[mask], order
        before = len(constraints)
        for i, k in combinations(range(q * q), 2):
            for j, ell in combinations(range(q * q), 2):
                keys = [(i, j), (k, ell), (i, ell), (k, j)]
                if not all(key in cells for key in keys):
                    continue
                entries = [cells[key] for key in keys]
                if entries[0][1] + entries[1][1] != entries[2][1] + entries[3][1]:
                    continue
                row = [0] * len(support)
                for (index, _), sign in zip(entries, (1, 1, -1, -1)):
                    row[index] += sign
                constraints.append(row)
        counts.append(len(constraints) - before)
    reduced, vectors = nullspace(constraints, len(support))
    leg_vectors = [[1] * len(support)]
    leg_vectors.extend([[int(((mask >> (count * leg)) & (q - 1)) == label) for mask in support]
                        for leg in range(4) for label in range(1, q)])
    assert all(sum(a * b for a, b in zip(row, vector)) == 0
               for row in constraints for vector in leg_vectors)
    span = list(leg_vectors)
    rank = len(rref(span, len(support)))
    leg_rank = rank
    extra = []
    for vector in vectors:
        following = len(rref(span + [vector], len(support)))
        if following > rank:
            span.append(vector)
            extra.append(vector)
            rank = following
    assert rank == len(vectors)
    assert all(sum(a * b for a, b in zip(row, vector)) == 0
               for row in constraints for vector in extra)
    return {'modes': count, 'alphabet': q, 'support': support,
            'constraints_by_place': counts, 'constraint_rank': len(reduced),
            'nullity': len(vectors), 'physical_leg_gauge_rank': leg_rank,
            'additional_integer_directions': extra,
            'scope': 'Exact leading two-by-two tie constraints only; higher minors are not certified.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', type=int, default=2)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = construct(args.modes)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in {'support', 'additional_integer_directions'}}, indent=2))


if __name__ == '__main__':
    main()
