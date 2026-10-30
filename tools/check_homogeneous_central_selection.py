#!/usr/bin/env python3
"""Homogeneous dual selection and complete raw-parameter Pfaffian profiles.

Exact finite evidence for homogeneous-central-selection.md.  The all-width
existence argument uses rank subadditivity, not the finite Gram-rank guess.
"""

import argparse
from functools import lru_cache
from itertools import product
import json

from check_angle_pfaffian import pfaffian_polynomial
from check_homogeneous_central_symbols import binomial_two, boolean_zeta_two
from check_modular_crossing import (
    column_bits, determinant_operator_two, inverse_two, rank_two,
    regular_partitions,
)
from check_partition_channels import operator


def regular_counts(bound: int) -> list[int]:
    counts = [0] * (bound + 1)
    counts[0] = 1
    for part in range(1, bound + 1, 2):
        for weight in range(bound, part - 1, -1):
            counts[weight] += counts[weight - part]
    return counts


def canonical_weight(bits: int) -> int:
    # The proof gives this elementary finite search bound.
    bound = 8 * (2 * bits + 16) ** 2
    counts = regular_counts(bound)
    return next(h for h in range(1, bound + 1)
                if counts[h] >= (h + 1) * (1 << bits))


def independent_rows(rows, wanted: int):
    """Select actual rows and pivot columns, preserving a nonsingular minor."""
    pivots = {}
    selected = []
    for index, original in enumerate(rows):
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                selected.append(index)
                break
        if len(selected) == wanted:
            return selected, tuple(pivots)
    return None


def choose_duals(weight: int, wanted: int):
    states = regular_partitions(2, weight)
    # Prefer fewer Lucas-supported angle indices, then larger operator weight.
    candidates = sorted(range(weight + 1),
                        key=lambda j: ((weight + j).bit_count(), -j))
    for j in candidates:
        parameters = regular_partitions(2, j)
        if len(parameters) < wanted or (weight - j) % 2:
            continue
        rows = [sum(determinant_operator_two(a, (), b) << k
                    for k, b in enumerate(states)) for a in parameters]
        chosen = independent_rows(rows, wanted)
        if chosen is None:
            continue
        row_indices, columns = chosen
        minor = [sum(((rows[i] >> col) & 1) << k
                     for k, col in enumerate(columns)) for i in row_indices]
        inverse = inverse_two(minor, wanted)
        basis = tuple(parameters[i] for i in row_indices)
        labels = tuple(states[i] for i in columns)
        duals = tuple(tuple(basis[k] for k in range(wanted) if row >> k & 1)
                      for row in inverse)
        for x, terms in enumerate(duals):
            for v, label in enumerate(labels):
                assert sum(determinant_operator_two(a, (), label)
                           for a in terms) % 2 == int(x == v)
        return j, labels, duals
    raise ValueError("no homogeneous vacuum minor of the requested size")


@lru_cache(None)
def raw_homogeneous_profile(a, u, v, y) -> int:
    """Pfaffian = Z^(h+j) sum_s Q_s(a,t^v)_{u,y} Z^s."""
    j, h = sum(a), sum(v)
    assert sum(u) == sum(y) == h and j <= h
    vertices = ([('U', t) for t in u] + [('A', n) for n in a]
                + [('B', m) for m in v] + [('Y', t) for t in y])
    matrix = [[0] * len(vertices) for _ in vertices]
    for row, (first, n) in enumerate(vertices):
        for col in range(row + 1, len(vertices)):
            second, m = vertices[col]
            value = 0
            if (first, second) == ('U', 'A'):
                if (m - 1) & ~(n - 1) == 0:
                    value = binomial_two(m) << n
            elif (first, second) == ('A', 'Y'):
                if n & m == 1:
                    value = binomial_two(n)
            elif (first, second) == ('A', 'B'):
                if n & m == 1:
                    value = 1 << n
            elif (first, second) == ('U', 'B'):
                if (m - 1) & ~(n - 1) == 0:
                    value = binomial_two(m)
            elif (first, second) == ('B', 'Y'):
                if n & m == 1:
                    value = binomial_two(n) << m
            elif (first, second) == ('U', 'Y') and n == m:
                value = 1 << m
            matrix[row][col] = matrix[col][row] = value
    result = pfaffian_polynomial(matrix)
    assert result.bit_length() <= 3 * h + j + 1
    return result


def reference_raw_profile(a, u, v, y) -> int:
    """Independent literal operator composition in the original order."""
    h, j = sum(v), sum(a)
    support = tuple(i for i in range(h + j + 1) if i & ~(h + j) == 0)
    result = 0
    for shift in range(-h - j, 2 * h + 1):
        value = 0
        for index in support:
            energy = 2 * h - shift - index
            for intermediate in regular_partitions(2, energy):
                value ^= ((operator(a, u, intermediate) % 2)
                          & (operator(v, intermediate, y) % 2))
        result |= value << (shift + h + j)
    return result


def complete_family_record(weight: int, wanted: int) -> dict:
    j, labels, duals = choose_duals(weight, wanted)
    matrices = [[0] * (wanted * wanted) for _ in range(3 * weight + j + 1)]
    for x, y, u, v in product(range(wanted), repeat=4):
        polynomial = 0
        for a in duals[x]:
            polynomial ^= raw_homogeneous_profile(a, labels[u], labels[v], labels[y])
        assert polynomial.bit_count() % 2 == int(x == v and y == u)
        while polynomial:
            bit = polynomial & -polynomial
            degree = bit.bit_length() - 1
            matrices[degree][x * wanted + y] ^= 1 << (v * wanted + u)
            polynomial ^= bit
    profile = []
    for degree, rows in enumerate(matrices):
        if any(rows):
            profile.append([degree - weight - j, rank_two(rows),
                            rank_two(column_bits(rows, wanted))])
    return {
        "state_weight": weight, "operator_weight": j, "labels": wanted,
        "angle_indices": [i for i in range(weight + j + 1)
                          if i & ~(weight + j) == 0],
        "state_partitions": [list(x) for x in labels],
        "homogeneous_duals": [[list(a) for a in terms] for terms in duals],
        "complete_crossing_entries": wanted ** 4,
        "nonzero_row_column_profile": profile,
        "maximum_rank_product": max(r * c for _, r, c in profile),
    }


def verify(canonical_bits: int = 1) -> dict:
    identity_entries = 0
    rank_budget_cases = 0
    for h in range(1, 17):
        labels = regular_partitions(2, h)
        if not labels:
            continue
        blocks = {}
        rank_sum = 0
        for j in range(h + 1):
            basis = regular_partitions(2, j)
            rows = [sum(determinant_operator_two(a, (), b) << k
                        for k, b in enumerate(labels)) for a in basis]
            blocks.update(zip(basis, rows))
            rank_sum += rank_two(rows)
        assert rank_sum >= len(labels)
        rank_budget_cases += 1
        for i, label in enumerate(labels):
            value = 0
            for monomial in boolean_zeta_two(label):
                value ^= blocks[monomial]
            assert value == 1 << i
            identity_entries += len(labels)

    raw_profile_checks = 0
    for h in range(1, 9):
        states = regular_partitions(2, h)
        for j in range(1, h + 1):
            for a in regular_partitions(2, j):
                for u, v, y in product(states, repeat=3):
                    assert raw_homogeneous_profile(a, u, v, y) == reference_raw_profile(a, u, v, y)
                    raw_profile_checks += 1

    records = [complete_family_record(12, 2), complete_family_record(17, 4)]
    if canonical_bits:
        h = canonical_weight(canonical_bits)
        record = complete_family_record(h, 1 << canonical_bits)
        record["canonical_label_bits"] = canonical_bits
        record["regular_state_dimension"] = len(regular_partitions(2, h))
        record["pigeonhole_threshold"] = (h + 1) * (1 << canonical_bits)
        assert record["regular_state_dimension"] >= record["pigeonhole_threshold"]
        records.append(record)
    return {
        "status": "pass",
        "scope": "uniform homogeneous dual selection and exact finite crossings; no factoring exponent above one half",
        "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "checks": {"source_duality_matrix_entries": identity_entries,
                   "homogeneous_rank_budget_cases": rank_budget_cases,
                   "raw_pfaffian_profiles_against_original_composition": raw_profile_checks},
        "complete_homogeneous_families": records,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--canonical-bits', type=int, default=1,
                        help='verify the canonical construction at this width; zero skips it')
    parser.add_argument('--output')
    args = parser.parse_args()
    report = verify(args.canonical_bits)
    encoded = json.dumps(report, sort_keys=True, indent=2) + '\n'
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end='')
