#!/usr/bin/env python3
"""Exact greedy Pascal-complement Smith profiles for arbitrary finite modes."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from check_block_smith_crossings import binary_rank, smith_profile, toeplitz_rank
from check_schur_gaussian_reduction import polynomial_multiply
from check_schur_one_species_interval import multiply, power


ROOT = Path(__file__).resolve().parent.parent


def greedy_complement(indices):
    """Greedy Pascal columns, skipping whole dependent binary intervals.

    The traversal depends on the bit length of the largest index rather
    than scanning every integer below it.
    """
    indices = tuple(sorted(indices))
    if not indices or len(set(indices)) != len(indices) or indices[0] < 0:
        raise ValueError("a nonempty finite set of nonnegative indices is required")
    excluded = set(indices)
    basis, chosen = {}, []
    def reduce(column):
        while column:
            pivot = column.bit_length() - 1
            if pivot in basis:
                column ^= basis[pivot]
            else:
                break
        return column

    def visit(prefix, free_bits):
        if len(chosen) == len(indices):
            return
        groups = {}
        low_mask = (1 << free_bits) - 1
        for i, a in enumerate(indices):
            if (a >> free_bits) & ~(prefix >> free_bits) == 0:
                low = a & low_mask
                groups[low] = groups.get(low, 0) | (1 << i)
        if all(reduce(column) == 0 for column in groups.values()):
            return
        if free_bits:
            visit(prefix, free_bits - 1)
            visit(prefix + (1 << (free_bits - 1)), free_bits - 1)
        elif prefix not in excluded:
            column = reduce(sum(int(a & ~prefix == 0) << i for i, a in enumerate(indices)))
            if column:
                basis[column.bit_length() - 1] = column
                chosen.append(prefix)

    visit(0, indices[-1].bit_length() + 1)
    assert len(chosen) == len(indices), "one repeated binary period must span every row"
    return tuple(chosen)


def greedy_complement_linear_scan(indices):
    """Independent small-input implementation used only by the checker."""
    basis, chosen = {}, []
    for c in range(1 << (max(indices).bit_length() + 1)):
        if c in indices:
            continue
        column = sum(int(a & ~c == 0) << i for i, a in enumerate(indices))
        while column:
            pivot = column.bit_length() - 1
            if pivot in basis:
                column ^= basis[pivot]
            else:
                basis[pivot] = column
                chosen.append(c)
                break
        if len(chosen) == len(indices):
            return tuple(chosen)
    raise AssertionError("the independent scan did not span")


def normalized_matrix(indices):
    indices = tuple(indices)
    n = len(indices)
    result = [[0] * (2 * n) for _ in range(2 * n)]
    for i, a in enumerate(indices):
        for j, b in enumerate(indices):
            result[i][j] = (power(3, a) << b) if a & ~b == 0 else 0
            result[i][n + j] = (1 << a) if b & ~a == 0 else 0
            result[n + i][j] = (1 << a) if i == j else 0
            result[n + i][n + j] = 3 << (a | b)
    return result


def check_set(indices, independent=False):
    indices = tuple(indices)
    n = len(indices)
    chosen = greedy_complement(indices)
    assert chosen == greedy_complement_linear_scan(indices)
    expected = tuple(sorted((*indices, *chosen)))
    assert len(set(expected)) == 2 * n
    matrix = normalized_matrix(indices)
    degree = max(value.bit_length() - 1 for row in matrix for value in row)
    actual, precision = smith_profile(matrix, degree)
    assert actual == expected

    # Independently verify the formal complement expansion through enough
    # precision to see every predicted invariant factor.
    jets = max(expected) + 1
    mask = (1 << jets) - 1
    p = [[int(b & ~a == 0) for b in indices] for a in indices]
    q = [[power(3, a) if i == j else 0 for j in range(n)] for i, a in enumerate(indices)]
    d = [[(1 << a) if i == j else 0 for j in range(n)] for i, a in enumerate(indices)]
    o = [[1 << (a | b) for b in indices] for a in indices]
    dp = multiply(d, p)
    other = multiply(multiply(q, list(map(list, zip(*p)))), o)
    schur = [[(dp[i][j] ^ polynomial_multiply(3, other[i][j])) & mask
              for j in range(n)] for i in range(n)]
    expanded = [[0] * n for _ in range(n)]
    excluded = set(indices)
    for c in range(jets):
        if c in excluded:
            continue
        for i, a in enumerate(indices):
            if a & ~c:
                continue
            for j, b in enumerate(indices):
                expanded[i][j] ^= (power(3, a + 1) << (c | b)) & mask
    assert expanded == schur

    # The chosen columns are independent, and each rejected earlier
    # column was dependent on the preceding chosen columns.
    selected_columns = [sum(int(a & ~c == 0) << i for i, a in enumerate(indices)) for c in chosen]
    assert binary_rank(selected_columns) == n
    toeplitz_checks = 0
    if independent:
        for count in range(1, max(expected) + 2):
            assert toeplitz_rank(matrix, count) == sum(max(0, count - order) for order in expected)
            toeplitz_checks += 1
    return {"indices": indices, "greedy_complement": chosen,
            "normalized_t_orders": actual, "exact_smith_precision": precision,
            "formal_series_checked_through": jets,
            "independent_toeplitz_ranks": toeplitz_checks}


def check_intervals():
    records = []
    for n in range(1, 65):
        indices = tuple(range(n))
        chosen = greedy_complement(indices)
        assert chosen == tuple(range(n, 2 * n))
        consecutive = [sum(int(a & ~(n + j) == 0) << a for a in range(n)) for j in range(n)]
        assert binary_rank(consecutive) == n
        records.append({"modes": n, "complement_start": n, "complement_end": 2 * n - 1})
    return records


def check_large_binary_indices():
    records = []
    for bits in (20, 64, 128, 256):
        high = 1 << bits
        cases = [((high - 1,), (2 * high - 1,)),
                 ((0, high), (1, high + 1)),
                 (tuple(high + a for a in range(13)), tuple(high + a for a in range(13, 26)))]
        for indices, expected in cases:
            actual = greedy_complement(indices)
            assert actual == expected
            records.append({"index_bits": bits, "indices": indices, "complement": actual})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cases = [tuple(i for i in range(8) if selector >> i & 1) for selector in range(1, 256)]
    cases += [(0, 3, 5, 8, 12), (1, 2, 4, 8), (7, 11, 13, 14), (0, 1, 2, 4, 8, 16)]
    independent = {(0,), (1,), (2,), (0, 2), (1, 2), (0, 1, 3), (1, 3, 5)}
    result = {
        "status": "exact arbitrary-mode Pascal-complement Smith formula; no stronger factoring exponent",
        "scope": "all finite nonempty sets of nonnegative mode indices; displayed normalized matrix",
        "sets": [check_set(indices, indices in independent) for indices in cases],
        "initial_intervals": check_intervals(),
        "large_index_binary_traversal_controls": check_large_binary_indices(),
        "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in (
            "tools/check_schur_pascal_complement.py", "tools/check_schur_one_species_interval.py",
            "tools/check_block_smith_crossings.py")},
    }
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
