"""Exact rational search for alternative bases of the crossing coefficients.

For t=2,4,6 this enumerates every one-dimensional intersection of t row/column
eigenvalue-zero hyperplanes. Such intersections suffice to minimize the rank
product subject to a nonzero central coefficient. The general obstruction is
proved in research/coefficient-basis.md; this is a finite independent audit.
"""

from __future__ import annotations

from functools import reduce
from itertools import combinations
from math import comb, gcd
import json


def eigenvalue_rows(t: int) -> list[tuple[int, ...]]:
    # Actual eigenvalues have a common nonzero factor 2^t, omitted here.
    return [
        tuple(0 if i < j else comb(t - j, i - j) * 4 ** (t - i) * (-3) ** (i - j)
              for i in range(t + 1))
        for j in range(t + 1)
    ]


def determinant(matrix: list[list[int]]) -> int:
    """Bareiss elimination with exact divisions and row pivoting."""
    n = len(matrix)
    if n == 0:
        return 1
    a = [row[:] for row in matrix]
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot = next((r for r in range(k, n) if a[r][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign *= -1
        value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = value * a[i][j] - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = value
    return sign * a[-1][-1]


def kernel_line(rows: tuple[tuple[int, ...], ...]) -> tuple[int, ...] | None:
    t = len(rows)
    coordinates = []
    for omitted in range(t + 1):
        minor = [list(row[:omitted] + row[omitted + 1:]) for row in rows]
        coordinates.append((-1) ** omitted * determinant(minor))
    if not any(coordinates):
        return None
    divisor = reduce(gcd, (abs(x) for x in coordinates))
    if next(x for x in coordinates if x) < 0:
        divisor *= -1
    result = tuple(x // divisor for x in coordinates)
    assert all(sum(x * y for x, y in zip(row, result)) == 0 for row in rows)
    return result


def rank_pair(t: int, coefficients: tuple[int, ...]) -> tuple[int, int]:
    rows = eigenvalue_rows(t)
    weights = [comb(t, j) * 8 ** (t - j) * 56 ** j for j in range(t + 1)]
    row_rank = sum(weight for functional, weight in zip(rows, weights)
                   if sum(a * b for a, b in zip(functional, coefficients)))
    column_rank = sum(weight for functional, weight in zip(rows, weights)
                      if sum(a * b for a, b in zip(reversed(functional), coefficients)))
    return row_rank, column_rank


def search(t: int) -> dict:
    if t < 2 or t % 2:
        raise ValueError("this exact minimax audit uses even t >= 2")
    m = t // 2
    rows = eigenvalue_rows(t)
    hyperplanes = tuple(rows + [tuple(reversed(row)) for row in rows])
    candidates = set()
    for selected in combinations(hyperplanes, t):
        line = kernel_line(selected)
        if line is not None and line[m]:
            candidates.add(line)
    central = tuple(int(i == m) for i in range(t + 1))
    expected_pair = rank_pair(t, central)
    expected = expected_pair[0] * expected_pair[1]
    best = min(rank_pair(t, c)[0] * rank_pair(t, c)[1] for c in candidates)
    minimizers = sorted(c for c in candidates if rank_pair(t, c)[0] * rank_pair(t, c)[1] == best)
    assert best == expected
    assert minimizers == [central]
    return {
        "t": t,
        "hyperplane_subsets": comb(2 * t + 2, t),
        "distinct_central_nonzero_lines": len(candidates),
        "minimum_rank_product": best,
        "central_row_rank": expected_pair[0],
        "unique_minimizer_coefficients": list(central),
    }


if __name__ == "__main__":
    print(json.dumps({
        "method": "exact integer minors for rational kernel lines",
        "scope": "linear span of the symmetric coefficient matrices A_i",
        "results": [search(t) for t in (2, 4, 6)],
    }, indent=2))
