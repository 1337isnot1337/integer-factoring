"""Exact finite checks for sparse rectangular crossing layouts.

Enumerate every 0/1 support of a square grid of side at most four.  Retain
regular supports and check rectangle counts, the exact channel-volume
formula, and the extremal union-of-bicliques examples.  No matrix-rank or
factoring-runtime claim is inferred from these finite checks.
"""
from fractions import Fraction
from itertools import combinations
from math import comb
import json


def inspect_grid(rows, k):
    degree = [row.bit_count() for row in rows]
    column_degree = [sum((row >> j) & 1 for row in rows) for j in range(k)]
    if len(set(degree + column_degree)) != 1:
        return None
    d = degree[0]
    if not d:
        return None
    vertices = [(a, b) for a in range(k) for b in range(k)
                if rows[a] >> b & 1]
    n = len(vertices)
    rect = sum(comb((rows[a] & rows[b]).bit_count(), 2)
               for a, b in combinations(range(k), 2))
    assert 4 * rect <= n * (d - 1) ** 2
    row_side = sum(comb(x, 2) for x in degree)
    col_side = sum(comb(x, 2) for x in column_degree)
    unmatched = 0
    covered = 0
    for (a, b), (c, e) in combinations(vertices, 2):
        if a == c or b == e:
            continue
        if rows[a] >> e & 1 and rows[c] >> b & 1:
            covered += 1
        else:
            unmatched += 1
    assert covered == 2 * rect
    assert comb(n, 2) == row_side + col_side + unmatched + covered
    cases = 0
    # c is a common declared row/column width; lambda is the old edge width.
    for lam in range(1, 5):
        for c in range(2 * lam + 1):
            row_volume = 2 * (lam * (col_side + unmatched) + c * rect)
            col_volume = 2 * (lam * (row_side + unmatched) + c * rect)
            formula = n * (n - 1) * lam - 2 * lam * k * comb(d, 2) \
                - 2 * rect * (2 * lam - c)
            assert row_volume == col_volume == formula
            lower_numerator = n * (n - 1) * lam \
                - n * lam * (d - 1) \
                - Fraction(n * (d - 1) ** 2 * (2 * lam - c), 2)
            assert formula >= lower_numerator
            if k > 1 and n > 1:
                assert Fraction(d - 1, n - 1) <= Fraction(1, k)
                assert Fraction((d - 1) ** 2, n - 1) \
                    <= Fraction(d - 1, k - 1)
            cases += 1
    return d, rect, cases


def main():
    regular = 0
    volume_cases = 0
    per_side = []
    for k in range(1, 5):
        counts = {}
        maxima = {}
        for mask in range(1 << (k * k)):
            rows = [(mask >> (k * a)) & ((1 << k) - 1) for a in range(k)]
            data = inspect_grid(rows, k)
            if data is None:
                continue
            d, rect, cases = data
            regular += 1
            volume_cases += cases
            counts[d] = counts.get(d, 0) + 1
            maxima[d] = max(maxima.get(d, 0), rect)
        per_side.append({"side": k, "nonempty_regular_supports": counts,
                         "maximum_rectangles": maxima})
    sharp = []
    for k in range(1, 13):
        for d in range(1, k + 1):
            if k % d:
                continue
            rows = [((1 << d) - 1) << ((a // d) * d) for a in range(k)]
            degree, rect, cases = inspect_grid(rows, k)
            assert degree == d
            assert 4 * rect == k * d * (d - 1) ** 2
            sharp.append({"side": k, "degree": d, "rectangles": rect})
    # Negative control: three occupied corners do not make a crossing.
    rows = [0b11, 0b01]
    assert sum(comb((rows[a] & rows[b]).bit_count(), 2)
               for a, b in combinations(range(2), 2)) == 0
    print(json.dumps({"all_binary_supports_through_side": 4,
                      "nonempty_regular_supports": regular,
                      "exact_channel_volume_cases": volume_cases,
                      "per_side": per_side,
                      "sharp_union_of_bicliques_cases": sharp,
                      "missing_corner_negative_control": True},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
