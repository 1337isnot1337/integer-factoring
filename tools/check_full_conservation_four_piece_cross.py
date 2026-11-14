#!/usr/bin/env python3
"""Close the first four-piece full-conservation windows."""
from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from math import floor, sqrt
from pathlib import Path

from check_full_conservation_three_piece_cross import (
    minimum_column_for,
    polygon_vertices,
    quotient_rank,
    two_piece_bound,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'evidence' / 'full-conservation-four-piece-cross.json'
LABELS = (29, 31, 33, 35, 37, 39)


def three_piece_minimum(q: int) -> tuple[float, tuple[int, int, int, int, int]]:
    n = q * q
    target = q ** sqrt(2)
    cap = floor(target - 2 * q)
    best = (float('inf'), None)
    for r1 in range(2, cap + 1):
        for r2 in range(2, cap + 1):
            r3 = n - r1 - r2
            assert r3 > q
            for c1, c2 in polygon_vertices(r1, r2, q):
                c3 = n - c1 - c2
                if c3 < 2:
                    continue
                cost = (sqrt(r1 * float(c1))
                        + sqrt(r2 * float(c2))
                        + sqrt(r3 * float(c3)))
                if cost < best[0]:
                    best = (cost, (r1, r2, int(c1), int(c2), int(c3)))
    return best


def group_minimum(q: int, n: int, small: tuple[int, ...], total: int) -> float:
    count = len(small)
    if count == 0:
        return 0.0 if total == 0 else float('inf')
    lower = list(small)
    upper = [n - quotient_rank(q, m) for m in small]
    if total < sum(lower) or total > sum(upper):
        return float('inf')
    best = float('inf')
    for free in range(count):
        for mask in range(1 << (count - 1)):
            subtotal = 0.0
            cost = 0.0
            bit = 0
            for i in range(count):
                if i == free:
                    continue
                value = upper[i] if (mask >> bit) & 1 else lower[i]
                subtotal += value
                cost += sqrt(small[i] * value)
                bit += 1
            free_value = total - subtotal
            if lower[free] <= free_value <= upper[free]:
                candidate = cost + sqrt(small[free] * free_value)
                if candidate < best:
                    best = candidate
    return best


def four_piece_minimum(q: int) -> tuple[float, tuple[tuple[int, ...], int, tuple[int, ...], tuple[int, ...]]]:
    n = q * q
    target = q ** sqrt(2)
    cap = floor(target - 3 * q)
    assert 5 * q > target
    best = (float('inf'), None)
    profiles = 0
    for small in combinations_with_replacement(range(2, cap + 1), 4):
        if sum(small) < q:
            continue
        profiles += 1
        for mask in range(8):
            rows = tuple(small[i] for i in range(4) if (mask >> i) & 1)
            columns = tuple(small[i] for i in range(4) if not ((mask >> i) & 1))
            row_small_sum = sum(rows)
            column_small_sum = sum(columns)
            row_large = group_minimum(q, n, rows, n - column_small_sum)
            if row_large == float('inf'):
                continue
            column_large = group_minimum(q, n, columns, n - row_small_sum)
            if column_large == float('inf'):
                continue
            cost = row_large + column_large
            if cost < best[0]:
                best = (cost, (small, mask, rows, columns))
    return best, profiles


def main() -> None:
    rows = []
    for q in LABELS:
        target = q ** sqrt(2)
        one_piece = q * q
        two_piece = two_piece_bound(q)
        three_cost, three_profile = three_piece_minimum(q)
        four_result, profiles = four_piece_minimum(q)
        four_cost, four_profile = four_result
        assert one_piece > target
        assert two_piece > target
        assert three_cost > target + 1
        assert four_cost > target + 1
        assert 5 * q > target
        rows.append({
            'q': q,
            'target_q_sqrt2': target,
            'one_piece_cost': one_piece,
            'two_piece_bound': two_piece,
            'three_piece_lower_bound': three_cost,
            'three_piece_profile': three_profile,
            'four_piece_lower_bound': four_cost,
            'four_piece_profile': {
                'small_ranks': list(four_profile[0]),
                'row_small_mask': four_profile[1],
                'row_small': list(four_profile[2]),
                'column_small': list(four_profile[3]),
            },
            'four_piece_profiles': profiles,
            'small_rank_cap': floor(target - 3 * q),
            'five_piece_floor': 5 * q,
            'margin': four_cost - target,
        })
    receipt = {
        'status': 'pass',
        'labels': list(LABELS),
        'rows': rows,
        'scope': (
            'For odd q=29,31,33,35,37,39 an exact crossing decomposition below '
            'q^sqrt(2) has at most four pieces. One, two, and three pieces are '
            'excluded by q^2, the all-field two-piece bound, and the three-piece '
            'cross relaxation. For four pieces, every piece has both ranks at least '
            'two and smaller rank at most floor(q^sqrt(2)-3q). Writing the larger '
            'rank L and smaller rank m, the quotient bound gives L<=q^2-f(m), and '
            'the two total-rank budgets fix the sums of L in the row-small and '
            'column-small groups. The checker enumerates every small-rank multiset '
            'and orientation partition, and minimizes the concave group costs at '
            'the bound vertices. All cases clear the target with margin >1. This is '
            'a finite cross-constraint closure, not a score change and not a claim '
            'for even q or q>=41.'
        ),
        'source_sha256': {
            'three_piece_tool': sha256((ROOT / 'tools' /
                                        'check_full_conservation_three_piece_cross.py').read_bytes()).hexdigest(),
            'two_piece_note': sha256((ROOT / 'research' /
                                      'two-piece-rank-sum-obstruction.md').read_bytes()).hexdigest(),
            'minimum_cut_note': sha256((ROOT / 'research' /
                                        'minimum-cut-rank-budget.md').read_bytes()).hexdigest(),
        },
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({
        'passed': True,
        'closed_q': list(LABELS),
        'minimum_margins': {row['q']: row['margin'] for row in rows},
        'four_piece_profiles': {row['q']: row['four_piece_profile'] for row in rows},
    }))


if __name__ == '__main__':
    main()
