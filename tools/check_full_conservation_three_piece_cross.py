#!/usr/bin/env python3
"""Close the three-piece full-conservation windows at q=23,25,27."""
from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from math import floor, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'evidence' / 'full-conservation-three-piece-cross.json'
LABELS = (23, 25, 27)


def quotient_rank(q: int, k: int) -> int:
    return max(s * (q - k * s) for s in range(q + 1))


def minimum_column_for(q: int, budget: int) -> int:
    for c in range(2, q + 1):
        if quotient_rank(q, c) <= budget:
            return c
    raise AssertionError('no admissible minimum column')


def two_piece_bound(q: int) -> float:
    return (sqrt(q * q - q + 1)
            + sqrt((q * q - 1) * (q - 1)))


def polygon_vertices(r1: int, r2: int, q: int) -> list[tuple[Fraction, Fraction]]:
    n = q * q
    a = n - quotient_rank(q, r1)
    b = n - quotient_rank(q, r2)
    c = n - minimum_column_for(q, r1 + r2)
    lines = [
        (Fraction(1), Fraction(0), Fraction(2)),
        (Fraction(0), Fraction(1), Fraction(2)),
        (Fraction(1), Fraction(0), Fraction(a)),
        (Fraction(0), Fraction(1), Fraction(b)),
        (Fraction(1), Fraction(1), Fraction(c)),
    ]
    vertices: set[tuple[Fraction, Fraction]] = set()
    for i, (a1, b1, c1) in enumerate(lines):
        for a2, b2, c2 in lines[i + 1:]:
            det = a1 * b2 - a2 * b1
            if det == 0:
                continue
            x = (c1 * b2 - c2 * b1) / det
            y = (a1 * c2 - a2 * c1) / det
            if (x >= 2 and y >= 2 and x <= a and y <= b
                    and x + y <= c):
                vertices.add((x, y))
    if not vertices:
        raise AssertionError(f'empty polygon at q={q}, r1={r1}, r2={r2}')
    return sorted(vertices)


def minimum_relaxed_cost(q: int) -> tuple[float, tuple[int, int, int, int, int]]:
    n = q * q
    target = q ** sqrt(2)
    cap = floor(target - 2 * q)
    assert 4 * q > target
    assert q * q > target
    assert two_piece_bound(q) > target
    best = (float('inf'), None)
    for r1 in range(2, cap + 1):
        for r2 in range(2, cap + 1):
            r3 = n - r1 - r2
            assert r3 > q
            for c1, c2 in polygon_vertices(r1, r2, q):
                c3 = n - c1 - c2
                assert c3 >= 2
                cost = (sqrt(r1 * float(c1))
                        + sqrt(r2 * float(c2))
                        + sqrt(r3 * float(c3)))
                if cost < best[0]:
                    best = (cost, (r1, r2, int(c1), int(c2), int(c3)))
    return best


def main() -> None:
    rows = []
    for q in LABELS:
        target = q ** sqrt(2)
        cap = floor(target - 2 * q)
        best_cost, profile = minimum_relaxed_cost(q)
        margin = best_cost - target
        assert margin > 1.0
        rows.append({
            'q': q,
            'target_q_sqrt2': target,
            'piece_floor_q': q,
            'three_piece_lower_bound': best_cost,
            'three_piece_profile': profile,
            'margin': margin,
            'small_rank_cap': cap,
            'two_piece_bound': two_piece_bound(q),
            'four_piece_floor': 4 * q,
        })
    receipt = {
        'status': 'pass',
        'labels': list(LABELS),
        'rows': rows,
        'scope': (
            'For q=23,25,27 an exact crossing decomposition below q^sqrt(2) '
            'has at most three pieces. The identity costs q^2, the floor q '
            'excludes at least four pieces, the all-field two-piece bound '
            'excludes two pieces, '
            'and the cross-quotient vertex relaxation excludes three pieces. '
            'The three-piece check transposes if necessary so two pieces have '
            'row rank at most floor(q^sqrt(2)-2q); all q=23,25,27 cases clear '
            'the target with margin >1. No other q and no score change.'
        ),
        'source_sha256': {
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
        'profiles': {row['q']: row['three_piece_profile'] for row in rows},
    }))


if __name__ == '__main__':
    main()
