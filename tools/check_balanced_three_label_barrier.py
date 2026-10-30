#!/usr/bin/env python3
"""Exact degeneration and entropy certificates for the original PR24 laser."""

import argparse
from collections import Counter
from itertools import product
import json
from math import comb, factorial
from pathlib import Path

from check_binary_coefficient_entropy import coefficient_rank
from check_finite_seed_amplification import (
    FLAGS, PIECES, SOURCE_COMMIT, column, conjugate, crossing, multiply, rank,
)


HEIGHTS = {
    "row": (0, 1, 2, 3, 4, 6, 7, 8, 9),
    "column": (0, 2, 3, 4, 8, 6, 10, 11, 12),
}
EXPECTED = {
    "row": {
        4: (1, ((2, 3), (7, 8))),
        5: (2, ((0, 2), (5, 7))),
        12: (0, ((4, 4),)),
        15: (3, ((3, 5),)),
        16: (0, ((0, 0), (3, 3), (5, 5), (8, 8))),
        17: (0, ((1, 1), (6, 6))),
    },
    "column": {
        4: (2, ((1, 3), (6, 8))),
        5: (2, ((0, 1), (4, 6))),
        12: (0, ((5, 5),)),
        15: (0, ((0, 0), (3, 3), (4, 4), (8, 8))),
        16: (4, ((3, 4),)),
        17: (0, ((2, 2), (7, 7))),
    },
}


def multiset_words(counts):
    """Yield distinct words, without factorial division in the base field."""
    remaining = dict(counts)
    size = sum(remaining.values())
    keys = sorted(remaining)

    def visit(prefix):
        if len(prefix) == size:
            yield tuple(prefix)
            return
        for key in keys:
            if remaining[key]:
                remaining[key] -= 1
                prefix.append(key)
                yield from visit(prefix)
                prefix.pop()
                remaining[key] += 1

    yield from visit([])


def base_n(word, base):
    result = 0
    for value in word:
        result = result * base + value
    return result


def restricted_initial_matrix(side, counts, w_count, patterns):
    copies = sum(counts.values())
    allowed = []
    for word in product(range(9), repeat=copies):
        if all(sum(patterns[k][j] for k in word) == w_count
               for j in range(2)):
            allowed.append(base_n(word, 9))
    index = {word: i for i, word in enumerate(allowed)}
    rows = [0] * len(allowed)
    for word in multiset_words(counts):
        entries = [EXPECTED[side][piece][1] for piece in word]
        for cells in product(*entries):
            output = index.get(base_n((r for r, _ in cells), 9))
            input_ = index.get(base_n((c for _, c in cells), 9))
            if output is not None and input_ is not None:
                rows[output] ^= 1 << input_
    return rows


def verify():
    tables = {}
    patterns = {}
    checked_cells = 0
    for side, heights in HEIGHTS.items():
        raw = PIECES if side == "row" else tuple(column(A, 3) for A in PIECES)
        patterns[side] = []
        for vector in FLAGS[side]:
            support = {(x == 2, y == 2)
                       for x in range(3) for y in range(3)
                       if (vector >> (3*x+y)) & 1}
            assert len(support) == 1
            patterns[side].append(support.pop())
        table = {}
        for piece, (valuation, expected) in EXPECTED[side].items():
            A = conjugate(multiply(raw[piece], crossing(3)), FLAGS[side])
            cells = [(heights[c]-heights[r], r, c)
                     for r, row in enumerate(A) for c in range(9)
                     if (row >> c) & 1]
            assert min(v for v, _, _ in cells) == valuation
            leading = tuple((r, c) for v, r, c in cells if v == valuation)
            assert leading == expected
            initial = [sum(1 << c for rr, c in leading if rr == r)
                       for r in range(9)]
            expected_rank = {4: 2, 5: 2, 12: 1, 17: 2}.get(piece)
            if expected_rank is not None:
                assert rank(initial) == expected_rank
            if piece in (4, 5, 12, 17):
                out_pattern, in_pattern = {
                    4: ((True, False), (False, False)),
                    5: ((False, False), (True, False)),
                    12: ((True, True), (True, True)),
                    17: ((False, True), (False, True)),
                }[piece]
                assert all(patterns[side][r] == out_pattern and
                           patterns[side][c] == in_pattern for r, c in leading)
            else:
                assert all(patterns[side][r] == (False, False) and
                           patterns[side][c] == (False, False)
                           for r, c in leading)
                identity_piece = 16 if side == "row" else 15
                if piece == identity_piece:
                    uu = [j for j, pattern in enumerate(patterns[side])
                          if pattern == (False, False)]
                    assert leading == tuple((j, j) for j in uu)
                else:
                    assert rank(initial) == 1
                    assert multiply(initial, initial) == (0,)*9
            checked_cells += 81
            table[str(piece)] = {"valuation": valuation, "initial_edges": leading}
        tables[side] = {"heights": heights, "flag": FLAGS[side],
                        "w_patterns": patterns[side], "pieces": table}

    finite = []
    # (common-w positions, core copies, binary positions, nilpotent count).
    cases = ((0, 1, 0, 0), (1, 1, 0, 0),
             (1, 0, 3, 1), (1, 0, 3, 2),
             (1, 0, 4, 2), (0, 1, 2, 1))
    for common, core, binary, nil_count in cases:
        copies = common + 3*core + binary
        w_count = common + core
        counts = Counter({12: common, 4: core, 5: core, 17: core,
                          15: nil_count, 16: binary-nil_count})
        counts = {i: amount for i, amount in counts.items() if amount}
        label_dimension = comb(copies, w_count)*2**(copies-w_count)
        core_factor = (comb(copies, common)*factorial(3*core)
                       // factorial(core)**3)*8**core
        # N^i I^(binary-i) has rank r_(binary-i); the other cut reverses it.
        binary_ranks = (coefficient_rank(binary, binary-nil_count),
                        coefficient_rank(binary, nil_count))
        result = []
        for j, side in enumerate(("row", "column")):
            rows = restricted_initial_matrix(side, counts, w_count, patterns[side])
            assert len(rows) == label_dimension**2
            actual = rank(rows)
            lower = core_factor*binary_ranks[j]
            assert actual >= lower
            result.append({"rank": actual, "direct_sum_lower_bound": lower})
        finite.append({"copies": copies, "w_count": w_count,
                       "label_dimension": label_dimension, "piece_counts": counts,
                       "cuts": result})

    # These integer comparisons certify the entropy tangent argument, with no
    # floating-point logarithms or numerical minimization in the certificate.
    assert 3**375 > 2**594
    assert 2**242776 * 761**38050 > 3**382959
    assert 18**2 * 2 > 25**2  # 18/25 > 1/sqrt(2).
    binary_bounds = 0
    for copies in range(65):
        index = copies//2
        lower = comb(2*copies, index)
        for side_index in {index, copies-index}:
            assert (copies+1)*coefficient_rank(copies, side_index) >= lower
            binary_bounds += 1

    return {
        "source_commit": SOURCE_COMMIT,
        "scope": "Literal type coefficients; original U=span(e0,e1), w=e2 only.",
        "degeneration_cells_checked": checked_cells,
        "degenerations": tables,
        "finite_restricted_initial_profiles": finite,
        "binary_middle_bounds": binary_bounds,
        "entropy_certificate": {
            "rho": [18, 25], "common_fraction": [32, 761],
            "slope_comparison": "3^375 > 2^594",
            "intercept_comparison": "2^242776 * 761^38050 > 3^382959",
            "comparisons_exact": True,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": "pass",
                      "degeneration_cells": result["degeneration_cells_checked"],
                      "finite_profiles": len(result["finite_restricted_initial_profiles"]),
                      "binary_middle_bounds": result["binary_middle_bounds"],
                      "rho_certificate": "18/25"}))


if __name__ == "__main__":
    main()
