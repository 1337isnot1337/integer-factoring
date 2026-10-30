#!/usr/bin/env python3
"""Finite exact checks for research/graded-space-barrier.md.

Matrix ranks are over Q. Large integer inequalities check finite rank
budgets, not asymptotic runtimes or a new crossing construction.
"""

from fractions import Fraction
import json


def rational_rank(matrix: list[list[int]]) -> int:
    rows = [[Fraction(value) for value in row] for row in matrix]
    if not rows:
        return 0
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("ragged matrix")
    pivot = 0
    for column in range(width):
        selected = next((r for r in range(pivot, len(rows)) if rows[r][column]), None)
        if selected is None:
            continue
        rows[pivot], rows[selected] = rows[selected], rows[pivot]
        divisor = rows[pivot][column]
        rows[pivot] = [value / divisor for value in rows[pivot]]
        for r in range(pivot + 1, len(rows)):
            scale = rows[r][column]
            if scale:
                rows[r] = [value - scale * other
                           for value, other in zip(rows[r], rows[pivot])]
        pivot += 1
        if pivot == len(rows):
            break
    return pivot


def crossing_matrix(labels: int, scale: int) -> list[list[int]]:
    if labels < 1:
        raise ValueError("a positive label count is required")
    return [[scale * (x == v and y == u)
             for u in range(labels) for v in range(labels)]
            for x in range(labels) for y in range(labels)]


def growth_allowed_by_rank(window: int, numerator: int, denominator: int) -> bool:
    if min(window, numerator, denominator) < 1:
        raise ValueError("positive integer parameters are required")
    return window**numerator >= 2**denominator


def verify() -> dict:
    rank_cases = 0
    for labels in range(1, 10):
        for scale in (-3, 1, 7):
            assert rational_rank(crossing_matrix(labels, scale)) == labels * labels
            rank_cases += 1
        assert rational_rank(crossing_matrix(labels, 0)) == 0

    examples = []
    cases = [(4, 1, 4, False), (4, 1, 3, False), (4, 2, 5, False),
             (4, 3, 7, False), (4, 49, 100, False), (4, 1, 2, True),
             (4, 3, 5, True), (8, 1, 3, True), (8, 1, 4, False),
             (9, 1, 3, True)]
    for window, numerator, denominator, expected in cases:
        actual = growth_allowed_by_rank(window, numerator, denominator)
        assert actual is expected
        examples.append({
            "window": window,
            "gamma": f"{numerator}/{denominator}",
            "left_integer": str(window**numerator),
            "right_integer": str(2**denominator),
            "passes_necessary_condition_only": actual,
        })

    # A first-order deficit eventually defeats a polynomial number of channels.
    # Here M=(2*lambda+1)^2, r_i<=2^ceil(1.9*lambda), m=2^lambda.
    rejected_budgets = []
    for bits in (256, 512, 1024, 4096):
        pieces = (2 * bits + 1)**2
        channels = 1 << ((19 * bits + 9) // 10)
        required = 1 << (2 * bits)
        assert pieces * channels < required
        # The same count with full row rank is not excluded by this criterion.
        assert pieces * required >= required
        rejected_budgets.append({
            "lambda": bits,
            "pieces": pieces,
            "channel_exponent": (19 * bits + 9) // 10,
            "required_rank_exponent": 2 * bits,
            "budget_rejected": True,
        })

    return {
        "scope": "finite rank and integer checks; the asymptotic restriction has a written proof",
        "nonzero_scaled_crossing_rank_cases": rank_cases,
        "zero_scale_negative_control_cases": 9,
        "dimension_growth_cases": examples,
        "finite_polynomial_piece_budgets": rejected_budgets,
        "claims_new_factoring_exponent": False,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
