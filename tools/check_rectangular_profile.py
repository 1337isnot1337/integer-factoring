#!/usr/bin/env python3
"""Exact arithmetic for the declared-profile minimax theorem.

This checks an analysis model, not tensor ranks or factoring runtimes.
No floating-point logarithms, square roots, or powers are used.
"""

from fractions import Fraction
from functools import lru_cache
import json


def combine_squared(a: Fraction, b: Fraction) -> Fraction:
    """Square of the minimax cost when the two child squares are a,b."""
    if a <= 0 or b <= 0:
        raise ValueError("positive child cost squares required")
    return a * b / (a + b)


@lru_cache(maxsize=None)
def tree_values(leaves: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    """Evaluate every ordered full binary tree with these leaf squares."""
    if not leaves:
        raise ValueError("at least one leaf required")
    if len(leaves) == 1:
        return leaves
    return tuple(combine_squared(a, b)
                 for cut in range(1, len(leaves))
                 for a in tree_values(leaves[:cut])
                 for b in tree_values(leaves[cut:]))


def verify() -> dict:
    half_cases = interior_cases = exponent_one_cases = three_quarter_cases = 0
    for a in range(1, 65):
        for b in range(1, 65):
            # p=2^a, d=2^b, alpha=1/2. These are cost SQUARES.
            ar, bc = Fraction(1, a), Fraction(1, b)
            r2, c2 = Fraction(a, a + b), Fraction(b, a + b)
            value = combine_squared(ar, bc)
            assert r2 + c2 == 1
            assert ar * r2 == bc * c2 == value == Fraction(1, a + b)
            # An exact lattice point gives this limiting Appendix B profile.
            h, s = a + b, 2 * (a - b)
            assert -2 * h < s < 2 * h
            assert Fraction(2 * h + s, 4 * h) == r2
            assert Fraction(2 * h - s, 4 * h) == c2
            half_cases += 1

            for scale2 in (Fraction(0), Fraction(1, 7), Fraction(1, 2), Fraction(1)):
                assert min(ar * r2 * scale2, bc * c2 * scale2) <= value
                interior_cases += 1

            # Negative control: alpha=1 exceeds the parent target.
            wrong = combine_squared(Fraction(1, a * a), Fraction(1, b * b))
            assert wrong == Fraction(1, a * a + b * b)
            assert wrong > Fraction(1, (a + b) ** 2)
            exponent_one_cases += 1

            # alpha=3/4, log p=a^2, log d=b^2. Remove the last square root
            # by squaring positive denominators: strict inequality is exact.
            assert (a**3 + b**3) ** 2 < (a*a + b*b) ** 3
            three_quarter_cases += 1

    tree_cases = 0
    for length in range(1, 10):
        leaves = tuple(Fraction(i + 2, i + 7) for i in range(length))
        expected = 1 / sum((1 / item for item in leaves), start=Fraction(0))
        values = tree_values(leaves)
        assert all(value == expected for value in values)
        tree_cases += len(values)

    # Uniform leaves attain the model's lower bound when K=K0^m.
    for count in range(1, 101):
        leaf_square = Fraction(1, 16)
        root_square = leaf_square / count
        assert root_square == 1 / (count / leaf_square)

    thin_cases = finite_thin_cases = 0
    for denominator in range(1, 33):
        a_cost = Fraction(1, 4)
        for numerator in range(1, denominator + 1):
            b_cost = a_cost * Fraction(numerator, denominator)
            z = b_cost / (1 + (b_cost / a_cost) ** 2)
            r, c = 2 * z / a_cost, 2 * z / b_cost - 1
            assert r >= 0 and c >= 0 and r * r + c * c == 1
            assert a_cost * r / 2 == b_cost * (1 + c) / 2 == z
            x = (a_cost / b_cost) ** 2
            x_next = (a_cost / z) ** 2
            assert x_next == x + 2 + 1 / x
            assert x + 2 <= x_next <= x + 3
            thin_cases += 1
            for d in (2, 3, 4, 16, 257, 2**30):
                factor = 1 - Fraction(1, d)
                row = a_cost * (Fraction(1, 2 * d) + factor * r / 2)
                column = b_cost * factor * (1 + c) / 2
                assert min(row, column) == factor * z
                finite_thin_cases += 1

    product = Fraction(1)
    recurrence_upper = Fraction(1)
    for step in range(1, 65):
        factor = 1 - Fraction(1, 2**step)
        product *= factor
        recurrence_upper = (recurrence_upper + 3) / factor**2
        assert product >= Fraction(1, 4)
        assert recurrence_upper <= (1 + 3 * step) / product**2
        assert recurrence_upper <= 16 * (1 + 3 * step)

    for bad in (Fraction(0), Fraction(-1)):
        try:
            combine_squared(bad, Fraction(1))
        except ValueError:
            pass
        else:
            raise AssertionError("nonpositive coefficient accepted")

    return {
        "scope": "declared channel-profile minimax model; no runtime or rank lower bound",
        "half_exponent_and_profile_cases": half_cases,
        "disk_interior_cases": interior_cases,
        "alpha_one_negative_controls": exponent_one_cases,
        "alpha_three_quarters_negative_controls": three_quarter_cases,
        "ordered_tree_cases": tree_cases,
        "uniform_leaf_cases": 100,
        "two_row_extremizer_and_recurrence_cases": thin_cases,
        "finite_two_row_profile_cases": finite_thin_cases,
        "finite_two_row_product_and_iteration_cases": 64,
        "invalid_cost_negative_controls": 2,
        "arithmetic": "integer and Fraction only",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
