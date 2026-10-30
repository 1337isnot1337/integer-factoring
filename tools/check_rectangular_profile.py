#!/usr/bin/env python3
"""Exact arithmetic for the declared-profile minimax theorem.

This checks an analysis model, not tensor ranks or factoring runtimes.
No floating-point logarithms, square roots, or powers are used.
"""

from fractions import Fraction
from functools import lru_cache
import json
from math import isqrt


def log_interval(n: int) -> tuple[Fraction, Fraction]:
    """Rational enclosure of natural log(n), using a positive atanh series."""
    if n < 1:
        raise ValueError("positive integer required")

    def series(z: Fraction) -> tuple[Fraction, Fraction]:
        assert 0 <= z <= Fraction(1, 3)
        lower = 2 * sum((z ** (2 * j + 1) / (2 * j + 1)
                         for j in range(20)), start=Fraction(0))
        tail = 2 * z**41 / (41 * (1 - z*z))
        return lower, lower + tail

    power = n.bit_length() - 1
    two_lower, two_upper = series(Fraction(1, 3))
    z = Fraction(n - (1 << power), n + (1 << power))
    lower, upper = series(z)
    return power * two_lower + lower, power * two_upper + upper


def sqrt_interval(value: Fraction) -> tuple[Fraction, Fraction]:
    """Enclose a nonnegative rational square root using integers only."""
    if value < 0:
        raise ValueError("nonnegative rational required")
    scale = 10**10
    root = isqrt(value.numerator * scale * scale // value.denominator)
    lower, upper = Fraction(root, scale), Fraction(root + 1, scale)
    assert lower * lower <= value < upper * upper
    return lower, upper


def geometry_interval(p: int, d: int, delta: Fraction,
                      logs: dict[int, tuple[Fraction, Fraction]]) -> tuple[Fraction, Fraction]:
    """Enclose R(p,d) from equation (15) with the correct endpoint directions."""
    pl, ph = logs[p]
    dl, dh = logs[d]
    low_l, high_l = pl + dl - delta, ph + dh - delta
    assert min(pl, dl) > delta and low_l > 0
    x, y = Fraction(p - 1, p * d), Fraction(d - 1, p * d)
    row_lower = sqrt_interval((pl - delta) / high_l)[0]
    row_upper = sqrt_interval((ph - delta) / low_l)[1]
    col_lower = sqrt_interval((dl - delta) / high_l)[0]
    col_upper = sqrt_interval((dh - delta) / low_l)[1]
    lower = 2*x*dl/(1 + row_upper) + 2*y*pl/(1 + col_upper) + low_l/(p*p*d*d)
    upper = 2*x*dh/(1 + row_lower) + 2*y*ph/(1 + col_lower) + high_l/(p*p*d*d)
    assert lower <= upper
    return lower, upper


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

    delta = Fraction(2, 3)
    logs = {n: log_interval(n) for n in range(2, 64)}
    geometry_cases = 0
    worst_upper = Fraction(0)
    worst_pair = None
    for p in range(2, 64):
        for d in range(2, 64):
            lower, upper = geometry_interval(p, d, delta, logs)
            assert upper <= delta - Fraction(1, 500)
            x, y = Fraction(p - 1, p * d), Fraction(d - 1, p * d)
            q = Fraction((p - 1) * (d - 1), p * d)
            assert 1 + x*x + y*y - q*q == 2*x + 2*y + Fraction(1, p*p*d*d)
            if upper > worst_upper:
                worst_upper, worst_pair = upper, [p, d]
            geometry_cases += 1
    # Analytic tail for max(p,d)>=64. The transcendental facts are proved
    # in the note; these comparisons are its exact rational arithmetic.
    assert Fraction(11, 4)**4 < 64 < Fraction(8, 3)**5
    assert Fraction(5, 11) > Fraction(4, 9)
    tail = Fraction(9, 20) + Fraction(5, 32) + Fraction(5, 8192)
    assert tail == Fraction(24857, 40960) < delta
    assert logs[2][0] - delta >= Fraction(2, 81)
    assert Fraction(81, 2048) < Fraction(1, 16)
    failed_shift = Fraction(33, 50)
    assert geometry_interval(4, 4, failed_shift, logs)[0] > failed_shift

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
        "all_finite_geometry_certificate_cases": geometry_cases,
        "finite_geometry_worst_upper_pair": worst_pair,
        "finite_geometry_uniform_margin": "1/500",
        "finite_geometry_tail_bound": "24857/40960 < 2/3",
        "finite_geometry_failed_shift_control": "33/50 at p=d=4",
        "invalid_cost_negative_controls": 2,
        "arithmetic": "integer and Fraction only",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
