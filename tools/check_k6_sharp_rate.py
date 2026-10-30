#!/usr/bin/env python3
"""Finite exact checks supporting research/k6-sharp-rate.md.

These checks are arithmetic regression evidence, not an asymptotic proof.
"""

from fractions import Fraction
import json


def next_volume(k: int, g: int, width: int) -> int:
    upper_total = (12 * g**3 * width + 11 * g**4 * width + 24 * g**4) // 12
    return k * k * ((k * upper_total + g**4 * width) // k**3)


def verify() -> dict:
    assert 12**28 > 11**29  # beta_12 > 1/29
    assert 12**27 < 11**28  # negative control: beta_12 is not >= 1/28
    assert Fraction(1, 29) > Fraction(1, 30)
    assert Fraction(14 * 11, 10 * 48) / Fraction(11, 12) == Fraction(7, 20)
    assert Fraction(700, 48) < 15

    steps = 0
    for k in [48, 49, 50, 127, 600, 1024]:
        for g in [k, k + 1, 5 * k**6]:
            for width in list(range(48)) + [127, 1024, 10**6 + 3]:
                old = g**4 * width
                new = next_volume(k, g, width)
                assert 12 * k * new <= (11 * k + 24) * old + 24 * k * g**4
                steps += 1

    products = 0
    for start in range(12):
        product = Fraction(1)
        relative_sum = Fraction(0)
        for stop in range(start, 64):
            error = Fraction(1, 22 * 2**stop)
            product *= 1 + error
            relative_sum += error
            assert relative_sum < Fraction(1, 11)
            assert product <= 1 / (1 - relative_sum) <= Fraction(11, 10)
            products += 1

    # Exact closed form, without computing real logarithms.
    k = 48
    for d in range(4):
        assert k == 48 ** (12**d) * 25 ** ((12**d - 1) // 11)
        assert k >= 48 * 2**d
        k = 25 * k**12

    return {
        "scope": "finite exact arithmetic; ordinary asymptotic proof is separate",
        "one_step_cases": steps,
        "product_cases": products,
        "schedule_levels": 4,
        "beta_gt_one_over_29": True,
        "beta_lt_one_over_28_negative_control": True,
        "integer_margin_above_one_over_29": str(12**28 - 11**29),
        "integer_margin_below_one_over_28": str(11**28 - 12**27),
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
