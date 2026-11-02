#!/usr/bin/env python3
"""A seven-letter asymmetric Smith profile compatible with inverse degree."""
import json
from fractions import Fraction

Q = 7
M = 12
DIMENSION = Q * Q
MU = [0] + [7] * 6 + [13] * 42
NU = [0] + [7] * 6 + [12] * 3 + [13] * 38 + [16]


def pressure(s):
    return sum(s ** (M - order) for order in MU) * sum(s ** (-order) for order in NU)


def main():
    assert len(MU) == len(NU) == DIMENSION
    assert all(MU[i] <= MU[i + 1] for i in range(DIMENSION - 1))
    assert all(NU[i] <= NU[i + 1] for i in range(DIMENSION - 1))
    assert MU[0] == NU[0] == 0
    assert MU[1] == 7
    assert sum(MU) == sum(NU) == DIMENSION * M
    assert sum(MU[:Q]) + sum(NU[:Q]) == Q * M
    assert MU[1] <= 2 * (NU[-1] - M)
    s = Fraction(34, 25)
    value = pressure(s)
    assert value.numerator ** 169 < 7 ** 478 * value.denominator ** 169
    assert 478 ** 2 < 8 * 169 ** 2
    report = {
        "alphabet": Q,
        "dimension": DIMENSION,
        "actual_degree": M,
        "infinity_orders": MU,
        "zero_orders": NU,
        "sums": [sum(MU), sum(NU)],
        "shared_prefix_sum": sum(MU[:Q]) + sum(NU[:Q]),
        "inverse_degree_bound": 2 * (NU[-1] - M),
        "second_infinity_order": MU[1],
        "sample_scale": f"{s.numerator}/{s.denominator}",
        "sample_pressure": f"{value.numerator}/{value.denominator}",
        "threshold_exponent": "478/169",
        "exact_threshold_comparison": "pressure^169 < 7^478",
        "threshold_below_two_sqrt_two": "478^2 < 8*169^2",
        "crosses_profile_threshold": True,
        "scope": "necessary Smith profile only; no algebraic crossing is constructed",
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
