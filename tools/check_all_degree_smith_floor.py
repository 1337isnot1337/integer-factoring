#!/usr/bin/env python3
"""Exact controls for the arbitrary-degree small-alphabet pressure theorem."""
import json
from math import log

from check_quadratic_crossings import ROOT, convolution


def check_profile(row, regularization_degree):
    mu, nu = row["infinity"], row["zero"]
    assert len(mu) == len(nu) == 4
    assert sorted(mu) == mu and sorted(nu) == nu
    common_order = min(nu)
    infinity_degrees = [regularization_degree - order - common_order for order in mu]
    zero_orders = [order-common_order for order in nu]
    normalized_degree = max(infinity_degrees)
    assert normalized_degree >= 0
    assert sum(infinity_degrees) >= 0
    assert min(zero_orders) == 0
    assert sum(zero_orders) <= 4*normalized_degree


def main():
    a, b = 577, 816
    assert 2*a*a-b*b == 2
    comparisons = []
    for q in (2, 3, 4):
        n = q*q
        left = (4*n*n)**(n*b)
        right = (n+1)**((n+1)*b) * q**(4*a*n)
        assert left > right
        floor = 4*n*n/(n+1)**((n+1)/n)
        comparisons.append({"alphabet": q, "matrix_rank": n,
                            "exact_comparison": "F_n > q^(4*577/816)",
                            "left_bit_length": left.bit_length(),
                            "right_bit_length": right.bit_length(),
                            "pressure_floor_approximate": floor,
                            "rho_floor_approximate": log(floor)/(4*log(q))})
    paths = ["quadratic-crossing-catalogue.json", "parity-cubic-crossing-catalogue.json",
             "parity-quartic-crossing-catalogue.json", "cubic-crossing-million-sample.json",
             "rank-one-constant-quadratic-catalogues.json"]
    checks = []
    for name in paths:
        data = json.loads((ROOT / "evidence" / name).read_text())
        catalogues = data if isinstance(data, list) else [data]
        for catalogue in catalogues:
            for row in catalogue["profiles"]:
                check_profile(row, catalogue["degree"])
        checks.append({"catalogue": name,
                       "profiles_checked": sum(len(d["profiles"]) for d in catalogues),
                       "normalization_and_four_constraints": "passed"})
    # Exact sharp scalar parity inequality, valid for every y in (0,1].
    assert convolution([1, -6, 9], [27, 14, 3]) == [27, -148, 162, 108, 27]
    print(json.dumps({"scope": "Ordinary all-degree proof plus exact finite controls; no factoring improvement",
                      "critical_rational": {"numerator": a, "denominator": b,
                                            "strict_square_comparison": 2,
                                            "approximate": a/b},
                      "comparisons": comparisons, "profile_controls": checks,
                      "all_field_scalar_parity_floor": "256/27"}, indent=2))


if __name__ == "__main__":
    main()
