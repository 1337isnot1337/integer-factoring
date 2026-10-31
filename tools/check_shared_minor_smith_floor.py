#!/usr/bin/env python3
"""Exact all-scale checks for the shared-minor small-alphabet theorem."""
import json
from math import comb

from check_quadratic_crossings import bernstein_certificate


def main():
    # U_q(x) is at least these rational numbers for every x>=1.
    targets = {2: (297, 100), 3: (526, 100), 4: (773, 100),
               5: (1029, 100), 6: (12933, 1000)}
    result = []
    for q, (numerator, denominator) in targets.items():
        # x^q(U_q(x)-T), cleared of denominators.
        polynomial = [0] * (2*q)
        polynomial[2*q-1] += denominator
        polynomial[q-1] += (q-1)*denominator
        polynomial[0] += q*(q-1)*denominator
        polynomial[q] -= numerator
        # Substitute x=1+2z, so z in [0,1] covers x in [1,3].
        transformed = [sum(polynomial[i]*comb(i, j)*2**j for i in range(j, 2*q))
                       for j in range(2*q)]
        certificate = bernstein_certificate(transformed)
        # The remaining ray is controlled by the first positive term alone.
        assert denominator*3**(q-1) > numerator
        # T^2 > q^(20/7) gives rho > 5/7 > 1/sqrt(2).
        assert numerator**14 > q**20 * denominator**14
        result.append({"alphabet": q, "sqrt_pressure_lower_bound": f"{numerator}/{denominator}",
                       "compact_interval": [1, 3], "bernstein_certificate": certificate,
                       "tail": "x^(q-1)>T for x>=3",
                       "exact_threshold_comparison": "T^14>q^20"})
    assert 2*5**2 > 7**2
    print(json.dumps({"scope": "Exact all-scale certificates supporting the ordinary shared-minor proof",
                      "rho_lower_bound": "5/7", "certificates": result}, indent=2))


if __name__ == "__main__":
    main()
