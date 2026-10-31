#!/usr/bin/env python3
"""Exact entropy duals for a repeated minimum Smith order on 7--10 letters."""
from fractions import Fraction as Q
import json
from math import lcm, log


def main():
    records = []
    alpha = (Q(1, 2), Q(1, 5), Q(3, 10))
    for q in range(7, 11):
        n = q*q
        multiplier = alpha[1]/(q-2)-alpha[2]/(n-q)
        beta2 = Q(n-q, n)*(1-q*multiplier)-alpha[2]
        beta1 = (q-1)*(multiplier+beta2/(n-q))
        beta = (1-beta1-beta2, beta1, beta2)
        assert sum(alpha) == sum(beta) == 1 and min(alpha+beta) > 0
        assert multiplier > 0
        assert beta[1]/(q-1)-beta[2]/(n-q) == multiplier
        assert 1-Q(n, n-q)*(alpha[2]+beta[2]) == q*multiplier
        denominator = lcm(*(weight.denominator for weight in alpha+beta))
        entropy_power = Q(1)
        approximation = 0
        for weights, multiplicities in ((alpha, (2, q-2, n-q)), (beta, (1, q-1, n-q))):
            for weight, multiplicity in zip(weights, multiplicities):
                exponent = weight*denominator
                assert exponent.denominator == 1
                entropy_power *= (Q(multiplicity)/weight)**int(exponent)
                approximation += float(weight)*log(multiplicity/float(weight))
        # exp(C)>q^(17/6), so rho>C/(4 log q)>17/24>1/sqrt(2).
        assert entropy_power**6 > q**(17*denominator)
        records.append({"alphabet": q, "alpha": list(map(str, alpha)), "beta": list(map(str, beta)),
                        "shared_prefix_multiplier": str(multiplier),
                        "common_weight_denominator": denominator,
                        "exact_entropy_comparison": "exp(6*D*C)>q^(17*D)",
                        "approximate_rho_lower_bound": approximation/(4*log(q))})
    assert 2*17**2-24**2 == 2
    print(json.dumps({"scope": "Bare profile, either minimum-order multiplicity at least two",
                      "rho_lower_bound": "17/24", "certificates": records}, indent=2))


if __name__ == '__main__':
    main()
