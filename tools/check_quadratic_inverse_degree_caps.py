#!/usr/bin/env python3
"""Exact entropy-slack degree caps for scalar quadratic-inverse crossings."""
from fractions import Fraction as Q
import json

from check_leading_rank_entropy_floor import entropy_dual


def main():
    cutoffs = {7: Q(3, 2), 8: Q(269, 200), 9: Q(239, 200), 10: Q(130, 121)}
    records = []
    for q, cutoff in cutoffs.items():
        alpha, _, _, denominator, entropy_power, _ = entropy_dual(q)
        penalty = alpha[0]-2*alpha[1]/(q-2)
        assert penalty > 0 and (denominator*penalty).denominator == 1
        adjusted = entropy_power/cutoff**int(denominator*penalty)
        assert adjusted**6 > q**(17*denominator)
        minimum_degree = 1
        while cutoff**(6*minimum_degree) <= q**17:
            minimum_degree += 1
        assert cutoff**(6*(minimum_degree-1)) <= q**17
        assert q*q-10 > 0  # c=0 forces a second-order channel.
        records.append({'alphabet': q, 'scale_split': str(cutoff),
                        'entropy_penalty': str(penalty),
                        'minimum_excluded_degree': minimum_degree,
                        'general_inverse_degree_ratio_cap': minimum_degree,
                        'rho_lower_bound': '17/24',
                        'small_scale_exact_comparison': '(exp(D*H)/a^(D*kappa))^6>q^(17*D)',
                        'large_scale_exact_comparison': 'a^(6*M)>q^17'})
    assert [r['minimum_excluded_degree'] for r in records] == [14, 20, 35, 91]
    print(json.dumps({'scope': 'Scalar polynomial crossing with polynomial inverse; general aspect cap and stronger quadratic cap; bare profiles',
                      'certificates': records}, indent=2))


if __name__ == '__main__':
    main()
