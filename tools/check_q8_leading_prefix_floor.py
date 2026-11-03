#!/usr/bin/env python3
"""Exact scalar certificate for a q8 barrier from only fifty protected orders."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = ROOT/'evidence/q8-leading-prefix-floor.json'
DEGREE = 4096


def polynomial():
    # [1+14u^30+35u^60+4046u^285]^7 - [1+7u^56]^20.
    result = {0: 1}
    for _ in range(7):
        following = {}
        for exponent, coefficient in result.items():
            for shift, weight in ((0, 1), (30, 14), (60, 35), (285, 4046)):
                key = exponent+shift
                following[key] = following.get(key, 0)+coefficient*weight
        result = following
    for j in range(21):
        key = 56*j
        result[key] = result.get(key, 0)-comb(20, j)*7**j
    result = {j: c for j, c in result.items() if c}
    assert min(result) == 30
    return {j-30: c for j, c in result.items()}


def compute():
    dimension, determinant_degree = 4096, 56*4096
    tail_count = dimension-(1+14+35)
    tail_mean = Fraction(determinant_degree, tail_count)
    assert tail_count == 4046 and tail_mean < 57
    rho = Fraction(5, 7)
    assert 2*rho*rho > 1 and 4*rho == Fraction(20, 7)
    assert 8/rho == Fraction(56, 5) and 4*rho*(8/rho) == 32

    reference = json.loads((ROOT/'evidence/spatial-mode-colour-fusion.json').read_text())
    for axis in reference['axes']:
        prefix = axis['whole_prefix']
        assert prefix['0'] >= 1 and prefix['6'] >= 14 and prefix['12'] >= 35
        assert axis['pressure']['dimension'] == dimension
        assert axis['pressure']['determinant_degree'] == determinant_degree

    coefficients = polynomial()
    assert max(coefficients) == 1965 and coefficients[0] == 98
    # Integer-scaled Bernstein coefficients, with exact signs throughout.
    denominator = lcm(*(comb(DEGREE, j) for j in coefficients))
    bernstein = [0]*(DEGREE+1)
    for j, coefficient in coefficients.items():
        weight = coefficient*(denominator//comb(DEGREE, j))
        value = 1
        for k in range(j, DEGREE+1):
            bernstein[k] += weight*value
            if k < DEGREE:
                value = value*(k+1)//(k+1-j)
    assert min(bernstein) > 0

    # Independent monomial reconstruction by forward differences, including
    # every zero coefficient above the actual polynomial degree.
    row = bernstein[:]
    for j in range(DEGREE+1):
        value = comb(DEGREE, j)*row[0]
        assert value % denominator == 0
        assert value//denominator == coefficients.get(j, 0)
        row = [right-left for left, right in zip(row, row[1:])]

    minimum = Fraction(min(bernstein), denominator)
    names = ['tools/check_q8_leading_prefix_floor.py',
             'tools/check_spatial_mode_colour_fusion.py',
             'evidence/spatial-mode-colour-fusion.json']
    return dict(
        scope='Ordinary protected-prefix profile-envelope barrier at q8, with an exact scalar certificate. No new tensor profile, independent mathematical review, factoring lower bound or stronger complete factoring exponent.',
        dimension=dimension, determinant_degree=determinant_degree,
        protected_orders=[[0, 1], [6, 14], [12, 35]],
        protected_order_interpretation='Distinct selected orders bounded ABOVE by these values; exact multiplicities are not required.',
        tail_count=tail_count, tail_mean_upper=str(tail_mean), tail_exponent=57,
        pressure_floor=str(rho), target='1/sqrt(2)',
        paired_word_degree=32, entropy_partition='(1+7*t^(56/5))^2',
        polynomial=sorted(coefficients.items()), actual_polynomial_degree=max(coefficients),
        bernstein_degree=DEGREE, positive_coefficients=len(bernstein),
        monomial_reconstructions=DEGREE+1, minimum_index=bernstein.index(min(bernstein)),
        minimum=str(minimum),
        integer_coefficient_sha256=sha256(json.dumps(bernstein, separators=(',', ':')).encode()).hexdigest(),
        source_sha256={name: sha256((ROOT/name).read_bytes()).hexdigest() for name in names})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(compute()))
    if args.write:
        EXPECTED.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert result == json.loads(EXPECTED.read_text())
    print('Q8 leading prefix: 4097 positive Bernstein coefficients and reconstructions; 5/7 floor PASS')
