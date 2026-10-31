#!/usr/bin/env python3
"""Exact scalar certificate and algebraic controls for the global binary floor.

The finite certificate uses outward integer intervals only.  The global
extension is the analytic argument recorded in the accompanying theorem.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from charged_gaussian_floor import Dyadic


def scalar_kernel_cover():
    """Certify f(u)+f(3u/2)-sqrt(2)f(3u/(2sqrt(2))) on [1/4,20]."""
    arithmetic = Dyadic(80)
    pending = [('', Fraction(1, 4), Fraction(20))]
    cover = []
    evaluated = 0
    while pending:
        path, left, right = pending.pop()
        argument = (arithmetic.rational(left.numerator, left.denominator)[0],
                    arithmetic.rational(right.numerator, right.denominator)[1])
        scaled = arithmetic.multiply(arithmetic.rational(3, 2), argument)
        positive = arithmetic.add(arithmetic.softplus(argument),
                                  arithmetic.softplus(scaled))
        negative = arithmetic.multiply(
            arithmetic.root_two,
            arithmetic.softplus(arithmetic.divide(scaled, arithmetic.root_two)))
        lower = arithmetic.subtract(positive, negative)[0]
        evaluated += 1
        if lower > 0:
            cover.append([path, str(lower)])
        else:
            assert len(path) < 24, (path, str(left), str(right), lower)
            middle = (left + right) / 2
            pending.extend(((path + '1', middle, right),
                            (path + '0', left, middle)))

    # Analytic tails: h(0)>8/21, h'>=-5/4; and the large-u
    # negative/positive ratio is <(15/7)exp(-u/20).
    assert Fraction(8, 21) - Fraction(5, 16) > 0
    assert Fraction(1) + 1 + Fraction(1, 2) + Fraction(1, 6) > Fraction(15, 7)
    assert Fraction(10, 7) ** 2 > 2
    assert Fraction(3, 2) / Fraction(10, 7) == Fraction(21, 20)
    assert evaluated == 2 * len(cover) - 1
    return {
        'status': 'Complete outward-integer positive cover; analytic tails are proved separately.',
        'kernel': 'f(u)+f(3*u/2)-sqrt(2)*f(3*u/(2*sqrt(2)))',
        'domain': ['1/4', '20'], 'dyadic_precision': 80,
        'intervals_evaluated': evaluated,
        'positive_cover_intervals': len(cover),
        'minimum_lower_numerator': str(min(int(value) for _, value in cover)),
        'common_denominator': str(arithmetic.scale), 'cover': cover,
    }


def exact_extension_controls():
    folding = edge = 0
    # The continuous folding proof requires M>=2.  At the right endpoint
    # use delta=2, not the floor representation with M=1.
    for n in range(2, 129):
        for numerator in range(0, 4 * (2 * n - 2) + 1):
            charge = Fraction(numerator, 4)
            k = min(charge // 2, n - 2)
            delta = charge - 2 * k
            lower, upper = n - k, n + k
            assert 0 <= delta <= 2 and 2 <= lower <= upper
            # For lengths >=4 the proved wider band contains [0,2]:
            # 2m-2-2sqrt(2)>2, verified by a rational upper bound.
            for length in (lower, upper):
                assert length in (2, 3) or 2 * length - 2 - 2 * Fraction(3, 2) >= 2
            folding += 1
    for i in range(2, 257):
        b = 2 * i + 1
        for a in (Fraction(-2), Fraction(-7, 4), Fraction(-1), Fraction(-1, 4), Fraction(0)):
            alpha = 2 + a / b
            assert alpha >= Fraction(8, 5) > Fraction(3, 2)
            assert 2 * b + a - 1 > 0
            # Convex averaging arguments and the argument of the negative
            # term agree exactly before any transcendental evaluation.
            assert ((2 * b + a - 1) + (2 * b + a + 1)) / 2 == 2 * b + a
            assert alpha * b == 2 * b + a
            edge += 1
    digit_identities = 0
    coefficients = [1]
    for exponent in range(10):
        shift = 1 << exponent
        expanded = [0] * (len(coefficients) + shift)
        for degree, coefficient in enumerate(coefficients):
            expanded[degree] += coefficient
            expanded[degree + shift] += coefficient
        coefficients = expanded
        assert coefficients == [1] * (1 << (exponent + 1))
        digit_identities += 1
    return {'folding_base_range_controls': folding,
            'edge_scalar_reduction_controls': edge,
            'exact_power_of_two_digit_factorizations': digit_identities}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=root / 'evidence' / 'global-charged-gaussian-floor.json')
    args = parser.parse_args()
    result = {
        'status': 'PASS: exact scalar certificate and global-extension hypotheses.',
        'scope': 'Binary interval Gaussian: every integer n>=1, real charge J, and real scale t>0.',
        'scalar_kernel': scalar_kernel_cover(),
        **exact_extension_controls(),
        'source_sha256': {
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (Path(__file__).resolve(), root / 'tools' / 'charged_gaussian_floor.py')
        },
        'dependencies': [
            'evidence/exterior-charged-gaussian-floor.json',
            'evidence/folded-charged-gaussian-floor.json',
            'evidence/small-charged-gaussian-bases.json',
        ],
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    summary = {key: value for key, value in result.items() if key != 'scalar_kernel'}
    summary['scalar_kernel'] = {key: value for key, value in result['scalar_kernel'].items()
                                if key != 'cover'}
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
