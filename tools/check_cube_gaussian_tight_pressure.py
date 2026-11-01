#!/usr/bin/env python3
"""A sharper exact rational pressure certificate for the published cube kernel."""
from fractions import Fraction
from hashlib import sha256
import argparse
import json
from math import log, log2, prod
from pathlib import Path
import time

from cube_homogeneous_gaussian import augmentation_crossblock, exact_monomial_smith, kernel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='evidence/cube-gaussian-tight-pressure.json')
    args = parser.parse_args()
    started = time.monotonic()
    source = kernel()
    profiles = [exact_monomial_smith(augmentation_crossblock(source, axis), axis)['orders']
                for axis in range(3)]
    assert profiles == [[2, 4, 8, 10]] * 3
    p, r, rho = Fraction(5, 24), Fraction(23, 53), Fraction(278, 445)
    assert rho < Fraction(5, 8)
    numerator = prod(r.denominator ** degree + r.numerator ** degree for degree in (1, 2, 4, 5))
    # Z(r) = numerator / 53^12. Clear all entropy and logarithm denominators.
    left_powers = [(numerator, p.denominator * rho.denominator),
                   (p.numerator, 4 * rho.numerator * p.numerator),
                   (p.denominator - p.numerator, 4 * rho.numerator * (p.denominator - p.numerator))]
    right_powers = [(r.denominator, 12 * p.denominator * rho.denominator - 4 * p.numerator * rho.denominator),
                    (r.numerator, 4 * p.numerator * rho.denominator),
                    (p.denominator, 4 * rho.numerator * p.denominator)]
    left = prod(base ** exponent for base, exponent in left_powers)
    right = prod(base ** exponent for base, exponent in right_powers)
    assert left < right
    entropy = -float(p) * log2(float(p)) - (1 - float(p)) * log2(1 - float(p))
    partition = prod(1 + r ** degree for degree in (1, 2, 4, 5))
    actual = (log2(float(partition)) - 4 * float(p) * log2(float(r))) / (4 * entropy)
    assert actual < float(rho)
    def digest(value):
        return sha256(value.to_bytes((value.bit_length() + 7) // 8, 'big')).hexdigest()
    paths = [Path(__file__), Path(__file__).with_name('cube_homogeneous_gaussian.py')]
    evidence = {'schema': 'cube-gaussian-tight-pressure-v1',
                'physical_frequency': str(p), 'scale_after_deflation': str(r),
                'certified_uniform_rho': str(rho), 'previous_uniform_rho': '5/8',
                'three_common_smith_profiles': profiles,
                'partition_numerator': numerator, 'partition_denominator': r.denominator ** 12,
                'left_integer_powers': left_powers, 'right_integer_powers': right_powers,
                'strict_integer_inequality_verified': True,
                'left_bits': left.bit_length(), 'right_bits': right.bit_length(),
                'left_sha256': digest(left), 'right_sha256': digest(right),
                'numerical_rho_at_selected_parameters': actual,
                'numerical_certified_rho': float(rho),
                'cube_transfer_alpha': log(float(1 / rho)) / log(3),
                'scope': 'A sharper uniform local cube certificate, still below the existing one-half factoring exponent.',
                'new_factoring_score': False,
                'source_sha256': {path.name: sha256(path.read_bytes()).hexdigest() for path in paths},
                'elapsed_seconds': time.monotonic() - started}
    output = Path(args.output); output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps({'rho': str(rho), 'parameter_value': actual,
                      'integer_bits': left.bit_length(), 'seconds': evidence['elapsed_seconds']}))


if __name__ == '__main__':
    main()
