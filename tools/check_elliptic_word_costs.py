#!/usr/bin/env python3
"""Finite group controls for the elliptic word multiplication cost proofs."""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/elliptic-word-cost-controls.json'


def product_count(order, value):
    if not value:
        return order*(2+order.bit_length()-1)//2
    valuation = (value & -value).bit_length()-1
    return (valuation+1)*order//2


def run():
    cyclic = []
    for width in range(1, 9):
        order = 1 << width
        histogram = Counter(a*b % order for a, b in product(range(order), repeat=2))
        assert sum(histogram.values()) == order*order
        assert all(histogram[value] == product_count(order, value) for value in range(order))
        assert max(histogram.values()) == sum(gcd(a, order) for a in range(order))
        assert 2*max(histogram.values()) == order*(width+2)
        cyclic.append(dict(width=width, order=order, grid_points=order*order,
                           maximum_agreement=max(histogram.values()),
                           histogram_sha256=sha256(json.dumps(sorted(histogram.items())).encode()).hexdigest()))
    outside = []
    for width in (1, 2, 3):
        order = 1 << width
        for extra in (1, 2):
            ambient, step = order << extra, 1 << extra
            cases, maximum = 0, 0
            for A, B, C in product(range(ambient), repeat=3):
                if B % step == 0:
                    continue
                count = sum((step*a*b-A*a-B*b-C) % ambient == 0
                            for a, b in product(range(order), repeat=2))
                assert count <= order
                maximum = max(maximum, count)
                cases += 1
            outside.append(dict(width=width, ambient_order=ambient,
                                outside_slope_cases=cases, maximum_agreement=maximum,
                                agreement_upper=order))
    # A geometric divisor of separate degree D cannot contain these many
    # entire slices. These controls audit the constant bookkeeping only.
    slice_controls = []
    for degree in (1, 2, 3, 8, 31, 256):
        slopes = list(range(2*degree+1, 8*degree+2, 2))
        assert len(slopes) == 3*degree+1
        assert all(a % 2 and a*a > 3*degree for a in slopes)
        assert slopes[-1] <= 9*degree
        pullback_bound = 3*degree*(1+slopes[-1]**2)
        assert pullback_bound <= 246*degree**3
        affine_slopes = list(range(2*degree+1, 8*degree+4, 2))
        assert len(affine_slopes) == 3*degree+2
        assert affine_slopes[-1] <= 11*degree
        affine_pullback = 3*degree*(1+affine_slopes[-1]**2)
        assert affine_pullback <= 366*degree**3
        slice_controls.append(dict(degree=degree, distinct_slices=len(slopes),
                                   maximum_slope=slopes[-1],
                                   maximum_pullback_degree=pullback_bound,
                                   declared_bound=246*degree**3,
                                   affine_slices=len(affine_slopes),
                                   affine_pullback_degree=affine_pullback,
                                   affine_declared_bound=366*degree**3+2))
    return dict(scope='Finite cyclic-group and integer-constant controls; the geometric all-size statements use the ordinary proofs and cited standard theorems.',
                cyclic_product_controls=cyclic, outside_subgroup_controls=outside,
                divisor_slice_controls=slice_controls,
                source=dict(author='J. S. Milne', title='Abelian Varieties', version='2.00, March 2008',
                            url='https://www.jmilne.org/math/CourseNotes/AVc.pdf',
                            standard_results=['I Theorem 3.2', 'I Corollary 3.7', 'I Theorem 7.2']),
                checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    args = parser.parse_args()
    receipt = run()
    if args.generate:
        RECEIPT.write_text(json.dumps(receipt, indent=2)+'\n')
    else:
        assert receipt == json.loads(RECEIPT.read_text())
    print(json.dumps(dict(scope=receipt['scope'],
                          cyclic_grids=len(receipt['cyclic_product_controls']),
                          outside_slope_cases=sum(r['outside_slope_cases'] for r in receipt['outside_subgroup_controls']),
                          slice_controls=len(receipt['divisor_slice_controls']))))


if __name__ == '__main__':
    main()
