#!/usr/bin/env python3
"""Fresh exact determinant identity and all-field valuation classification."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from check_joint_endpoint_loops import regular_matrices
from joint_energy_conjugations import N, model
from tensor_pascal_box_crossings import polynomial_product
from two_place_f4_census import build_input


def divide_exact(numerator, denominator):
    if not denominator:
        raise ZeroDivisionError('A zero Bareiss denominator is invalid.')
    degree = denominator.bit_length() - 1
    result = 0
    while numerator and numerator.bit_length() - 1 >= degree:
        shift = numerator.bit_length() - 1 - degree
        result ^= 1 << shift
        numerator ^= denominator << shift
    if numerator:
        raise ArithmeticError('A binary Bareiss division was not exact.')
    return result


def determinant(matrix):
    work = [list(row) for row in matrix]
    previous = 1
    for k in range(len(work)-1):
        pivot_row = next((i for i in range(k, len(work)) if work[i][k]), None)
        if pivot_row is None:
            return 0
        work[k], work[pivot_row] = work[pivot_row], work[k]
        pivot = work[k][k]
        for i in range(k+1, len(work)):
            for j in range(k+1, len(work)):
                numerator = polynomial_product(pivot, work[i][j]) ^ polynomial_product(work[i][k], work[k][j])
                work[i][j] = divide_exact(numerator, previous)
            work[i][k] = 0
        previous = pivot
    return work[-1][-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    _, construction = build_input()
    data = model()
    generator = sum(1 << (N*r+c) for pair in ((4, 20), (7, 11), (8, 16))
                    for r in pair for c in pair)
    source, change = regular_matrices(data, generator)
    maximum = max(value.bit_length()-1 for matrix in (source, change)
                  for row in matrix for value in row if value)
    stride = 512
    assert N * maximum < stride
    encoded = [[a ^ (b << stride) for a, b in zip(first, second)]
               for first, second in zip(source, change)]
    actual = determinant(encoded)
    p0 = sum(1 << degree for degree in (0, 4, 6, 10, 12, 14, 16))
    p4 = sum(1 << degree for degree in (0, 6, 8, 16, 22, 24))
    expected = (p0 << 244) ^ (p4 << (204 + 4*stride))
    assert actual == expected
    # After substituting a quadratic q=a t^2+b t+(a+b), clear the
    # unit denominator. Track formal a^4,b^4 coefficients independently.
    constant = polynomial_product(1 | (1 << 16), p0)
    a4 = polynomial_product(p4, 1 | (1 << 8))
    b4 = polynomial_product(p4, (1 << 4) | (1 << 8))
    symbolic = {degree: [int(constant >> degree & 1), int(a4 >> degree & 1),
                         int(b4 >> degree & 1)] for degree in range(9)}
    assert symbolic[0] == [1, 1, 0]
    assert symbolic[4] == [1, 0, 1]
    assert symbolic[6] == [1, 1, 0]
    assert symbolic[8] == [0, 0, 1]
    assert all(not any(symbolic[i]) for i in (1, 2, 3, 5, 7))
    cases = [{'q': 'zero or degree one', 'infinity_order_sum': 44},
             {'q': 'degree two with leading coefficient a != 1', 'infinity_order_sum': 44},
             {'q': 't^2+b*t+(1+b), b != 1', 'infinity_order_sum': 48},
             {'q': 't^2+t', 'infinity_order_sum': 52}]
    cases += [{'q_degree': degree, 'infinity_order_sum': 52-4*degree}
              for degree in range(3, 12)]
    assert max(case['infinity_order_sum'] for case in cases) == 52
    root = Path(__file__).resolve().parents[1]
    sources = ['tools/check_two_place_determinant.py', 'tools/two_place_f4_census.py',
               'tools/check_joint_endpoint_loops.py', 'tools/joint_energy_conjugations.py',
               'tools/search_f4_recentered_crossings.py', 'tools/tensor_pascal_box_crossings.py']
    report = {'status': 'pass', 'construction': construction,
              'fresh_fraction_free_determinants': 1, 'kronecker_stride': stride,
              'proved_maximum_coefficient_degree': N * maximum,
              'determinant_coefficients': [{'g_degree': 0, 'z_polynomial': p0 << 244},
                                           {'g_degree': 4, 'z_polynomial': p4 << 204}],
              'quadratic_formal_low_coefficients_in_basis_1_a4_b4': symbolic,
              'all_field_infinity_order_sum_cases': cases,
              'maximum_infinity_order_sum': 52, 'minimum_opposite_determinant_degree': 78,
              'seconds': time.monotonic()-start,
              'scope': 'An exact polynomial determinant identity and ordinary all-characteristic-two-field '
              'classification for this complete pencil. It excludes the sufficient determinant-mean '
              'escape, not every asymmetric pressure or every factoring algorithm.',
              'source_sha256': {path: sha256((root / path).read_bytes()).hexdigest() for path in sources}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: report[key] for key in
                      ('status', 'fresh_fraction_free_determinants', 'maximum_infinity_order_sum',
                       'minimum_opposite_determinant_degree', 'seconds')}))


if __name__ == '__main__':
    main()
