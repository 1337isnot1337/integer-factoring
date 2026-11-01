#!/usr/bin/env python3
"""Exact formal identities proving all-degree row fusion and pair symmetry."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from compressed_fused_six_vertex import rectangular_row
from explicit_fused_row_formula import explicit_row


# A Laurent polynomial over F2 is its finite set of nonzero monomials.
# Exponents are ordered (theta, N, B, z), with N and B independent variables.
def monomial(theta=0, n=0, b=0, z=0):
    return {(theta, n, b, z)}


def add(*polynomials):
    result = set()
    for polynomial in polynomials:
        result.symmetric_difference_update(polynomial)
    return result


def multiply(*polynomials):
    result = monomial()
    for polynomial in polynomials:
        updated = set()
        for first in result:
            for second in polynomial:
                term = tuple(a + b for a, b in zip(first, second))
                if term in updated:
                    updated.remove(term)
                else:
                    updated.add(term)
        result = updated
    return result


def shift_weight(polynomial, n_shift=0, b_shift=0, z_shift=0):
    return {(q + n * n_shift + b * b_shift + z * z_shift, n, b, z)
            for q, n, b, z in polynomial}


def row_entry(input_bit, output_bit):
    if (input_bit, output_bit) == (0, 0):
        return add(monomial(-1, 0, 1), monomial(1, 0, -1, 1))
    if (input_bit, output_bit) == (0, 1):
        return add(monomial(1, 0, 0, 1), monomial(-1, -2, 2, 1))
    if (input_bit, output_bit) == (1, 0):
        return add(monomial(0, 1), monomial(-2, 1, -2))
    return add(monomial(-1, 1, -1), monomial(1, -1, 1, 1))


def final_gate(first, second, output_first, output_second):
    if first + second != output_first + output_second:
        return set()
    if first == second:
        return add(monomial(-1), monomial(3, -2, 0, 1))
    if output_first == first:
        c = add(monomial(1), monomial(-1))
        return multiply(c, monomial(2, -2, 0, 1)) if first == 0 else c
    return add(monomial(), monomial(2, -2, 0, 1))


def formal_row_checks():
    result = []
    for a in range(2):
        for d in range(2):
            for e in range(2):
                source = set()
                for input_last in range(2):
                    carrier = d + e - input_last
                    if not 0 <= carrier <= 1:
                        continue
                    split = monomial(0, 0, -1) if input_last == 0 else monomial()
                    previous = shift_weight(row_entry(a, carrier), n_shift=-1, b_shift=-input_last)
                    source = add(source, multiply(split, previous,
                                                   final_gate(carrier, input_last, e, d)))
                new_factor = add(monomial(-1), monomial(3, -2, 0, 1))
                split = monomial(d - a, 0, -1) if e == 0 else monomial()
                target = multiply(new_factor, split, row_entry(a, d))
                assert source == target
                result.append({'input_bit': a, 'output_bit': d, 'last_symmetric_bit': e,
                               'source_monomials': sorted(source), 'target_monomials': sorted(target),
                               'residual_monomials': sorted(add(source, target))})
    return result


def formal_pair_checks():
    result = []
    for input_weight in range(3):
        coefficients = []
        for output_low, output_high in ((0, 1), (1, 0)):
            value = set()
            for input_low in range(2):
                input_high = input_weight - input_low
                if not 0 <= input_high <= 1:
                    continue
                split = monomial(-1) if (input_low, input_high) == (1, 0) else monomial()
                first = shift_weight(row_entry(input_high, output_high), z_shift=2)
                second = shift_weight(row_entry(input_low, output_low), b_shift=input_high - output_high)
                value = add(value, multiply(split, first, second))
            coefficients.append(value)
        target = multiply(monomial(-1), coefficients[0])
        assert coefficients[1] == target
        result.append({'input_weight': input_weight,
                       'coefficient_01_monomials': sorted(coefficients[0]),
                       'coefficient_10_monomials': sorted(coefficients[1]),
                       'theta_inverse_times_01_monomials': sorted(target),
                       'residual_monomials': sorted(add(coefficients[1], target))})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    row_identities = formal_row_checks()
    pair_identities = formal_pair_checks()
    finite = []
    for theta in (2, 8, 4):
        for m in range(1, 19):
            for shift in sorted({0, m - 1}):
                explicit = explicit_row(m, theta, shift)
                recursive = rectangular_row(m, theta, shift)
                assert explicit['matrix'] == recursive['matrix']
                assert explicit['common_factor'] == recursive['common_factor']
                finite.append({'m': m, 'theta': theta, 'shift': shift,
                               'admissible_matrix_entries': len(explicit['matrix']),
                               'common_degree': max(explicit['common_factor'])})
    directory = Path(__file__).parent
    sources = (Path(__file__).name, 'explicit_fused_row_formula.py', 'compressed_fused_six_vertex.py',
               'fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py', 'f4_lattice_linear.py')
    result = {
        'scope': 'Eight formal Laurent induction identities and three adjacent-carrier symmetry identities over F2; finite endpoint regularity is a separate condition.',
        'ring': 'F2[theta^±1,N^±1,B^±1,z], with independent N and B before specialization.',
        'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in sources},
        'row_identities': row_identities, 'pair_identities': pair_identities,
        'finite_comparisons': finite,
        'totals': {'formal_row_identities': len(row_identities),
                   'formal_pair_identities': len(pair_identities),
                   'complete_finite_row_comparisons': len(finite),
                   'finite_matrix_entries': sum(row['admissible_matrix_entries'] for row in finite)}}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['totals']))


if __name__ == '__main__':
    main()
