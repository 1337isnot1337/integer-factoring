#!/usr/bin/env python3
"""Exact finite certificates for an all-multiple fused spectator theorem.

Two bivariate polynomial period identities at each fixed root order prove
the infinite family. Complete larger matrices provide independent finite
controls; no extrapolation from those controls is used in the proof.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from compressed_fused_six_vertex import construct_compressed
from explicit_fused_row_formula import explicit_row
from f16_identity_operator_profiles import MUL, INV, evaluate_one
from fused_six_vertex_crossings import power, pscale


PARAMETERS = {2: 3, 8: 5, 4: 15}


def add(first, second):
    result = dict(first)
    for monomial, value in second.items():
        result[monomial] = result.get(monomial, 0) ^ value
        if result[monomial] == 0:
            del result[monomial]
    return result


def multiply(first, second):
    result = {}
    for (x1, z1), a in first.items():
        for (x2, z2), b in second.items():
            monomial = x1 + x2, z1 + z2
            result[monomial] = result.get(monomial, 0) ^ MUL[a][b]
            if result[monomial] == 0:
                del result[monomial]
    return result


def encoded(polynomial):
    return [[x, z, coefficient] for (x, z), coefficient in sorted(polynomial.items())]


def digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def transition(theta, ell, index, output, middle):
    """Literal nonzero column entries of the displayed T_i^e polynomial."""
    b, i = middle, index
    if output == 0:
        diagonal = {(0, 0): power(theta, b - 1),
                    (0, 1): power(theta, 1 - b + 2 * i)}
        coefficient = MUL[power(theta, i + 2)][power(theta, -1) ^ power(theta, -2 * b - 3)]
        return [(b, diagonal)] + ([(b + 1, {(1, 0): coefficient})] if coefficient else [])
    diagonal = {(1, 0): power(theta, i - b),
                (1, 1): power(theta, 3 * i + b + 4)}
    coefficient = MUL[power(theta, 2 * i)][theta ^ power(theta, 2 * b + 1)]
    return [(b, diagonal)] + ([(b - 1, {(0, 1): coefficient})] if coefficient else [])


def period(theta, ell, output, shift=0):
    matrix = [[{(0, 0): 1} if i == j else {} for j in range(ell)] for i in range(ell)]
    for index in range(shift, shift + ell):
        updated = [[{} for _ in range(ell)] for _ in range(ell)]
        for row, column in product(range(ell), repeat=2):
            for inner, coefficient in transition(theta, ell, index, output, column):
                assert 0 <= inner < ell
                updated[row][column] = add(updated[row][column], multiply(matrix[row][inner], coefficient))
        matrix = updated
    target = {(ell * output, 0): 1, (ell * output, ell): 1}
    assert all(value == (target if i == j else {})
               for i, row in enumerate(matrix) for j, value in enumerate(row))
    return matrix


def row_interpretation(theta, ell):
    checked = 0
    m = ell - 1
    for index in range(ell):
        literal = explicit_row(m, theta, index)
        # explicit_row makes the scalar gcd monic and moves its old leading
        # scalar onto L_m. Undo that field scalar to obtain unscaled L_m.
        leading = 1
        for j in range(1, m):
            leading = MUL[leading][power(theta, 1 + 2 * index - 2 * j)]
        for output, middle in product(range(2), range(ell)):
            expected = {}
            for input_bit in range(2):
                new_middle = middle + input_bit - output
                polynomial = literal['matrix'].get((new_middle, output, input_bit, middle), {})
                polynomial = pscale(polynomial, INV[leading])
                coefficient = power(theta, (index + 2) * input_bit)
                value = {(input_bit, z): MUL[coefficient][a] for z, a in polynomial.items()}
                if value:
                    expected[new_middle] = value
            actual = dict(transition(theta, ell, index, output, middle))
            assert actual == expected
            checked += 1
    return checked


def generating_vector_controls(theta, ell):
    checked = 0
    for m in range(1, 9):
        for word in product(range(2), repeat=m):
            weight = sum(word)
            inversions = sum(word[i] and not word[j] for i in range(m) for j in range(i + 1, m))
            exponent = weight * (weight - 1) // 2 + sum(i - (m - 1) for i, bit in enumerate(word) if bit)
            assert power(theta, exponent) == power(theta, -inversions)
            checked += 1
    for multiple in (1, 2, 3, 7, 16):
        for remainder in range(ell):
            weight = multiple * ell + remainder
            assert power(theta, weight * (weight - 1) // 2) == power(theta, remainder * (remainder - 1) // 2)
            checked += 1
    return checked


def complete_multiple_control(theta, multiple, base):
    ell = PARAMETERS[theta]
    m, q = multiple * ell - 1, multiple * ell
    larger = construct_compressed(m, theta, normalize_joint=True)
    assert larger['endpoint_identity'] and larger['degree'] == base['degree'] == ell - 1
    matrix, endpoint = larger['matrix'], larger['endpoint_before_normalization']
    scalar = endpoint[0][0]
    assert scalar
    normalized_entries = raw_endpoint_entries = 0
    for c, d, a, b in product(range(q), repeat=4):
        fine_key = (ell * (c % ell) + d % ell, ell * (a % ell) + b % ell)
        expected = (base['matrix'][fine_key[0]][fine_key[1]]
                    if c // ell == a // ell and d // ell == b // ell else {})
        assert matrix[q * c + d][q * a + b] == expected
        normalized_entries += 1
        expected_endpoint = scalar if (c % ell == a % ell and d % ell == b % ell
                                       and c // ell == b // ell and d // ell == a // ell) else 0
        assert endpoint[q * c + d][q * a + b] == expected_endpoint
        raw_endpoint_entries += 1
    return {'theta': theta, 'root_order': ell, 'multiple': multiple, 'm': m, 'alphabet': q,
            'degree': larger['degree'], 'normalized_matrix_entries': normalized_entries,
            'raw_endpoint_entries': raw_endpoint_entries, 'raw_endpoint_scalar': scalar,
            'normalized_matrix_sha256': digest(matrix), 'raw_endpoint_sha256': digest(endpoint),
            'implied_smith_multiplicity': multiple ** 2}


def parse_multiple(value):
    try:
        theta, multiple = map(int, value.split(':'))
    except ValueError as exc:
        raise argparse.ArgumentTypeError('Use theta:positive-multiple.') from exc
    if theta not in PARAMETERS or multiple < 1:
        raise argparse.ArgumentTypeError('Theta must be 2, 8, or 4, and the multiple positive.')
    return theta, multiple


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--literal-multiple', type=parse_multiple, action='append', default=[],
                        help='add a complete compressed-matrix control theta:t')
    args = parser.parse_args()
    periods, bases = [], {}
    totals = {'certified_period_identities': 0, 'extra_cyclic_shift_controls': 0,
              'literal_row_interpretations': 0, 'generating_vector_coefficients': 0,
              'base_endpoint_entries': 0}
    base_records = []
    for theta, ell in PARAMETERS.items():
        assert power(theta, ell) == 1 and all(power(theta, i) != 1 for i in range(1, ell))
        totals['literal_row_interpretations'] += row_interpretation(theta, ell)
        totals['generating_vector_coefficients'] += generating_vector_controls(theta, ell)
        for output in range(2):
            matrix = period(theta, ell, output)
            periods.append({'theta': theta, 'root_order': ell, 'output_bit': output,
                            'matrix': [[encoded(value) for value in row] for row in matrix],
                            'off_diagonal_zero': True,
                            'diagonal_terms': encoded(matrix[0][0])})
            totals['certified_period_identities'] += 1
            for shift in sorted({1, ell - 1}):
                assert period(theta, ell, output, shift) == matrix
                totals['extra_cyclic_shift_controls'] += 1
        base = construct_compressed(ell - 1, theta, normalize_joint=True)
        assert base['endpoint_identity'] and base['degree'] == ell - 1
        assert base['scalar_identity_endpoint_before_normalization']
        assert all(evaluate_one(value) == int(i == j)
                   for i, row in enumerate(base['matrix']) for j, value in enumerate(row))
        bases[theta] = base
        totals['base_endpoint_entries'] += ell ** 4
        base_records.append({'theta': theta, 'root_order': ell, 'm': ell - 1,
                             'degree': base['degree'],
                             'normalized_matrix_sha256': digest(base['matrix']),
                             'endpoint_sha256': digest(base['endpoint_before_normalization']),
                             'strict_endpoint_entries': ell ** 4})
    requested = sorted(set([(2, 2), (2, 3), (8, 2)] + args.literal_multiple))
    controls = [complete_multiple_control(theta, multiple, bases[theta]) for theta, multiple in requested]
    totals['complete_multiple_cases'] = len(controls)
    totals['normalized_multiple_entries'] = sum(row['normalized_matrix_entries'] for row in controls)
    totals['raw_multiple_endpoint_entries'] = sum(row['raw_endpoint_entries'] for row in controls)
    directory = Path(__file__).parent
    sources = (Path(__file__).name, 'compressed_fused_six_vertex.py', 'explicit_fused_row_formula.py',
               'fused_six_vertex_crossings.py', 'f16_identity_operator_profiles.py', 'f4_lattice_linear.py')
    result = {'schema_version': 1, 'status': 'pass',
              'claim': 'Six exact transfer-period identities establish all-t spectator factorization for m=t*ell-1 at the three specified F16 parameters; no universal weighted barrier or factoring improvement claimed.',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in sources},
              'period_identities': periods, 'base_crossings': base_records,
              'complete_multiple_controls': controls, 'totals': totals}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(totals, sort_keys=True))


if __name__ == '__main__':
    main()
