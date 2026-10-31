#!/usr/bin/env python3
"""Exact controls for Walsh twirling, dyadic energies, and scalar decoders."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
from random import Random

from boolean_cauchy_operator import kernel_matrices, smith_orders, value_operator
from search_interval_permanent_decoder import feature_matrix
from walsh_cauchy_permanents import (
    interval_first_three_jets, local_kernel_jets, local_kernel_polynomials,
    minimum_assignment, permanent, principal_energies, principal_permanents,
    principal_permanent_polynomials, reflect_interval_jets, sparse_interval_jets,
    walsh_transform,
)


def order(value, precision):
    return (value & -value).bit_length() - 1 if value else precision


def check():
    root = Path(__file__).resolve().parent.parent
    precision = 128
    modulus = 1 << precision
    energy_cases = assignment_cases = determinant_cases = balanced_cases = interval_cases = 0
    for length in (2, 4, 8):
        a = length.bit_length() - 1
        for center in range(1 - length, length):
            if not center:
                continue
            physical, _ = kernel_matrices(length, center, 1, precision)
            values = principal_permanents(physical, modulus)
            for subset in range(1 << length):
                selected = [x for x in range(length) if subset >> x & 1]
                matrix = [[physical[x][y] for y in selected] for x in selected]
                weights = [[order(value, precision) for value in row] for row in matrix]
                absolute, quadratic = principal_energies(length, center, subset)
                assert minimum_assignment(weights) == absolute
                assignment_cases += 1
                orders = smith_orders(matrix, precision)
                assert all(value is not None for value in orders)
                assert sum(orders) == quadratic
                determinant_cases += 1
                observed = order(values[subset], precision)
                assert observed >= absolute
                if absolute == quadratic:
                    assert observed == absolute
                    balanced_cases += 1
                if selected and selected == list(range(selected[0], selected[-1] + 1)):
                    assert absolute == quadratic == observed
                    interval_cases += 1
                    if center == 1:
                        size = len(selected)
                        assert observed == a + 1 - order(size, precision)
                energy_cases += 1

    twirl_cases = twirl_entries = 0
    for length in (2, 4):
        count = 1 << length
        centers = sorted(set([1 - length, -1, 1, length - 1]))
        for center in centers:
            for point in (0, 1, 2, 3):
                physical, _ = kernel_matrices(length, center, point, 24)
                small_modulus = 1 << 24
                kernel = value_operator(physical, small_modulus)
                lambdas = principal_permanents(physical, small_modulus)
                assert lambdas[0] == 1 and all(not value % 2 for value in lambdas[1:])
                expected_by_difference = [value % small_modulus for value in walsh_transform(lambdas)]
                for row in range(count):
                    for column in range(count):
                        literal = sum(kernel[row ^ translate][column ^ translate]
                                      for translate in range(count)) % small_modulus
                        assert literal == expected_by_difference[row ^ column]
                        assert literal & 1
                        if not point:
                            assert literal == 1
                        twirl_entries += 1
                twirl_cases += 1

    kernel_coefficients = 0
    for length in (2, 4, 8):
        for center in range(1 - length, length):
            if not center:
                continue
            for bits in (5, 9):
                fast = local_kernel_jets(length, center, bits)
                independent = local_kernel_polynomials(length, center, bits)
                for difference, polynomial in independent.items():
                    expected = (polynomial + [0] * bits)[:bits]
                    assert fast[difference] == expected
                    assert not any(polynomial[bits:])
                    for degree, value in enumerate(polynomial):
                        assert order(value, bits) >= min(degree, bits)
                    kernel_coefficients += len(polynomial)

    low_jet_coefficients = sparse_interval_coefficients = 0
    for length in (4, 8):
        bits = 32
        fast = interval_first_three_jets(length, bits)
        selected = [(1 << k) - 1 for k in range(length + 1)]
        for center, rows in fast.items():
            independent = principal_permanent_polynomials(length, center, bits, selected)
            for k, subset in enumerate(selected):
                assert rows[k] == (independent[subset] + [0] * 4)[:4]
                low_jet_coefficients += 4
        for center in (1, 3):
            bits = length.bit_length() + 2
            sparse, _ = sparse_interval_jets(length, center, bits)
            reflected = reflect_interval_jets(length, center, bits, sparse)
            for signed, actual in ((center, sparse), (-center, reflected)):
                independent = principal_permanent_polynomials(length, signed, bits, selected)
                for k, subset in enumerate(selected):
                    assert actual[k] == (independent[subset] + [0] * bits)[:bits]
                    assert not any(independent[subset][bits:])
                    sparse_interval_coefficients += bits

    decoder_data = json.loads((root / 'evidence/interval-permanent-decoders.json').read_text())
    decoder_records = []
    for record in decoder_data['records']:
        length, bits = record['L'], record['precision']
        features, matrix = feature_matrix(length, bits)
        assert [list(feature) for feature in features] == record['features']
        mod = 1 << bits
        if record['certified']:
            feature_index = {feature: index for index, feature in enumerate(features)}
            dual = [(feature_index[(term['center'], term['order'])], term['coefficient'])
                    for term in record['dual']]
            for column in range(length + 1):
                assert sum(value * matrix[row][column] for row, value in dual) % mod == 1 << record['division_loss']
            decoder_records.append({'L': length, 'precision': bits, 'status': 'dual_verified',
                                    'basis_directions_checked': length + 1})
        else:
            collision = [record['collision'].get(str(k), 0) for k in range(length + 1)]
            assert sum(collision) & 1
            assert all(sum(a * b for a, b in zip(row, collision)) % mod == 0 for row in matrix)
            decoder_records.append({'L': length, 'precision': bits, 'status': 'six_query_collision_verified',
                                    'nonzero_coefficients': sum(bool(value) for value in collision)})

    witness = json.loads((root / 'evidence/walsh-permanent-witness.json').read_text())
    length, hidden = witness['L'], witness['hidden_bits']
    coefficients = {int(key): int(value) for key, value in witness['coefficient_nonzero'].items()}
    assert length == 8 and hidden == 112 and sum(coefficients.values()) & 1
    assert all(0 <= value < 1 << hidden for value in coefficients.values())
    full_coefficients = [coefficients.get(subset, 0) for subset in range(1 << length)]
    column = list(map(int, witness['canonical_column']))
    creation = witness['canonical_column_precision']
    assert creation == hidden + length == 120
    assert column == [value % (1 << creation) for value in walsh_transform(full_coefficients)]
    assert column[0] == int(witness['target']) and column[0] & 1
    transformed = walsh_transform(column)
    assert all(value % (1 << length) == 0 for value in transformed)
    lifted_coefficients = [value // (1 << length) for value in transformed]
    assert all((left - right) % (1 << hidden) == 0
               for left, right in zip(lifted_coefficients, full_coefficients))
    witness_records = []
    coefficient_controls = 0
    attained_nonconstant = False
    witness_modulus = 1 << (hidden + 1)
    for center in range(-length, length):
        if not center % 2:
            continue
        polynomials = principal_permanent_polynomials(length, center, hidden + 1, coefficients)
        combined = [0] * (length * (2 * length - 1) + 1)
        for subset, polynomial in polynomials.items():
            assert subset and polynomial[0] == 0
            assert all(value % 2 == 0 for value in polynomial)
            coefficient = coefficients[subset]
            for degree, value in enumerate(polynomial):
                combined[degree] = (combined[degree] + coefficient * value) % witness_modulus
                coefficient_controls += 1
        assert all(value % (1 << hidden) == 0 for value in combined)
        local_minimum = min((order(value, hidden + 1) for value in combined if value), default=hidden + 1)
        if any(value for value in combined[1:]):
            attained_nonconstant = True
        witness_records.append({'center': center, 'minimum_content_or_precision_cutoff': local_minimum})
    assert attained_nonconstant
    assert min(record['minimum_content_or_precision_cutoff'] for record in witness_records) == hidden

    # These independently evaluate the congruences in the all-size proof.
    # The cofactor bound is deliberately tested only at size 2^b-1.
    random = Random(20261031)
    circulant_controls = odd_minor_controls = cauchy_step_controls = 0
    for n in (1, 2, 4, 8, 16):
        bits = n + 3
        mod = 1 << bits
        trials = 8 if n <= 8 else 2
        for _ in range(trials):
            row = [1 + 2 * random.randrange(1 << n) for _ in range(n)]
            circulant = [[row[(i - j) % n] for j in range(n)] for i in range(n)]
            value = permanent(circulant, mod)
            assert (value - factorial(n)) % (1 << n) == 0
            assert order(value, bits) == n - 1
            circulant_controls += 1
            m = n - 1
            odd_minor = [[1 + 2 * random.randrange(1 << n) for _ in range(m)]
                         for _ in range(m)]
            assert permanent(odd_minor, mod) % (1 << (m - m.bit_count())) == 0
            odd_minor_controls += 1
        steps = (2, 4, 6, 10) if n <= 8 else (2, 6)
        points = (-7, 1, 3, 11) if n <= 8 else (3, 11)
        for step in steps:
            for point in points:
                matrix = [[pow(point + step * (i - j), -1, mod)
                           for j in range(n)] for i in range(n)]
                value = permanent(matrix, mod)
                expected = factorial(n) * pow(point, -n, mod)
                assert (value - expected) % (1 << n) == 0
                assert order(value, bits) == n - 1
                cauchy_step_controls += 1

    spaced_polynomial_controls = spaced_polynomial_coefficients = 0
    for packet_length in (2, 4, 8, 16):
        a = packet_length.bit_length() - 1
        selected = set()
        if packet_length <= 8:
            for n in (1, 2, 4):
                for step in (2, 4, 6):
                    for start in range(packet_length):
                        positions = [start + step * i for i in range(n)]
                        if positions[-1] < packet_length:
                            selected.add(sum(1 << x for x in positions))
            centers = range(1 - packet_length, packet_length, 2)
        else:
            selected = {sum(1 << (2 * i) for i in range(8))}
            centers = (-15, -7, -1, 1, 7, 15)
        for center in centers:
            for subset in selected:
                n = subset.bit_count()
                bits = n * (a + 2)
                polynomial = principal_permanent_polynomials(
                    packet_length, center, bits, [subset])[subset]
                expected = [0] * len(polynomial)
                expected[n] = 1 << (bits - 1)
                assert polynomial == expected
                spaced_polynomial_controls += 1
                spaced_polynomial_coefficients += len(polynomial)

    factorial_records = []
    direct_permanent_controls = 0
    for n in list(range(1, 17)) + [32, 64, 128]:
        points = [1, 3, 4 * n + 1] if n <= 16 else [4 * n + 1]
        for point in points:
            bits = 8 * n + 32
            mod = 1 << bits
            squared = [[pow(point + 2 * i - 2 * j, -2, mod) for j in range(n)] for i in range(n)]
            observed = smith_orders(squared, bits)
            expected = [2 * k + 2 * (k - k.bit_count()) + order(k + 1, bits) for k in range(n)]
            assert observed == expected
            ordinary_determinant = sum(2 * k + 2 * (k - k.bit_count()) for k in range(n))
            permanent_order = sum(observed) - ordinary_determinant
            assert permanent_order == n - n.bit_count()
            if n <= 8:
                matrix = [[pow(point + 2 * i - 2 * j, -1, mod) for j in range(n)] for i in range(n)]
                assert order(permanent(matrix, mod), bits) == permanent_order
                direct_permanent_controls += 1
            factorial_records.append({'n': n, 'odd_T': point, 'precision': bits,
                                      'squared_cauchy_orders': observed,
                                      'ordinary_determinant_order': ordinary_determinant,
                                      'permanent_order_by_borchardt': permanent_order})

    assert Fraction(1, 2 * 3) - Fraction(1, 1 * 6) == 0
    sources = ['tools/check_walsh_cauchy_permanents.py', 'tools/walsh_cauchy_permanents.py',
               'tools/search_interval_permanent_decoder.py', 'tools/boolean_cauchy_operator.py',
               'tools/cauchy_packet_polynomials.py', 'research/walsh-principal-permanents.md',
               'evidence/walsh-permanent-witness.json', 'evidence/interval-permanent-decoders.json']
    return {'status': 'pass',
            'proved_formula_controls': {
                'principal_energy_cases': energy_cases,
                'independent_minimum_assignments': assignment_cases,
                'independent_determinant_smith_checks': determinant_cases,
                'balanced_exact_permanent_cases': balanced_cases,
                'contiguous_interval_cases': interval_cases,
                'literal_twirl_cases': twirl_cases, 'literal_twirl_entries': twirl_entries,
                'independent_local_kernel_coefficients': kernel_coefficients,
                'low_jet_formula_coefficients': low_jet_coefficients,
                'sparse_interval_and_reflection_coefficients': sparse_interval_coefficients,
                'power_of_two_circulant_congruences': circulant_controls,
                'special_size_odd_minor_bounds': odd_minor_controls,
                'even_spacing_cauchy_congruences': cauchy_step_controls,
                'even_spacing_complete_polynomials': spaced_polynomial_controls,
                'even_spacing_polynomial_coefficients': spaced_polynomial_coefficients},
            'restricted_decoders': decoder_records,
            'high_height_witness': {'L': length, 'hidden_bits': hidden,
                                   'canonical_creation_precision': creation,
                                   'complete_polynomial_coefficients_checked': coefficient_controls,
                                   'odd_original_contraction': str(column[0]),
                                   'local_records': witness_records},
            'finite_factorial_pattern': {
                                        'status': 'Individual squared-Cauchy Smith orders and '
                                                  'non-power-of-two permanent formula remain '
                                                  'finite observations. The note proves the '
                                                  'power-of-two permanent formula for all sizes.',
                                        'direct_permanent_controls': direct_permanent_controls,
                                        'records': factorial_records},
            'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest() for name in sources},
            'scope': 'All-size twirl, discrepancy, and power-of-two even-spacing permanent '
                     'theorems; exact finite restricted decoders; an explicit high-height '
                     'collision; and separately labeled finite individual Smith patterns. '
                     'No general contraction improvement or new factoring score.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = check()
    except Exception as error:
        if args.output:
            args.output.write_text(json.dumps({'status': 'fail', 'error': repr(error)}, indent=2) + '\n')
        raise
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'],
                      'proved_formula_controls': result['proved_formula_controls'],
                      'restricted_decoders': result['restricted_decoders'],
                      'witness_hidden_bits': result['high_height_witness']['hidden_bits'],
                      'factorial_cases': len(result['finite_factorial_pattern']['records'])}, indent=2))


if __name__ == '__main__':
    main()
