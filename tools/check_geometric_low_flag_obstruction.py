#!/usr/bin/env python3
"""Replay literal geometric flag identities and independent finite-field controls."""

from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_pure_two_copy_flag_floor import (
    in_span, integer_pressure_certificate, profile_and_flags, span_basis,
)


def specialize_equation(form, first_left, first_right):
    result = set()
    for coordinate in range(25):
        if not form >> coordinate & 1:
            continue
        i, j = divmod(coordinate, 5)
        if i < first_left or j < first_right:
            continue
        powers = [0] * 10
        if i != first_left:
            powers[i] = 1
        if j != first_right:
            powers[5 + j] = 1
        monomial = tuple(powers)
        if monomial in result:
            result.remove(monomial)
        else:
            result.add(monomial)
    return result


def identity_residual(chart, forms):
    p, q = chart['first_left'], chart['first_right']
    expected_variables = [['a', i] for i in range(p + 1, 5)] + [['b', j] for j in range(q + 1, 5)]
    assert chart['variables'] == expected_variables
    positions = [index + (5 if letter == 'b' else 0) for letter, index in expected_variables]
    assert len(chart['multipliers']) == len(chart['equations']) == len(forms)
    residual = set()
    for form, claimed_equation, multiplier in zip(forms, chart['equations'], chart['multipliers']):
        equation = specialize_equation(form, p, q)
        assert sorted(tuple(term[index] for index in positions) for term in equation) == list(map(tuple, claimed_equation))
        for term in multiplier:
            assert len(term) == len(positions) and all(isinstance(power, int) and power >= 0 for power in term)
            full = [0] * 10
            for index, power in zip(positions, term):
                full[index] = power
            for monomial in equation:
                at = tuple(a + b for a, b in zip(full, monomial))
                assert sum(at) <= chart['identity_degree_bound']
                if at in residual:
                    residual.remove(at)
                else:
                    residual.add(at)
    return residual


def finite_multiply(a, b, modulus):
    value = 0
    while b:
        if b & 1:
            value ^= a
        a <<= 1
        b >>= 1
    while value.bit_length() >= modulus.bit_length():
        value ^= modulus << (value.bit_length() - modulus.bit_length())
    return value


def finite_control(forms, size, modulus):
    table = [[finite_multiply(a, b, modulus) for b in range(size)] for a in range(size)]
    inverses = [0] + [next(b for b in range(1, size) if table[a][b] == 1) for a in range(1, size)]
    supports = [[[i for i in range(5) if form >> (5 * i + j) & 1] for j in range(5)] for form in forms]
    checked = 0
    for first in range(5):
        for tail in product(range(size), repeat=4 - first):
            left = (0,) * first + (1,) + tail
            pivots = {}
            for support in supports:
                row = []
                for entries in support:
                    coefficient = 0
                    for i in entries:
                        coefficient ^= left[i]
                    row.append(coefficient)
                for column in range(5):
                    coefficient = row[column]
                    if not coefficient:
                        continue
                    if column not in pivots:
                        pivots[column] = [table[inverses[coefficient]][entry] for entry in row]
                        break
                    row = [entry ^ table[coefficient][previous] for entry, previous in zip(row, pivots[column])]
                if len(pivots) == 5:
                    break
            assert len(pivots) == 5
            checked += 1
    assert checked == (size ** 5 - 1) // (size - 1)
    return {'field_size': size, 'polynomial_modulus': modulus,
            'projective_left_factors_exhausted': checked, 'every_right_constraint_rank': 5}


def verify():
    path = Path('evidence/geometric-low-flag-obstruction.json')
    record = json.loads(path.read_text())
    data = profile_and_flags()
    assert record['reference_digest'] == data['reference_digest']
    basis, forms = record['physical_image_flag_basis'], record['annihilating_forms']
    assert basis == data['common_flag_basis']
    assert len(span_basis(basis)) == 11 and len(span_basis(forms)) == 14
    assert all((row & form).bit_count() % 2 == 0 for row in basis for form in forms)
    assert [(chart['first_left'], chart['first_right']) for chart in record['charts']] == list(product(range(5), repeat=2))
    for chart in record['charts']:
        assert identity_residual(chart, forms) == {(0,) * 10}
    assert record['projective_chart_count'] == len(record['charts']) == 25
    assert record['maximum_identity_degree'] == max(chart['identity_degree_bound'] for chart in record['charts']) == 3
    terms = sum(sum(map(len, chart['multipliers'])) for chart in record['charts'])
    assert terms == record['total_multiplier_monomials'] == 60
    damaged = json.loads(json.dumps(record['charts'][0]))
    damaged['multipliers'][0] = []
    assert identity_residual(damaged, forms) != {(0,) * 10}

    # The next flag really contains a pure matrix, already over the prime field.
    left, right = (0, 1, 1, 1, 0), (0, 0, 1, 0, 0)
    pure = sum((a * b) << (5 * i + j) for i, a in enumerate(left) for j, b in enumerate(right))
    assert in_span(pure, span_basis(data['full_image_basis'][:12]))
    assert not in_span(pure, span_basis(basis))
    pressure = integer_pressure_certificate(data['relative_profile_polynomial'])
    finite = [finite_control(forms, size, modulus) for size, modulus in ((2, 3), (4, 7), (16, 19))]
    return {'status': 'exact rank-one exclusion over every characteristic-two field certified',
            'certificate_sha256': sha256(path.read_bytes()).hexdigest(),
            'reference_digest': data['reference_digest'],
            'literal_projective_chart_identities': 25,
            'maximum_identity_degree': 3, 'total_multiplier_monomials': terms,
            'coefficient_scope': 'Every characteristic-two field and algebraic closure.',
            'independent_finite_field_controls': finite,
            'damaged_identity_rejected': True,
            'next_flag_rank_one_positive_control': {'relative_order': 10, 'left': left, 'right': right},
            'single_pure_two_copy_rate_floor': pressure['rate_floor'],
            'pressure_interval_count': pressure['interval_count'],
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
