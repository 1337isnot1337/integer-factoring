#!/usr/bin/env python3
"""Exact rational-span, bounded-space, and saturated-module certificates."""

from hashlib import sha256
import json
from pathlib import Path

from check_rectangular_defect_support import Rat, ZERO, ONE
from source_slab_module_saturation import (
    bounded_polynomial_space, direction_record, encode, insert_binary, row_degree,
)


def in_rational_span(row, certificate):
    weights = [Rat(row[column]) for column in certificate['pivot_columns']]
    for column, value in enumerate(row):
        predicted = ZERO
        for weight, source in zip(weights, certificate['rref']):
            numerator, denominator = source[column]
            predicted += weight * Rat(numerator, denominator)
        if predicted != Rat(value):
            return False
    return True


def verify():
    path = Path('evidence/bounded-dual-completion-h16.json')
    source = json.loads(path.read_text())
    data = json.loads(Path('evidence/source-slab-module-directions-h16.json').read_text())
    assert json.loads(json.dumps(direction_record(source, 18))) == data
    certificate = data['rational_row_certificate']
    rows = [source['kernel_slabs'][i] for i in certificate['source_indices']]
    rank = certificate['rank']
    assert rank == len(rows) == 10
    assert certificate['source_indices'] == [0, 3, 5, 7, 8, 11, 15, 16, 21, 23]
    rational_checks = 0
    for i, coefficients in enumerate(certificate['transformation']):
        for column in range(125):
            actual = ZERO
            for coefficient, row in zip(coefficients, rows):
                actual += Rat(*coefficient) * Rat(row[column])
            assert actual == Rat(*certificate['rref'][i][column])
            rational_checks += 1
    for i, row in enumerate(certificate['rref']):
        for j, column in enumerate(certificate['pivot_columns']):
            assert Rat(*row[column]) == (ONE if i == j else ZERO)
    generators = data['minimal_basis']['generators']
    degrees = data['minimal_basis']['degrees']
    assert degrees == [0, 0, 0, 8, 8, 8, 12, 16, 16, 16]
    assert list(map(row_degree, generators)) == degrees
    assert all(in_rational_span(row, certificate) for row in generators)
    assert all(in_rational_span(row, certificate) for row in source['kernel_slabs'])
    leading_pivots = {}
    for row, degree in zip(generators, degrees):
        leading = sum(((value >> degree) & 1) << j for j, value in enumerate(row))
        assert insert_binary(leading_pivots, leading)
    assert len(leading_pivots) == rank

    # Completeness through d-1 certifies the absence of polynomial torsion.
    saturation_degree = max(degrees) - 1
    space = bounded_polynomial_space(certificate, saturation_degree)
    span = {}
    for row, degree in zip(generators, degrees):
        for shift in range(max(0, saturation_degree - degree + 1)):
            assert insert_binary(span, encode([value << shift for value in row], saturation_degree + 1))
    assert len(span) == space['dimension'] == 76
    for row in space['slabs']:
        assert not insert_binary(span, encode(row, saturation_degree + 1))
    assert len(space['equation_pivots']) == 84
    assert data['saturation_certificate']['complete_polynomial_module']

    degree_controls = []
    for degree in (0, 1, 7, 8, 11, 12, 15, 16, 17, 18, 20):
        computed = bounded_polynomial_space(certificate, degree)['dimension']
        formula = sum(max(0, degree - value + 1) for value in degrees)
        assert computed == formula
        degree_controls.append({'degree': degree, 'dimension': computed})

    endpoint_controls = 0
    for row, primitive in zip(generators, data['primitives']):
        assert primitive['normalized_slab'] == [value ^ (value << 1) for value in row]
        for value in primitive['normalized_slab']:
            assert value.bit_count() % 2 == 0
            endpoint_controls += 1
    assert len(data['direction_descriptors']) == 480
    assert len(set(map(tuple, data['direction_descriptors']))) == 480
    for primitive, logical, shift in data['direction_descriptors']:
        assert 0 <= logical < 5 and shift >= 0
        assert data['primitives'][primitive]['normalized_degree'] + shift <= 18
    bounded_endpoint_dimensions = {
        degree: 5 * sum(max(0, degree - value) for value in degrees)
        for degree in (1, 9, 13, 16, 17, 18)}
    assert bounded_endpoint_dimensions == {1: 15, 9: 150, 13: 275, 16: 380, 17: 430, 18: 480}

    # A coordinate outside the pivot set cannot be changed independently.
    bad = generators[0][:]
    outside = next(i for i in range(125) if i not in certificate['pivot_columns'])
    bad[outside] ^= 1
    assert not in_rational_span(bad, certificate)
    # RREF identities over the rational function field are not polynomial bases.
    assert any(denominator != 1 for row in certificate['rref'] for _, denominator in row)
    # Omitting a top-degree row loses a direction at degree sixteen.
    assert sum(max(0, 17 - value) for value in degrees[:-1]) == 85
    assert bounded_polynomial_space(certificate, 16)['dimension'] == 86

    return {'status': 'exact full polynomial-module saturation certificate passed',
            'source_sha256': sha256(path.read_bytes()).hexdigest(),
            'rational_rank': rank, 'rational_reconstruction_entries': rational_checks,
            'minimal_basis_degrees': degrees, 'independent_leading_rows': len(leading_pivots),
            'saturation_degree_bound': saturation_degree, 'saturation_equation_rank': 84,
            'saturation_polynomial_space_dimension': 76,
            'bounded_degree_controls': degree_controls,
            'endpoint_coefficient_controls': endpoint_controls,
            'complete_endpoint_zero_tensor_dimensions': bounded_endpoint_dimensions,
            'negative_controls': 3, 'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
