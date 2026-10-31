#!/usr/bin/env python3
"""Exact full polynomial module of the actual five-letter Gaussian complements."""

from hashlib import sha256
import json
from pathlib import Path

from aligned_gaussian_module_directions import aligned_module_record
from check_rectangular_defect_support import Rat, ZERO
from check_source_slab_module_saturation import in_rational_span
from source_slab_module_saturation import bounded_polynomial_space, encode, insert_binary, row_degree


def verify():
    path = Path('evidence/aligned-gaussian-complement-directions-q5.json')
    source = json.loads(path.read_text())
    result = json.loads(Path('evidence/aligned-gaussian-module-directions-q5.json').read_text())
    expected = aligned_module_record(source, 18)
    expected['source_record_sha256'] = sha256(path.read_bytes()).hexdigest()
    assert json.loads(json.dumps(expected)) == result
    certificate = result['rational_row_certificate']
    raw = [item['raw_slab_in_z'] for item in source['raw_slabs']]
    rows = [raw[i] for i in certificate['source_indices']]
    assert certificate['rank'] == 37
    rational_checks = 0
    for i, coefficients in enumerate(certificate['transformation']):
        sparse = [(Rat(*coefficient), row) for coefficient, row in zip(coefficients, rows)
                  if coefficient[0]]
        for column in range(125):
            actual = ZERO
            for coefficient, row in sparse:
                actual += coefficient * Rat(row[column])
            assert actual == Rat(*certificate['rref'][i][column])
            rational_checks += 1
    for row in raw:
        assert in_rational_span(row, certificate)
        assert all(value.bit_count() % 2 == 0 for value in row)
    degrees = result['minimal_basis']['degrees']
    assert degrees == [0] + [3] * 8 + [4] * 7 + [5] * 3 + [8] * 5 + [9] * 6 + [10] * 2 + [11] * 2 + [13] * 3
    generators = result['minimal_basis']['generators']
    assert list(map(row_degree, generators)) == degrees
    assert all(in_rational_span(row, certificate) for row in generators)
    leading = {}
    endpoint = {}
    for row, degree in zip(generators, degrees):
        assert insert_binary(leading, sum(((value >> degree) & 1) << j for j, value in enumerate(row)))
        assert insert_binary(endpoint, sum((value.bit_count() % 2) << j for j, value in enumerate(row)))
    assert len(leading) == len(endpoint) == 37

    complete_space = bounded_polynomial_space(certificate, 12)
    shifted = {}
    for row, degree in zip(generators, degrees):
        for shift in range(max(0, 13 - degree)):
            assert insert_binary(shifted, encode([value << shift for value in row], 13))
    assert len(shifted) == complete_space['dimension'] == 239
    assert len(complete_space['equation_pivots']) == 242
    for row in complete_space['slabs']:
        assert not insert_binary(shifted, encode(row, 13))
    controls = {}
    for degree in (0, 3, 5, 8, 12, 13, 17, 20):
        controls[degree] = bounded_polynomial_space(certificate, degree)['dimension']
        assert controls[degree] == sum(max(0, degree - d + 1) for d in degrees)

    tensor_pivots = {}
    for row, primitive in zip(generators, result['primitives']):
        assert primitive['normalized_slab'] == [value ^ (value << 1) for value in row]
        assert all(value.bit_count() % 2 == 0 for value in primitive['normalized_slab'])
    for primitive, first, shift in result['direction_descriptors']:
        row = result['primitives'][primitive]['normalized_slab']
        vector = sum((value << shift) << (19 * (first + 5 * j)) for j, value in enumerate(row))
        assert insert_binary(tensor_pivots, vector)
    assert len(tensor_pivots) == 2120
    endpoint_dimensions = {degree: 5 * sum(max(0, degree - d) for d in degrees)
                           for degree in (1, 4, 5, 6, 9, 10, 11, 12, 14, 18)}

    # The crossing itself must not be inferred from the complementary module.
    excluded_endpoints = 0
    for x in range(5):
        crossing_slab = sum(int(v == x and y == u) << (y + 5 * u + 25 * v)
                            for v in range(5) for u in range(5) for y in range(5))
        assert insert_binary(dict(endpoint), crossing_slab)
        excluded_endpoints += 1
    assert sum(max(0, 14 - d) for d in degrees[:-1]) == 275
    assert controls[13] == 276
    return {'status': 'exact full aligned Gaussian polynomial module certified',
            'source_record_sha256': result['source_record_sha256'],
            'reference_digest': result['reference_digest'],
            'rational_rank': 37, 'rational_reconstruction_entries': rational_checks,
            'minimal_basis_degrees': degrees, 'independent_leading_rows': 37,
            'independent_endpoint_rows': 37,
            'saturation_checked_degree': 12, 'saturation_equation_rank': 242,
            'saturation_polynomial_dimension': 239,
            'bounded_degree_dimensions': controls,
            'complete_endpoint_zero_tensor_dimensions': endpoint_dimensions,
            'literal_degree18_direction_rank': len(tensor_pivots),
            'crossing_endpoint_slabs_outside_complementary_module': excluded_endpoints,
            'missing_generator_negative_controls': 1,
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
