#!/usr/bin/env python3
"""Reconstruct and certify the exact fourth-root direction family."""

from hashlib import sha256
import json
from pathlib import Path

from check_rectangular_defect_support import Rat, ZERO
from source_slab_module_saturation import bounded_polynomial_space, encode, insert_binary, row_degree
from source_slab_root_directions import compress_exponents, inflate_exponents, root_direction_record


def verify():
    source_path = Path('evidence/source-slab-module-directions-h16.json')
    output_path = Path('evidence/source-slab-fourth-root-directions-h16.json')
    source = json.loads(source_path.read_text())
    result = json.loads(output_path.read_text())
    rebuilt = root_direction_record(source, 4, 18)
    rebuilt['source_record_sha256'] = sha256(source_path.read_bytes()).hexdigest()
    assert json.loads(json.dumps(rebuilt)) == result
    compressed = result['compressed_source_generators']
    coefficient_checks = 0
    for original, root in zip(source['minimal_basis']['generators'], compressed):
        for value, root_value in zip(original, root):
            assert inflate_exponents(root_value, 4) == value
            coefficient_checks += 1

    certificate = result['rational_row_certificate']
    assert certificate['rank'] == 10
    rational_checks = 0
    for i, coefficients in enumerate(certificate['transformation']):
        for column in range(125):
            actual = ZERO
            for coefficient, row in zip(coefficients, compressed):
                actual += Rat(*coefficient) * Rat(row[column])
            assert actual == Rat(*certificate['rref'][i][column])
            rational_checks += 1
    degrees = result['minimal_basis']['degrees']
    assert degrees == [0, 0, 0, 2, 2, 2, 3, 4, 4, 4]
    generators = result['minimal_basis']['generators']
    assert list(map(row_degree, generators)) == degrees
    leading = {}
    for row, degree in zip(generators, degrees):
        assert insert_binary(leading, sum(((value >> degree) & 1) << j for j, value in enumerate(row)))

    # The previously proved bounded saturation lemma applies at d-1=3.
    complete_space = bounded_polynomial_space(certificate, 3)
    shifted = {}
    for row, degree in zip(generators, degrees):
        for shift in range(max(0, 4 - degree)):
            assert insert_binary(shifted, encode([value << shift for value in row], 4))
    assert len(shifted) == complete_space['dimension'] == 19
    for row in complete_space['slabs']:
        assert not insert_binary(shifted, encode(row, 4))
    assert len(complete_space['equation_pivots']) == 21
    dimensions = {}
    for degree in (0, 1, 2, 3, 4, 8, 17):
        dimensions[degree] = bounded_polynomial_space(certificate, degree)['dimension']
        assert dimensions[degree] == sum(max(0, degree - d + 1) for d in degrees)

    endpoint_checks = 0
    for row, primitive in zip(generators, result['primitives']):
        assert primitive['normalized_slab'] == [value ^ (value << 1) for value in row]
        for value in primitive['normalized_slab']:
            assert value.bit_count() % 2 == 0
            endpoint_checks += 1
    # Literal dense independence, including all allowed monomial shifts.
    tensor_basis = {}
    for primitive, logical, shift in result['direction_descriptors']:
        slab = result['primitives'][primitive]['normalized_slab']
        vector = sum((value << shift) << (19 * (logical + 5 * j)) for j, value in enumerate(slab))
        assert insert_binary(tensor_basis, vector)
    assert len(tensor_basis) == 795
    endpoint_dimensions = {degree: 5 * sum(max(0, degree - d) for d in degrees)
                           for degree in (1, 3, 4, 5, 9, 18)}
    assert endpoint_dimensions == {1: 15, 3: 60, 4: 95, 5: 145, 9: 345, 18: 795}

    for value, factor in ((2, 4), (1, 0)):
        try:
            compress_exponents(value, factor)
        except ValueError:
            pass
        else:
            raise AssertionError('Malformed exponent root accepted.')
    return {'status': 'exact literal fourth-root family certified',
            'source_record_sha256': result['source_record_sha256'],
            'literal_root_coefficient_checks': coefficient_checks,
            'rational_reconstruction_entries': rational_checks,
            'minimal_basis_degrees': degrees,
            'saturation_checked_degree': 3, 'saturation_equation_rank': 21,
            'saturation_polynomial_dimension': 19, 'bounded_dimensions': dimensions,
            'endpoint_coefficient_checks': endpoint_checks,
            'endpoint_zero_tensor_dimensions': endpoint_dimensions,
            'literal_independent_directions_at_degree_18': len(tensor_basis),
            'negative_controls': 2, 'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
