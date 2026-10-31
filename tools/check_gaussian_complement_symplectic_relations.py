#!/usr/bin/env python3
"""Symplectic contraction relations and the two complementary source sectors."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import random

from source_slab_module_saturation import (
    bounded_polynomial_space, encode, insert_binary, minimal_generators,
    rational_row_basis, row_degree,
)


def relations_and_parts(source):
    code = source['code_masks']
    pairs = [(1 << i) | (1 << (7 - i)) for i in range(4)]
    relations = []
    for positions in combinations(range(8), 2):
        beta = sum(1 << i for i in positions)
        relations.append([beta | pair for pair in pairs if not beta & pair])
    relations.append([first | second for first, second in combinations(pairs, 2)])
    parts = {'admissible': [], 'nonadmissible': []}
    for item in source['raw_slabs']:
        mask = item['complement_mask']
        admissible = all((mask & pair).bit_count() == 1 for pair in pairs)
        parts['admissible' if admissible else 'nonadmissible'].append(item)
    assert all(mask not in code for relation in relations for mask in relation)
    return pairs, relations, parts


def weak_popov_degrees(rows, width=125):
    encoded = []
    for row in rows:
        value = 0
        for coordinate, polynomial in enumerate(row):
            while polynomial:
                bit = polynomial & -polynomial
                polynomial ^= bit
                value ^= 1 << (width * (bit.bit_length() - 1) + coordinate)
        encoded.append(value)
    pivots, operations = {}, 0
    for value in encoded:
        while value:
            degree, place = divmod(value.bit_length() - 1, width)
            if place not in pivots:
                pivots[place] = value
                break
            other = pivots[place]
            other_degree = (other.bit_length() - 1) // width
            if degree < other_degree:
                value, pivots[place] = other, value
                degree, other_degree = other_degree, degree
                other = pivots[place]
            value ^= other << (width * (degree - other_degree))
            operations += 1
        assert value, 'The supplied sector generators were not independent.'
    return sorted((value.bit_length() - 1) // width for value in pivots.values()), operations


def transpose(rows, size):
    return [sum(((row >> column) & 1) << i for i, row in enumerate(rows))
            for column in range(size)]


def reflected(value):
    return sum(((value >> i) & 1) << (8 * (i // 8) + 7 - i % 8) for i in range(16))


def binary_product(left, right):
    result = []
    for row in left:
        value = 0
        for j, other in enumerate(right):
            if row >> j & 1:
                value ^= other
        result.append(value)
    return result


def minor(rows, row_mask, column_mask):
    columns = [j for j in range(16) if column_mask >> j & 1]
    matrix = [sum(((row >> j) & 1) << k for k, j in enumerate(columns))
              for i, row in enumerate(rows) if row_mask >> i & 1]
    assert len(matrix) == len(columns) == 8
    pivots = {}
    for row in matrix:
        insert_binary(pivots, row)
    return int(len(pivots) == 8)


def generic_symplectic_controls(code, labels, relations):
    rng = random.Random(90751)
    controls = 0
    nonorthogonal = 0
    for sample in range(2):
        matrix = [1 << i for i in range(16)]
        for step in range(20):
            vector = 1 if step == 0 else rng.randrange(1, 1 << 16)
            covector = reflected(vector)
            transvection = [(1 << i) ^ (covector if vector >> i & 1 else 0) for i in range(16)]
            matrix = binary_product(transvection, matrix)
        columns = transpose(matrix, 16)
        assert all(((columns[i] & reflected(columns[j])).bit_count() % 2)
                   == int(j == 8 * (i // 8) + 7 - i % 8)
                   for i in range(16) for j in range(16))
        nonorthogonal += any(sum(((column >> (8 * block + i)) & 1)
                                * ((column >> (8 * block + 7 - i)) & 1)
                                for block in range(2) for i in range(4)) % 2
                             for column in columns)
        for _ in range(16):
            u, v, y = (rng.choice(code) for _ in range(3))
            values = {alpha: minor(matrix, v | (y << 8), u | (alpha << 8)) for alpha in labels}
            for relation in relations:
                assert sum(values[alpha] for alpha in relation) % 2 == 0
                controls += 1
    assert nonorthogonal == 2
    return controls


def verify():
    path = Path('evidence/aligned-gaussian-complement-directions-q5.json')
    source = json.loads(path.read_text())
    module = json.loads(Path('evidence/aligned-gaussian-module-directions-q5.json').read_text())
    pairs, relations, parts = relations_and_parts(source)
    masks = [item['complement_mask'] for item in source['raw_slabs']]
    raw = {item['complement_mask']: item['raw_slab_in_z'] for item in source['raw_slabs']}
    relation_pivots = {}
    polynomial_checks = 0
    for relation in relations:
        for coordinate in range(125):
            value = 0
            for mask in relation:
                value ^= raw[mask][coordinate]
            assert value == 0
            polynomial_checks += 1
        insert_binary(relation_pivots, sum(1 << masks.index(mask) for mask in relation))
    first_rank = {}
    for relation in relations[:-1]:
        insert_binary(first_rank, sum(1 << masks.index(mask) for mask in relation))
    assert len(first_rank) == 27 and len(relation_pivots) == 28

    sector_records = {}
    generators = []
    constant_total = {}
    width = max(value.bit_length() for row in raw.values() for value in row)
    for name, items in parts.items():
        rows = [item['raw_slab_in_z'] for item in items]
        constant_rank = {}
        for row in rows:
            insert_binary(constant_rank, encode(row, width))
            insert_binary(constant_total, encode(row, width))
        rational = rational_row_basis(rows)
        minimal = minimal_generators(rational, 24)
        degrees = minimal['degrees']
        assert len(degrees) == rational['rank'] == len(constant_rank)
        bound = max(degrees) - 1
        complete = bounded_polynomial_space(rational, bound)
        span = {}
        for row, degree in zip(minimal['generators'], degrees):
            for shift in range(max(0, bound - degree + 1)):
                assert insert_binary(span, encode([value << shift for value in row], bound + 1))
        assert len(span) == complete['dimension']
        for row in complete['slabs']:
            assert not insert_binary(span, encode(row, bound + 1))
        sector_records[name] = {'raw_complements': len(items), 'constant_and_rational_rank': rational['rank'],
                               'saturated_minimal_degrees': degrees,
                               'saturation_checked_degree': bound,
                               'saturation_bounded_dimension': complete['dimension'],
                               'endpoint_tensor_dimension_at_degree18': 5 * sum(max(0, 18 - d) for d in degrees)}
        generators.extend(minimal['generators'])
    assert len(constant_total) == 37
    assert sector_records['admissible']['constant_and_rational_rank'] == 11
    assert sector_records['nonadmissible']['constant_and_rational_rank'] == 26
    combined_degrees, operations = weak_popov_degrees(generators)
    assert len(combined_degrees) == 37 and sum(combined_degrees) == 335
    assert sum(module['minimal_basis']['degrees']) == 242
    combined_dimension = 5 * sum(max(0, 18 - d) for d in combined_degrees)
    assert combined_dimension == 1655
    combined_endpoint = {}
    for row in generators:
        insert_binary(combined_endpoint, sum((value.bit_count() % 2) << j for j, value in enumerate(row)))
    assert len(combined_endpoint) == 28
    combined_endpoint_zero_dimension = 5 * (sum(max(0, 19 - d) for d in combined_degrees)
                                             - len(combined_endpoint))
    assert combined_endpoint_zero_dimension == 1700
    full_endpoint = {}
    for row in module['minimal_basis']['generators']:
        insert_binary(full_endpoint, sum((value.bit_count() % 2) << j for j, value in enumerate(row)))
    for x in range(5):
        crossing = sum(int(v == x and y == u) << (y + 5 * u + 25 * v)
                       for v in range(5) for u in range(5) for y in range(5))
        assert insert_binary(full_endpoint, crossing)
    assert len(full_endpoint) == 42
    labels = sorted(masks + source['code_masks'])
    generic_checks = generic_symplectic_controls(source['code_masks'], labels, relations)
    return {'status': 'exact symplectic contraction relations and complementary sectors certified',
            'source_record_sha256': sha256(path.read_bytes()).hexdigest(),
            'reflection_pairs': pairs, 'relation_count': 29,
            'first_contraction_rank': 27, 'total_relation_rank': 28,
            'complete_polynomial_relation_entries': polynomial_checks,
            'raw_complement_constant_rank': len(constant_total),
            'sector_modules': sector_records,
            'combined_unsaturated_weak_popov_degrees': combined_degrees,
            'combined_unimodular_row_operations': operations,
            'separately_degree18_bounded_endpoint_dimension': 1250,
            'combined_restored_factor_module_dimension_at_degree18': combined_dimension,
            'combined_module_endpoint_rank': len(combined_endpoint),
            'complete_combined_module_endpoint_zero_dimension_at_degree18': combined_endpoint_zero_dimension,
            'fully_saturated_endpoint_dimension_at_degree18': 2120,
            'full_primitive_endpoint_rank': len(full_endpoint),
            'polynomial_quotient_dimension_of_full_saturation': 93,
            'generic_nonorthogonal_symplectic_contraction_controls': generic_checks,
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
