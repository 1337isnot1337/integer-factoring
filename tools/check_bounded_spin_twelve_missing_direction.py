#!/usr/bin/env python3
"""Test the twelve-site analogue of the ten-site completing coefficient direction."""
from __future__ import annotations

from hashlib import sha256
from math import comb
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

from bounded_spin_crossing import paths
from bounded_spin_six_vertex import monomial_entry
from build_bounded_spin_twelve_coefficient_table import solve_mod
from check_bounded_spin_transposition_subalgebra import (
    frame, transposition_matrices,
)
from native_gf5_flint import rank

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / 'evidence/bounded-spin-twelve-coefficient-table.json'
DIMENSION = ROOT / 'evidence/bounded-spin-twelve-coefficient-dimension.json'
IDEMPOTENTS = ROOT / 'evidence/bounded-spin-twelve-coefficient-idempotents.json'
TEN_MODEL = ROOT / 'evidence/bounded-spin-ten-algebra-model.json'
OUTPUT = ROOT / 'evidence/bounded-spin-twelve-missing-direction.json'
TEN_BASIS = [0, 2, 3, 4, 5, 6, 7, 10]
TEN_WITNESS = [(27, 26, 4, 15), (14, 2, 8, 18), (29, 23, 19, 15),
               (11, 15, 5, 19), (11, 9, 14, 14), (2, 4, 2, 4),
               (10, 10, 9, 9), (5, 16, 16, 5)]
LENGTH = 12
PRIME = 5


def original_form(four_paths):
    length = len(four_paths[0])
    result = [0] * (length + 1)
    for (_, b, c), value in monomial_entry(four_paths).items():
        for j in range(2 * b + 1):
            result[2 * c + j] = (result[2 * c + j]
                                 + value * comb(2 * b, j) * (-1) ** j) % PRIME
    return result


def generated_algebra(unit, regular, generators):
    def left(vector):
        return (np.tensordot(vector, regular, axes=(0, 0)) % PRIME).astype(np.uint8)

    def product(left_vector, right_vector):
        return (left(left_vector) @ right_vector) % PRIME

    space = [unit]
    frontier = [unit]
    while frontier:
        new = []
        for vector in frontier:
            for generator in generators:
                candidate = product(vector, generator)
                if not any(np.array_equal(candidate, known) for known in space):
                    if rank(np.array(space + [candidate], dtype=np.uint8)) > len(space):
                        space.append(candidate)
                        new.append(candidate)
        frontier = new
    return np.array(space, dtype=np.uint8)


def coordinate_product(left_vector, right_vector, regular):
    return (np.tensordot(left_vector, regular, axes=(0, 0)) @ right_vector) % PRIME


def trace_weights(table):
    return np.array([sum(int(table[a][b][b]) for b in range(len(table))) % PRIME
                     for a in range(len(table))], dtype=np.uint8)


def trace_left(element, weights):
    return int(element @ weights % PRIME)


def block_character(element, idempotent, weights, regular):
    product = coordinate_product(element, idempotent, regular)
    numerator = trace_left(product, weights)
    denominator = trace_left(idempotent, weights)
    assert denominator % PRIME
    return numerator * pow(denominator, -1, PRIME) % PRIME


def ten_site_control():
    all_paths = paths(10)
    rows = [[original_form(tuple(all_paths[i] for i in ids))[k] for k in TEN_BASIS]
            for ids in TEN_WITNESS]
    entry_matrix = np.array(rows, dtype=np.uint8)
    assert rank(entry_matrix) == 8
    _, words, vectors, inverse_norms = frame(10)
    transpositions = transposition_matrices(10, words, vectors, inverse_norms)
    replica_rhs = [int(x == u and y == v) for x, y, u, v in TEN_WITNESS]
    class_rhs = [sum(int(matrix[x, v]) * int(matrix[y, u])
                     for matrix in transpositions) % PRIME
                 for x, y, u, v in TEN_WITNESS]
    replica = solve_mod(entry_matrix.tolist(), replica_rhs)
    class_sum = solve_mod(entry_matrix.tolist(), class_rhs)
    assert replica == [1, 0, 0, 0, 0, 0, 0, 0]
    assert class_sum == [0, 2, 0, 1, 0, 2, 2, 0]

    model = json.loads(TEN_MODEL.read_text())
    model_unit = np.array(model['unit'], dtype=np.uint8)
    model_regular = np.array(model['regular'], dtype=np.uint8)
    model_basis = model['basis_indices']
    assert model_basis == TEN_BASIS
    b2 = np.zeros(8, dtype=np.uint8); b2[2] = 1
    b4 = np.zeros(8, dtype=np.uint8); b4[4] = 1
    rt_dimension = len(generated_algebra(model_unit, model_regular,
                                         [np.array(class_sum, dtype=np.uint8),
                                          np.array(replica, dtype=np.uint8)]))
    b4_dimension = len(generated_algebra(model_unit, model_regular,
                                         [np.array(class_sum, dtype=np.uint8),
                                          np.array(replica, dtype=np.uint8), b4]))
    b2_b4_sum_dimension = len(generated_algebra(model_unit, model_regular,
                                                [np.array(class_sum, dtype=np.uint8),
                                                 np.array(replica, dtype=np.uint8),
                                                 (b2 + b4) % PRIME]))
    b2_b4_pair_dimension = len(generated_algebra(model_unit, model_regular,
                                                 [np.array(class_sum, dtype=np.uint8),
                                                  np.array(replica, dtype=np.uint8), b2, b4]))
    assert (rt_dimension, b4_dimension, b2_b4_sum_dimension,
            b2_b4_pair_dimension) == (6, 8, 7, 8)
    return {'witness_count': len(TEN_WITNESS), 'replica_coordinates': replica,
            'class_sum_coordinates': class_sum,
            'matches_published_ten_site': True,
            'generated_dimension': rt_dimension,
            'b4_generated_dimension': b4_dimension,
            'b2_b4_sum_generated_dimension': b2_b4_sum_dimension,
            'b2_b4_pair_generated_dimension': b2_b4_pair_dimension}


def main() -> None:
    started = time.monotonic()
    saved = json.loads(TABLE.read_text())
    dimension_saved = json.loads(DIMENSION.read_text())
    assert saved['passed'] and saved['length'] == LENGTH
    assert saved['basis_indices'] == [0, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    basis_indices = saved['basis_indices']
    entry_matrix = saved['entry_matrix']
    unit = saved['unit']
    table = np.array(saved['table'], dtype=np.uint8)
    regular = np.array(saved['regular'], dtype=np.uint8)
    ten_control = ten_site_control()

    all_paths = paths(LENGTH)
    q = len(all_paths)
    index = {tuple(path): i for i, path in enumerate(all_paths)}
    witness = [tuple(index[tuple(path)] for path in quadruple)
               for quadruple in saved['witness_quadruples']]
    assert len(witness) == 10

    _, words, vectors, inverse_norms = frame(LENGTH)
    transpositions = transposition_matrices(LENGTH, words, vectors, inverse_norms)
    unit_rhs = []
    replica_rhs = []
    class_rhs = []
    for x, y, u, v in witness:
        unit_rhs.append(int(x == u and y == v))
        replica_rhs.append(int(x == u and y == v))
        class_rhs.append(sum(int(matrix[x, v]) * int(matrix[y, u])
                             for matrix in transpositions) % PRIME)
    assert rank(np.array(entry_matrix, dtype=np.uint8)) == 10
    assert solve_mod(entry_matrix, unit_rhs) == unit
    replica_coordinates = np.array(solve_mod(entry_matrix, replica_rhs), dtype=np.uint8)
    class_coordinates = np.array(solve_mod(entry_matrix, class_rhs), dtype=np.uint8)

    idempotent_saved = json.loads(IDEMPOTENTS.read_text())
    weights = trace_weights(table)
    separator = np.zeros(10, dtype=np.uint8); separator[2] = 1; separator[4] = 1
    primitive_values = []
    for item in idempotent_saved['primitive_central_idempotents']:
        idempotent = np.array(item['idempotent'], dtype=np.uint8)
        assert trace_left(idempotent, weights) == item['regular_rank']
        values = {
            'index': item['index'],
            'kind': item['kind'],
            'regular_rank': item['regular_rank'],
            'class_sum': int(block_character(class_coordinates, idempotent, weights, regular)),
            'separator': int(block_character(separator, idempotent, weights, regular)),
            'replica': int(block_character(replica_coordinates, idempotent, weights, regular)),
        }
        primitive_values.append(values)
    assert [item['class_sum'] for item in primitive_values] == [0, 2, 3, 3, 0, 0, 1, 3]
    assert [item['separator'] for item in primitive_values] == [2, 0, 2, 3, 0, 3, 4, 1]
    assert all(item['replica'] == 1 for item in primitive_values)
    joint_pairs = [(item['class_sum'], item['separator']) for item in primitive_values]
    assert len(set(joint_pairs)) == 8
    character_separation = {
        'primitive_count': len(primitive_values),
        'joint_pairs': joint_pairs,
        'joint_pairs_distinct': True,
        'class_sum_distinct_values': len({item['class_sum'] for item in primitive_values}),
        'separator_distinct_values': len({item['separator'] for item in primitive_values}),
        'replica_is_unit': True,
    }

    subalgebra = generated_algebra(unit, regular, [class_coordinates, replica_coordinates])
    exact_dimension = len(subalgebra)
    single = []
    for coordinate_index in range(10):
        direction = np.zeros(10, dtype=np.uint8)
        direction[coordinate_index] = 1
        dimension = len(generated_algebra(unit, regular, [class_coordinates, replica_coordinates, direction]))
        single.append({'coordinate_index': coordinate_index,
                       'polynomial_coefficient': basis_indices[coordinate_index],
                       'generated_dimension': dimension})
    completing = [item for item in single if item['generated_dimension'] == 10]
    analogue = next(item for item in single if item['coordinate_index'] == 4)
    pairs = []
    sums = []
    for i in range(10):
        for j in range(i + 1, 10):
            first = np.zeros(10, dtype=np.uint8); first[i] = 1
            second = np.zeros(10, dtype=np.uint8); second[j] = 1
            pair_dimension = len(generated_algebra(unit, regular,
                                                   [class_coordinates, replica_coordinates,
                                                    first, second]))
            if pair_dimension == 10:
                pairs.append([i, j])
            sum_dimension = len(generated_algebra(unit, regular,
                                                  [class_coordinates, replica_coordinates,
                                                   (first + second) % PRIME]))
            if sum_dimension == 10:
                sums.append([i, j])
    separator_dimension = len(generated_algebra(unit, regular,
                                                [class_coordinates, replica_coordinates,
                                                 separator]))
    assert separator_dimension == 10
    assert [2, 4] in sums and [1, 5] in sums
    receipt = {
        'passed': True,
        'length': LENGTH,
        'prime': PRIME,
        'basis_indices': basis_indices,
        'ten_site_control': ten_control,
        'entry_witness_count': len(witness),
        'unit_coordinates_control': True,
        'primitive_character_values': primitive_values,
        'joint_character_separation': character_separation,
        'replica_coordinates': replica_coordinates.tolist(),
        'class_sum_coordinates': class_coordinates.tolist(),
        'generated_dimension': exact_dimension,
        'single_direction_scan': single,
        'completing_coordinate_indices': [item['coordinate_index'] for item in completing],
        'completing_polynomial_coefficients': [item['polynomial_coefficient'] for item in completing],
        'ten_site_analogue': analogue,
        'ten_site_analogue_completes': analogue['generated_dimension'] == 10,
        'completing_pairs': pairs,
        'completing_sums': sums,
        'separator_coordinates': [0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        'separator_polynomial_coefficients': [basis_indices[2], basis_indices[4]],
        'separator_generated_dimension': separator_dimension,
        'completing_sum_polynomial_coefficients': [
            [basis_indices[i], basis_indices[j]] for i, j in sums],
        'completing_pair_polynomial_coefficients': [
            [basis_indices[i], basis_indices[j]] for i, j in pairs],
        'table_sha256': sha256(TABLE.read_bytes()).hexdigest(),
        'dimension_receipt_sha256': sha256(DIMENSION.read_bytes()).hexdigest(),
        'ten_site_model_sha256': sha256(TEN_MODEL.read_bytes()).hexdigest(),
        'idempotent_receipt_sha256': sha256(IDEMPOTENTS.read_bytes()).hexdigest(),
        'source_sha256': {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                          for name in ('bounded_spin_crossing.py',
                                       'bounded_spin_six_vertex.py',
                                       'check_bounded_spin_transposition_subalgebra.py',
                                       'build_bounded_spin_twelve_coefficient_table.py')},
        'scope': ('A twelve-site structural diagnostic: exact R,T coordinates in the '
                  'published ten-dimensional coefficient table, the exact generated '
                  'dimension, a single basis-direction scan, and a linear-combination '
                  'scan of the ten-site separator direction, and the joint T/separator '
                  'character values on the eight primitive central idempotents. No '
                  'cut-rank bound or score change is claimed.'),
        'elapsed_seconds': time.monotonic() - started,
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'dim<R,T>': exact_dimension,
                      'b4_dimension': analogue['generated_dimension'],
                      'separator_dimension': separator_dimension,
                      'completing_indices': receipt['completing_coordinate_indices'],
                      'completing_sums': receipt['completing_sums'],
                      'completing_pairs': receipt['completing_pairs'],
                      'joint_character_separation': character_separation['joint_pairs_distinct'],
                      'replica': replica_coordinates.tolist(),
                      'class_sum': class_coordinates.tolist(),
                      'seconds': receipt['elapsed_seconds']}))


if __name__ == '__main__':
    main()
