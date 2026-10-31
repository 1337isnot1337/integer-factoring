"""Constant coefficient closure of polynomial source slabs.

Repeated characteristic-two variable squaring and rational-span
saturation recover precisely this constant coefficient space.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from source_slab_module_saturation import insert_binary, row_degree


def coefficient_basis(rows):
    pivots, selected = {}, []
    for index, row in enumerate(rows):
        for power in range(row_degree(row) + 1):
            vector = sum(((value >> power) & 1) << column for column, value in enumerate(row))
            if insert_binary(pivots, vector):
                selected.append({'generator_index': index, 'coefficient_power': power,
                                 'constant_slab_bits': vector})
    return selected


def crossing_pair_permutations():
    generators = ((3, 1, 2, 0), (0, 2, 1, 3), (1, 0, 3, 2))
    group = {(0, 1, 2, 3)}
    while True:
        enlarged = group | {tuple(first[second[i]] for i in range(4))
                            for first in group for second in generators}
        if enlarged == group:
            return sorted(group)
        group = enlarged


def permute_tensor_bits(value, dimension, permutation):
    result = 0
    while value:
        bit = value & -value
        value ^= bit
        index = bit.bit_length() - 1
        coordinates = [(index // dimension ** leg) % dimension for leg in range(4)]
        target = sum(coordinates[permutation[leg]] * dimension ** leg for leg in range(4))
        result ^= 1 << target
    return result


def constant_tensor_basis(slabs, dimension):
    result = []
    for item in slabs:
        for first in range(dimension):
            value = item['constant_slab_bits']
            result.append(sum(((value >> j) & 1) << (first + dimension * j)
                              for j in range(dimension ** 3)))
    return result


def closure_record(source, degree=18):
    dimension = source['dimensions'][0]
    if source['dimensions'] != [dimension] * 4:
        raise ValueError('Use four equal physical alphabets.')
    rows = source['minimal_basis']['generators']
    selected = coefficient_basis(rows)
    primitives = []
    for item in selected:
        values = [3 * ((item['constant_slab_bits'] >> j) & 1) for j in range(dimension ** 3)]
        primitives.append(item | {'normalized_slab': values, 'normalized_degree': 1})
    original = constant_tensor_basis(selected, dimension)
    pivots, orbit_basis, prefix = {}, [], []
    group = crossing_pair_permutations()
    for permutation in group:
        for index, value in enumerate(original):
            transformed = permute_tensor_bits(value, dimension, permutation)
            if insert_binary(pivots, transformed):
                orbit_basis.append({'original_constant_direction': index,
                                    'permutation': permutation,
                                    'constant_tensor_bits': transformed})
        prefix.append(len(pivots))
    descriptors = [(primitive, first, shift)
                   for primitive in range(len(primitives))
                   for first in range(dimension) for shift in range(degree)]
    return {'constructor': 'source_slab_frobenius_coefficient_closure',
            'dimensions': source['dimensions'], 'degree': degree,
            'constant_slab_dimension': len(selected),
            'primitives': primitives, 'direction_descriptors': descriptors,
            'direction_count_before_rank_reduction': len(descriptors),
            'unshifted_tensor_dimension': len(original),
            'crossing_pair_permutations': group,
            'permutation_prefix_ranks': prefix,
            'permutation_constant_tensor_basis': orbit_basis,
            'permutation_constant_tensor_dimension': len(orbit_basis),
            'permutation_endpoint_dimension_at_degree_bound': degree * len(orbit_basis),
            'scope': 'Exact coefficient closure under finite Frobenius reparametrization and rational saturation; each displayed endpoint direction includes a factor1+w.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--degree', type=int, default=18)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = closure_record(json.loads(args.source.read_text()), args.degree)
    result['source_record_sha256'] = sha256(args.source.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in ('constant_slab_dimension',
        'unshifted_tensor_dimension', 'permutation_constant_tensor_dimension',
        'direction_count_before_rank_reduction', 'permutation_endpoint_dimension_at_degree_bound')}))


if __name__ == '__main__':
    main()
