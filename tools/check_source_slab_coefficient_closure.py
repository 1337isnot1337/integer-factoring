#!/usr/bin/env python3
"""Exact coefficient-closure and crossing-preserving permutation controls."""

from hashlib import sha256
import json
from pathlib import Path

from source_slab_coefficient_closure import (
    closure_record, coefficient_basis, constant_tensor_basis, permute_tensor_bits,
)
from source_slab_module_saturation import insert_binary


def verify():
    path = Path('evidence/source-slab-module-directions-h16.json')
    source = json.loads(path.read_text())
    result = json.loads(Path('evidence/source-slab-coefficient-directions-h16.json').read_text())
    expected = closure_record(source, 18)
    expected['source_record_sha256'] = sha256(path.read_bytes()).hexdigest()
    assert json.loads(json.dumps(expected)) == result
    selected = result['primitives']
    assert len(selected) == 22
    pivots = {}
    for item in selected:
        index, power = item['generator_index'], item['coefficient_power']
        row = source['minimal_basis']['generators'][index]
        literal = sum(((value >> power) & 1) << j for j, value in enumerate(row))
        assert literal == item['constant_slab_bits']
        assert insert_binary(pivots, literal)
        assert item['normalized_slab'] == [3 * ((literal >> j) & 1) for j in range(125)]
    # The raw bounded source and its saturated basis have identical closures.
    bounded = json.loads(Path('evidence/bounded-dual-completion-h16.json').read_text())
    raw = coefficient_basis(bounded['kernel_slabs'])
    assert len(raw) == 22
    for item in raw:
        assert not insert_binary(pivots, item['constant_slab_bits'])

    unshifted = constant_tensor_basis(selected, 5)
    tensor_pivots = {}
    for row in unshifted:
        assert insert_binary(tensor_pivots, row)
    assert len(tensor_pivots) == 110
    group = list(map(tuple, result['crossing_pair_permutations']))
    assert len(group) == 8 and len(set(group)) == 8
    crossing = sum(int(index % 5 == index // 125 and index // 5 % 5 == index // 25 % 5) << index
                   for index in range(625))
    orbit_pivots = {}
    for permutation in group:
        assert permute_tensor_bits(crossing, 5, permutation) == crossing
        for row in unshifted:
            insert_binary(orbit_pivots, permute_tensor_bits(row, 5, permutation))
    assert len(orbit_pivots) == 367
    retained_pivots = {}
    for item in result['permutation_constant_tensor_basis']:
        value = permute_tensor_bits(unshifted[item['original_constant_direction']], 5, item['permutation'])
        assert value == item['constant_tensor_bits']
        assert insert_binary(retained_pivots, value)
    assert len(retained_pivots) == 367

    # Degree independence follows from independent leading coefficient tensors.
    # Check the full first-leg export literally as an additional control.
    directions = {}
    for primitive, first, shift in result['direction_descriptors']:
        slab = result['primitives'][primitive]['normalized_slab']
        packed = sum((value << shift) << (19 * (first + 5 * j)) for j, value in enumerate(slab))
        assert insert_binary(directions, packed)
    assert len(directions) == 1980
    assert result['permutation_endpoint_dimension_at_degree_bound'] == 18 * 367 == 6606

    pure_controls = []
    for index, supports in ((8, ((0,), (2, 3), (1, 2, 3, 4))),
                            (19, ((2, 4), (2, 4), (2, 4)))):
        value = sum(1 << (y + 5 * u + 25 * v)
                    for y in supports[0] for u in supports[1] for v in supports[2])
        assert value == selected[index]['constant_slab_bits']
        pure_controls.append({'basis_index': index, 'three_leg_supports': supports})
    return {'status': 'exact coefficient closure and permutation orbit certified',
            'source_record_sha256': result['source_record_sha256'],
            'source_and_saturated_constant_slab_dimension': 22,
            'unshifted_first_leg_tensor_dimension': 110,
            'crossing_preserving_permutations': 8,
            'permutation_constant_tensor_dimension': 367,
            'literal_degree18_direction_rank': 1980,
            'permutation_endpoint_dimension_at_degree18': 6606,
            'pure_tensor_controls': pure_controls,
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
