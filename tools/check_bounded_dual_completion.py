#!/usr/bin/env python3
"""Independent source and affine-coordinate controls for bounded duals."""

from itertools import product
import json
from pathlib import Path
import random

from bounded_dual_completion import full_weight_completion, source_monomial_entry
from check_modular_crossing import regular_partitions, determinant_operator_two
from check_modular_parameter_descent import reduced_dual
from homogeneous_dual_completion import specialize


def kernel_basis(slabs):
    width = max((v.bit_length() for slab in slabs for v in slab), default=1)
    pivots, independent, relations = {}, [], []
    for index, slab in enumerate(slabs):
        vector = sum(value << (width * j) for j, value in enumerate(slab))
        relation = 1 << index
        while vector:
            pivot = vector.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = vector, relation
                independent.append(index)
                break
            previous, combination = pivots[pivot]
            vector ^= previous
            relation ^= combination
        if not vector:
            relations.append([i for i in range(len(slabs)) if relation >> i & 1])
    for relation in relations:
        for j in range(len(slabs[0])):
            value = 0
            for i in relation:
                value ^= slabs[i][j]
            assert value == 0
    return independent, relations


def verify():
    source_checks = 0
    h = 8
    labels = regular_partitions(2, h)
    for w in range(h + 1):
        for a in regular_partitions(2, w):
            for y, v, u in product(labels, repeat=3):
                start, bits = source_monomial_entry(16, a, y, v, u)
                assert source_monomial_entry(32, a, y, v, u) == (start, bits)
                for shift in range(-2 * h, 2 * h + 1):
                    expected = ((bits >> (shift - start)) & 1) if shift >= start else 0
                    actual = 0
                    for i in range(h + w + 1):
                        if i & ~(h + w):
                            continue
                        for middle in regular_partitions(2, 2 * h - shift - i):
                            actual ^= (determinant_operator_two(a, u, middle)
                                       & determinant_operator_two(v, middle, y))
                    assert actual == expected
                    source_checks += 1

    records = []
    endpoint_checks = vacuum_checks = 0
    rng = random.Random(90616)
    for weight in (8, 16):
        seed = full_weight_completion(weight)
        retained = json.loads(Path(f'evidence/bounded-dual-completion-h{weight}.json').read_text())
        retained.pop('kernel_tensor_rank')
        assert json.loads(json.dumps(seed)) == retained
        q = seed['dimensions'][0]
        labels = list(map(tuple, seed['labels']))
        independent, relations = kernel_basis(seed['kernel_slabs'])
        assert len(independent) == (2 if weight == 8 else 10)
        for terms in seed['kernel_parameters']:
            for b in labels:
                assert sum(determinant_operator_two(tuple(a), (), b) for a in terms) % 2 == 0
                vacuum_checks += 1
        lower_parameters = [m for w in range(weight) for m in regular_partitions(2, w)]
        changes = []
        for x, b in enumerate(labels):
            for a, coefficient in reduced_dual(2, b):
                if sum(a) < weight and coefficient:
                    changes.append((x, lower_parameters.index(a), 1))
        assert [pair[0] for pair in specialize(seed, changes)] == seed['integral_baseline']
        for _ in range(8):
            changes = [(x, i, rng.randrange(4)) for x in range(q) for i in independent]
            values = specialize(seed, changes)
            for index, (low, high) in enumerate(values):
                x = index % q
                y = index // q % q
                u = index // q ** 2 % q
                v = index // q ** 3
                assert ((low.bit_count() & 1) | ((high.bit_count() & 1) << 1)) == int(x == v and y == u)
                endpoint_checks += 1
        records.append({'weight': weight, 'alphabet': q,
                        'parameter_dimension': seed['full_parameter_dimension'],
                        'formal_dual_direction_dimension': seed['dual_parameter_dimension'],
                        'kernel_slab_rank': len(independent),
                        'actual_tensor_direction_dimension': q * len(independent),
                        'independent_kernel_indices': independent,
                        'coefficientwise_relations': relations,
                        'common_laurent_shift': seed['common_laurent_shift'],
                        'verified_exponent_step': seed['verified_exponent_step']})
    return {'status': 'exact bounded-source compiler and dual-completion controls passed',
            'literal_source_coefficients': source_checks,
            'kernel_vacuum_entries': vacuum_checks,
            'F4_endpoint_entries': endpoint_checks,
            'integral_duals_recovered': True, 'cases': records,
            'factoring_bound_improved': False}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
