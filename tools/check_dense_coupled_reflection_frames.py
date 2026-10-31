#!/usr/bin/env python3
"""Independent dense Laurent expansion controls for rank-one physical frames.

The reference represents each coefficient as an explicit exponent dictionary,
uses a separate F4 multiplication, and applies the four complete Laurent
matrices in Bell coordinates before restoring the original basis.  It never
requires an intermediate leg action to be polynomial.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from dense_coupled_reflection_frames import (
    rank_one_generators, rank_one_terms, transformed_from_terms, digest,
)
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import bell_normalize, coupled_admissible_poles


def multiply(a, b):
    product = 0
    while b:
        if b & 1:
            product ^= a
        a <<= 1
        b >>= 1
    if product & 4:
        product ^= 7
    return product


def coefficients(value):
    a, b = value
    return {exponent: ((a >> exponent) & 1) | (((b >> exponent) & 1) << 1)
            for exponent in range((a | b).bit_length())
            if ((a | b) >> exponent) & 1}


def add_term(target, exponent, value):
    result = target.get(exponent, 0) ^ value
    if result:
        target[exponent] = result
    elif exponent in target:
        del target[exponent]


def apply_matrix(values, matrix, leg, q, exponent=0, identity=False):
    """Apply I+w^exponent matrix, or just w^exponent matrix, on one leg."""
    stride = q ** leg
    output = [dict(value) if identity else {} for value in values]
    for index, polynomial in enumerate(values):
        source = index // stride % q
        for target in range(q):
            scalar = matrix[target][source]
            if scalar:
                destination = index + (target - source) * stride
                for power, coefficient in polynomial.items():
                    add_term(output[destination], power + exponent,
                             multiply(scalar, coefficient))
    return output


def matrix_product(first, second):
    rows, inner, cols = len(first), len(second), len(second[0])
    output = [[0] * cols for _ in range(rows)]
    for i in range(rows):
        for j in range(cols):
            for k in range(inner):
                output[i][j] ^= multiply(first[i][k], second[k][j])
    return output


def explicit_frame(seed, factors, left, right, pole):
    q = seed['q_x']
    values = [coefficients(value) for value in seed['values']]
    values = apply_matrix(values, factors['inverse_factor'], 0, q)
    values = apply_matrix(values, list(map(list, zip(*factors['first_factor']))), 3, q)
    generator = [[multiply(u, v) for v in right] for u in left]
    dual = list(map(list, zip(*generator)))
    for leg, matrix in enumerate((generator, dual, generator, dual)):
        values = apply_matrix(values, matrix, leg, q, -pole, identity=True)
    values = apply_matrix(values, factors['first_factor'], 0, q)
    values = apply_matrix(values, list(map(list, zip(*factors['inverse_factor']))), 3, q)
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    seed = tensor(4)
    _, factors = bell_normalize(seed)
    generators = list(rank_one_generators(4))
    assert len(generators) == 85 * 63 == 5355
    matrices = set()
    for left, right in generators:
        generator = [[multiply(u, v) for v in right] for u in left]
        flat = tuple(value for row in generator for value in row)
        assert any(flat) and flat not in matrices
        matrices.add(flat)
        assert matrix_product(generator, generator) == [[0] * 4 for _ in range(4)]
        # Every 2x2 minor vanishes.  Nonzero plus this condition is rank one.
        for i in range(4):
            for j in range(i):
                for a in range(4):
                    for b in range(a):
                        assert multiply(generator[i][a], generator[j][b]) == multiply(generator[i][b], generator[j][a])

    selected = [0, 16, 100, 252, 900, 1350, 2715, 4000, 5000, 5354]
    # Include an order-eight survivor, whose pole cancellation is not
    # explained by testing the separate paired actions.
    for index, (left, right) in enumerate(generators):
        terms = rank_one_terms(seed, left, right, factors)
        if 8 in coupled_admissible_poles(terms, 8):
            selected.append(index)
            break
    comparisons = endpoint_entries = rejected_poles = 0
    records = []
    for index in dict.fromkeys(selected):
        left, right = generators[index]
        terms = rank_one_terms(seed, left, right, factors)
        allowed = coupled_admissible_poles(terms, 12)
        assert 1 in allowed
        for pole in sorted(set([1, 2, max(allowed), max(allowed) + 1])):
            reference = explicit_frame(seed, factors, left, right, pole)
            is_polynomial = all(all(power >= 0 for power in value) for value in reference)
            assert is_polynomial == (pole in allowed)
            comparisons += len(reference)
            if not is_polynomial:
                rejected_poles += 1
                try:
                    transformed_from_terms(seed, terms, pole)
                except ArithmeticError:
                    pass
                else:
                    raise AssertionError('An independently witnessed pole was accepted.')
                continue
            actual = transformed_from_terms(seed, terms, pole)
            assert reference == [coefficients(value) for value in actual['values']]
            assert max(max(value, default=-1) for value in reference) == seed['degree']
            for offset, value in enumerate(reference):
                endpoint = 0
                for coefficient in value.values():
                    endpoint ^= coefficient
                x, y, u, v = (offset // 4 ** leg % 4 for leg in range(4))
                assert endpoint == int(x == v and y == u)
            endpoint_entries += len(reference)
        records.append({'generator_index': index, 'u': left, 'v': right,
                        'allowed_poles_through_12': allowed})

    catalogue_summary = None
    if args.catalogue:
        catalogue = json.loads(args.catalogue.read_text())
        assert catalogue['source_seed_sha256'] == digest(seed['values'])
        assert catalogue['full_q4_rank_one_catalogue']
        assert catalogue['generators_checked'] == len(generators)
        assert len(catalogue['records']) == len(generators)
        poles, hashes = Counter(), set()
        for index, ((left, right), record) in enumerate(zip(generators, catalogue['records'])):
            assert record['index'] == index and record['u'] == left and record['v'] == right
            terms = rank_one_terms(seed, left, right, factors)
            assert any(terms[order] for order in range(1, 5))
            allowed = coupled_admissible_poles(terms, catalogue['max_pole_checked'])
            assert allowed == record['pole_orders']
            assert len(record['candidates']) == len(allowed)
            for pole, candidate in zip(allowed, record['candidates']):
                actual = transformed_from_terms(seed, terms, pole)
                assert candidate == {'pole': pole, 'degree': 24, 'monomial_content': 0,
                                     'seed_sha256': digest(actual['values'])}
                poles[pole] += 1
                hashes.add(candidate['seed_sha256'])
        assert sum(poles.values()) == catalogue['polynomial_frames'] == 7176
        assert dict(poles) == {1: 5355, 2: 1755, 3: 21, 4: 21, 5: 6, 6: 6, 7: 6, 8: 6}
        assert len(hashes) == 7176
        catalogue_summary = {'sha256': sha256(args.catalogue.read_bytes()).hexdigest(),
                             'pole_counts': dict(sorted(poles.items())),
                             'distinct_polynomial_tensors': len(hashes),
                             'all_positive_poles_covered': catalogue['max_pole_checked'] >= seed['degree']}

    sources = [Path(__file__).name, 'dense_coupled_reflection_frames.py',
               'reflection_laurent_physical_frames.py', 'reflection_dual_rail_gaussian.py',
               'physical_split_crossings.py']
    result = {'status': 'passed', 'scope': 'exact generators, independent complete Laurent coefficient controls, optional exhaustive catalogue replay',
              'unique_nonzero_square_zero_rank_one_generators': len(matrices),
              'independent_laurent_tensor_entries': comparisons,
              'independent_endpoint_entries': endpoint_entries,
              'independent_rejected_poles': rejected_poles,
              'selected_generators': records, 'catalogue': catalogue_summary,
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in sources}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('selected_generators', 'source_sha256')}))


if __name__ == '__main__':
    main()
