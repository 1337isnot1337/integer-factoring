#!/usr/bin/env python3
"""Check the fan's two-state polynomial oracle and restricted-rank interface."""

import argparse
from collections import Counter
from functools import reduce
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import random

from check_modular_crossing import rank_two
from check_ward_fan_family import fan_types
from check_ward_pfaffian_queries import direct_coefficient
from complement_shell_tensor import multiply
from ward_fan_entry_queries import fan_entry_query
from ward_pfaffian_queries import BinaryField, squarefree_coefficient_stream


def columns(forms):
    types, _ = fan_types(forms)
    return [[(value >> a) & 1 for a in range(forms)] for value in types for _ in range(2)]


def oracle_controls():
    counts = Counter()
    rng = random.Random(792211)
    binary = BinaryField(3)
    for forms in (1, 2):
        weights = columns(forms)
        for values in product(product(range(2), repeat=2), repeat=len(weights)):
            direct = direct_coefficient(forms, [(a, m, b) for a, (m, b) in zip(weights, values)], binary)
            result = fan_entry_query(forms, iter(values), binary.mul)
            assert result == direct
            counts['exhaustive_binary_queries'] += 1
            counts['nonzero_binary_queries'] += bool(result)
    for modulus in (7, 11, 19):
        field = BinaryField(modulus)
        for forms in range(1, 9):
            weights = columns(forms)
            for _ in range(25):
                values = [(rng.randrange(field.size), rng.randrange(field.size)) for _ in weights]
                sites = [(a, m, b) for a, (m, b) in zip(weights, values)]
                result = fan_entry_query(forms, iter(values), field.mul)
                assert result == direct_coefficient(forms, sites, field)
                assert result == squarefree_coefficient_stream(forms, iter(sites), field).value
                counts['extension_field_three_way_queries'] += 1
                counts['nonzero_extension_queries'] += bool(result)

    class PolynomialRing:
        @staticmethod
        def mul(first, second):
            return multiply(first, second)

    for forms in range(1, 4):
        weights = columns(forms)
        degree = len(weights) + forms
        stride = degree + 1
        for _ in range(30):
            values = []
            for _ in weights:
                m = sum(rng.randrange(2) << exponent for exponent in (0, 1, stride))
                b = sum(rng.randrange(2) << exponent for exponent in (0, 1, stride, stride + 1))
                values.append((m, b))
            sites = [(a, m, b) for a, (m, b) in zip(weights, values)]
            result = fan_entry_query(forms, iter(values), PolynomialRing.mul)
            assert result == direct_coefficient(forms, sites, PolynomialRing())
            counts['full_bivariate_polynomial_queries'] += 1
    # A ring with nonzero nilpotents: no inversion or field hypothesis is used.
    class DualNumberRing:
        @staticmethod
        def mul(first, second):
            return multiply(first, second) & 3

    for forms in range(1, 5):
        weights = columns(forms)
        for _ in range(40):
            values = [(rng.randrange(4), rng.randrange(4)) for _ in weights]
            sites = [(a, m, b) for a, (m, b) in zip(weights, values)]
            assert fan_entry_query(forms, iter(values), DualNumberRing.mul) == direct_coefficient(forms, sites, DualNumberRing())
            counts['nonreduced_commutative_ring_queries'] += 1
    return dict(counts)


def wedge(linear_forms, size):
    result = {0}
    for form in linear_forms:
        following = set()
        for term in result:
            for i in range(size):
                if form >> i & 1 and not term >> i & 1:
                    value = term | (1 << i)
                    if value in following:
                        following.remove(value)
                    else:
                        following.add(value)
        result = following
    return result


def restricted_interface():
    forms = 3
    types, _ = fan_types(forms)
    labels, preimages, images = {}, {}, {}
    sectors = []
    for size in (3, 4):
        for chosen in combinations(range(len(types)), size):
            selected = [types[i] for i in chosen]
            if rank_two(selected) != forms:
                continue
            rows = [sum(((value >> a) & 1) << (2 * j + copy)
                        for j, value in enumerate(selected) for copy in range(2))
                    for a in range(forms)]
            sigma = ((1 << len(types)) - 1) ^ sum(1 << i for i in chosen)
            orientations = 1 << (2 * (len(types) - size) - 1)
            if size == forms:
                residuals = [None]
                actual = [wedge(rows, 2 * size)]
                preimage_sizes = [1]
                assert len(actual[0]) == 8
                pivots = [min(actual[0])]
            else:
                cycle = next(mask for mask in range(1, 1 << size)
                             if not reduce(int.__xor__, [selected[i] for i in range(size) if mask >> i & 1], 0))
                first = next(i for i in range(size) if cycle >> i & 1)
                residuals = [3 << (2 * first),
                             sum(1 << (2 * i + 1) for i in range(size) if cycle >> i & 1)]
                assert all((row & residual).bit_count() % 2 == 0 for row in rows for residual in residuals)
                actual = [wedge(rows + [residual], 2 * size) for residual in residuals]
                preimage_sizes = [residual.bit_count() for residual in residuals]
                assert actual[0] - actual[1] and actual[1] - actual[0]
                pivots = [min(actual[0] - actual[1]), min(actual[1] - actual[0])]
                assert len(actual[0]) == 16
                assert len(actual[1]) == (32 if cycle.bit_count() == 3 else 56)
            for a, image in enumerate(actual):
                for b, pivot in enumerate(pivots):
                    assert (pivot in image) == (a == b)
            labels[sigma] = orientations * len(actual)
            preimages[sigma] = orientations * sum(preimage_sizes)
            images[sigma] = orientations * sum(map(len, actual))
            sectors.append({'type_syndrome': sigma, 'selected_type_indices': chosen,
                            'orientation_count': orientations, 'labels': labels[sigma],
                            'preimage_support_sizes': preimage_sizes,
                            'input_support_sizes': list(map(len, actual)),
                            'coordinate_output_pivots': pivots})
    blocks = []
    for syndrome in range(1 << len(types)):
        size = sum(value * labels.get(sigma ^ syndrome, 0) for sigma, value in labels.items())
        left = sum(value * preimages.get(sigma ^ syndrome, 0) for sigma, value in labels.items())
        right = sum(value * images.get(sigma ^ syndrome, 0) for sigma, value in labels.items())
        if size:
            blocks.append({'type_pair_syndrome': syndrome, 'coefficient_matrix_dimension': size,
                           'uncached_physical_polynomial_queries': left * right})
    total = sum(labels.values())
    assert total == 84 and len(blocks) == 31
    assert sum(row['coefficient_matrix_dimension'] for row in blocks) == total ** 2
    assert max(row['coefficient_matrix_dimension'] for row in blocks) == 592
    return {'pairs': 10, 'forms': forms, 'length': 10, 'labels': total,
            'sectors': sectors, 'blocks': blocks,
            'input_preimage_support_sum': sum(preimages.values()),
            'input_image_support_sum': sum(images.values()),
            'dense_coefficient_entries_after_syndrome_split': sum(row['coefficient_matrix_dimension'] ** 2 for row in blocks),
            'dense_coefficient_entries_without_split': total ** 4,
            'uncached_physical_polynomial_queries_for_all_blocks': sum(row['uncached_physical_polynomial_queries'] for row in blocks),
            'scope': 'Exact sparse code and block interface; no restricted coefficient matrix has been assembled and no rank is inferred.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    names = (Path(__file__).name, 'ward_fan_entry_queries.py', 'ward_pfaffian_queries.py',
             'check_ward_pfaffian_queries.py', 'check_ward_fan_family.py', 'complement_shell_tensor.py')
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in names},
              'oracle_controls': oracle_controls(), 'restricted_rank_interface': restricted_interface(),
              'scope': 'Division-free formal physical-entry queries and exact sparse interface; no rank or factoring-score improvement asserted.'}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
