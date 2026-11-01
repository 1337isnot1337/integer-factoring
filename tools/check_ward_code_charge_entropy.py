#!/usr/bin/env python3
"""Exact dimension and positive-parameter bounds for binary Ward codes."""

import argparse
from fractions import Fraction
import hashlib
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path

from check_modular_crossing import rank_two


def charge(column, forms):
    result, bit = 0, 0
    for i in range(forms):
        for j in range(i, forms):
            result |= ((column >> i & 1) * (column >> j & 1)) << bit
            bit += 1
    return result


def exact_dimensions(columns, forms):
    pairs = len(columns)
    charges = [charge(column, forms) for column in columns]
    charge_rank = rank_two(charges)
    dimension = [0] * (2 * pairs + 1)
    kernel_count = admissible = 0
    for mask in range(1 << pairs):
        selected = [i for i in range(pairs) if mask >> i & 1]
        value = 0
        for i in selected:
            value ^= charges[i]
        if value:
            continue
        kernel_count += 1
        if rank_two([columns[i] for i in selected]) != forms:
            continue
        admissible += 1
        even = len(selected)
        assert even >= 2 * forms
        singles = pairs - even
        for residual in range(even - 2 * forms + 1):
            length = singles + 2 * forms + 2 * residual
            dimension[length] += (1 << singles) * comb(even - 2 * forms, residual)
    assert kernel_count == 1 << (pairs - charge_rank)
    assert sum(dimension) == Fraction(2) ** (pairs - 2 * forms) * admissible
    assert sum(dimension) <= Fraction(2) ** (2 * pairs - 2 * forms - charge_rank)
    assert charge_rank >= forms
    distinct = set(columns) - {0}
    assert (charge_rank == forms) == (len(distinct) == forms)
    return charge_rank, admissible, dimension


def positive_parameter_bound(pairs, forms, charge_rank, z):
    return z ** (2 * forms) * (1 + z * z) ** (charge_rank - 2 * forms) * (1 + z) ** (2 * (pairs - charge_rank))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = bounds = 0
    examples = []
    for forms, largest in ((1, 10), (2, 8), (3, 6)):
        current = 0
        for pairs in range(forms, largest + 1):
            for columns in combinations_with_replacement(range(1 << forms), pairs):
                if rank_two(columns) != forms:
                    continue
                q, admissible, dimensions = exact_dimensions(columns, forms)
                for z in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1), Fraction(4, 3), Fraction(2), Fraction(4)):
                    full_bound = positive_parameter_bound(pairs, forms, q, z)
                    actual = sum(value * z ** degree for degree, value in enumerate(dimensions))
                    assert actual <= full_bound
                    assert all(value * z ** degree <= full_bound for degree, value in enumerate(dimensions))
                    universal = positive_parameter_bound(pairs, forms, forms, z)
                    assert full_bound <= universal
                    bounds += 1
                cases += 1
                current += 1
        examples.append({'forms': forms, 'maximum_pairs': largest, 'full_rank_column_multisets': current})
    records = []
    for columns in ((1, 1, 2, 2, 3, 3), (1, 1, 1, 2, 2, 3, 3),
                    (1, 1, 2, 2, 4, 4), (1, 2, 3, 4, 5, 6, 7)):
        forms = max(columns).bit_length()
        q, admissible, dimensions = exact_dimensions(columns, forms)
        records.append({'columns': list(columns), 'forms': forms, 'charge_rank': q,
                        'admissible_even_subsets': admissible,
                        'homogeneous_dimensions': dimensions,
                        'total_homogeneous_dimension': sum(dimensions)})
    result = {'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'complete_small_multiset_ranges': examples,
              'exact_dimension_and_charge_classification_cases': cases,
              'exact_rational_generating_function_bounds': bounds,
              'examples': records}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
