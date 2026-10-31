#!/usr/bin/env python3
"""Exact source-oriented Walsh twirl profiles and antipodal-mode controls."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from boolean_cauchy_operator import kernel_matrices, smith_orders
from walsh_cauchy_permanents import (
    principal_energies, principal_permanents, walsh_transform,
)


def valuation(value, precision):
    return (value & -value).bit_length() - 1 if value else precision


def check():
    root = Path(__file__).resolve().parent.parent
    precision = 128
    modulus = 1 << precision
    records = []
    reshuffled_entries = antipodal_controls = determinant_controls = 0
    for length in (2, 4, 8):
        count = 1 << length
        for center in sorted({1, length // 2, length - 1}):
            coordinate = valuation(center, precision)
            first_rows = sum(1 << position for position in range(length)
                             if not (position >> coordinate) & 1)
            second_rows = (count - 1) ^ first_rows
            for point in (1, 2):
                physical, _ = kernel_matrices(length, center, point, precision)
                lambdas = principal_permanents(physical, modulus)
                values = [value % modulus for value in walsh_transform(lambdas)]
                ordinary = [[values[row ^ column] for column in range(count)]
                            for row in range(count)]
                # Build the actual grouping by assigning exactly one port at
                # each position to each side, with opposite symbol orientation.
                grouped = []
                for first in range(count):
                    row = []
                    for second in range(count):
                        original_row = (first & first_rows) | (second & second_rows)
                        original_column = (second & first_rows) | (first & second_rows)
                        row.append(ordinary[original_row][original_column])
                    grouped.append(row)
                assert grouped == ordinary
                reshuffled_entries += count * count
                observed = smith_orders(grouped, precision)
                assert all(order is not None for order in observed)
                eigen_orders = [valuation(value, precision) for value in lambdas]
                assert sum(observed) == length * count + sum(eigen_orders)
                determinant_controls += 1
                assert observed[0] == 0 and all(order > 0 for order in observed[1:])
                assert max(observed) <= 2 * length + max(eigen_orders)
                if point & 1:
                    low_subsets = []
                    half = length // 2
                    for choices in range(1 << half):
                        subset = sum(1 << (position + half * ((choices >> position) & 1))
                                     for position in range(half))
                        absolute, quadratic = principal_energies(length, center, subset)
                        assert absolute == quadratic <= length
                        assert eigen_orders[subset] == absolute
                        low_subsets.append(subset)
                        antipodal_controls += 1
                    assert len(set(low_subsets)) == 1 << half
                    assert sum(order <= 3 * length for order in observed) >= 1 << half
                records.append({
                    'L': length, 'center': center, 'local_point': point,
                    'precision': precision,
                    'source_group_coordinate': coordinate,
                    'grouped_smith_order_multiplicities': dict(sorted(Counter(observed).items())),
                    'principal_permanent_valuation_multiplicities':
                        dict(sorted(Counter(eigen_orders).items())),
                    'full_rank': count,
                    'largest_smith_order': max(observed),
                    'normalization_loss': length,
                })
    sources = [
        'tools/check_walsh_twirl_group.py', 'tools/walsh_cauchy_permanents.py',
        'tools/boolean_cauchy_operator.py', 'research/walsh-principal-permanents.md',
    ]
    return {
        'status': 'pass',
        'independent_source_group_reshuffled_entries': reshuffled_entries,
        'antipodal_balanced_mode_controls': antipodal_controls,
        'smith_determinant_controls': determinant_controls,
        'records': records,
        'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest()
                          for name in sources},
        'scope': 'Exact finite group profiles plus controls for the proved all-size '
                 'antipodal-mode lower bound. Applies to full packet operators and '
                 'the source orientation; no general scalar-decoder impossibility.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = check()
    except Exception as error:
        if args.output:
            args.output.write_text(json.dumps({'status': 'fail', 'error': repr(error)}, indent=2) + '\n')
        raise
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: result[key] for key in (
        'status', 'independent_source_group_reshuffled_entries',
        'antipodal_balanced_mode_controls', 'smith_determinant_controls')}, indent=2))


if __name__ == '__main__':
    main()
