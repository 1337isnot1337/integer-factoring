#!/usr/bin/env python3
"""Exact dimensions for a width-two overlapping Ward family."""

import argparse
import hashlib
import json
from math import comb, log2, sqrt
from pathlib import Path

from check_modular_crossing import rank_two
from check_ward_interval_transport import interval_profile, minimal_span_basis


def fan_types(forms):
    result, edges = [], []
    for i in range(forms):
        result.append(1 << i)
        edges.append((0, i + 1))
        if i + 1 < forms:
            result.append((1 << i) | (1 << (i + 1)))
            edges.append((i + 1, i + 2))
    return result, edges


def count_polynomial(forms):
    connected, pending = [0, 1], [1]
    for _ in range(1, forms):
        size = len(connected) + 2
        first, second = [0] * size, [0] * size
        for i, value in enumerate(connected):
            first[i + 1] += 2 * value
            first[i + 2] += value
            second[i] += value
        for i, value in enumerate(pending):
            first[i + 2] += value
            second[i + 1] += value
        connected, pending = first, second
    return connected


def charge(column, forms):
    result, bit = 0, 0
    for a in range(forms):
        for b in range(a, forms):
            result |= ((column >> a & 1) * (column >> b & 1)) << bit
            bit += 1
    return result


def family_record(forms, exhaustive=False):
    types, edges = fan_types(forms)
    columns = [column for column in types for _ in range(2)]
    pairs = len(columns)
    rows = [sum(((column >> a) & 1) << i for i, column in enumerate(columns)) for a in range(forms)]
    assert rank_two(columns) == forms
    assert rank_two([charge(column, forms) for column in types]) == 2 * forms - 1
    total_gram = 0
    for column in columns:
        total_gram ^= charge(column, forms)
    assert total_gram == 0
    optimized, _, _ = minimal_span_basis(rows, pairs)
    width = 1 if forms == 1 else 2
    assert max(interval_profile(rows, pairs)) == max(interval_profile(optimized, pairs)) == width
    counts = count_polynomial(forms)
    subsets = 0
    if exhaustive:
        direct = [0] * len(counts)
        for mask in range(1 << len(types)):
            selected = [i for i in range(len(types)) if mask >> i & 1]
            full = rank_two([types[i] for i in selected]) == forms
            reached = {0}
            while True:
                expanded = reached | {b for i in selected for a, b in (edges[i], edges[i][::-1]) if a in reached}
                if expanded == reached:
                    break
                reached = expanded
            assert full == (len(reached) == forms + 1)
            if full:
                direct[len(selected)] += 1
            subsets += 1
        assert direct == counts
    dimensions = [0] * (2 * pairs + 1)
    for selected, count in enumerate(counts):
        if not count:
            continue
        assert forms <= selected <= 2 * forms - 1
        singles = pairs - 2 * selected
        for residual in range(2 * selected - 2 * forms + 1):
            length = singles + 2 * forms + 2 * residual
            dimensions[length] += count * (1 << singles) * comb(2 * selected - 2 * forms, residual)
    total = (1 << (pairs - 2 * forms)) * sum(counts)
    assert sum(dimensions) == total
    central = dimensions[pairs]
    assert central * (2 * forms - 1) >= total
    zero_syndrome = comb(2 * forms - 2, forms - 1)
    assert (central - zero_syndrome) % 2 == 0
    return {'forms': forms, 'pairs': pairs, 'columns': columns, 'quadratic_charge_rank': 2 * forms - 1,
            'optimal_interval_width': width, 'operator_channel_bound': 8 * 2 ** width,
            'full_rank_type_subset_counts_by_size': counts,
            'connected_spanning_subgraphs': sum(counts), 'type_subsets_exhausted': subsets,
            'total_homogeneous_dimension': total, 'central_grade_dimension': central,
            'zero_syndrome_central_dimension': zero_syndrome,
            'orientation_selected_central_labels': (central - zero_syndrome) // 2,
            'homogeneous_dimensions': dimensions}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records = [family_record(k, exhaustive=k <= 7) for k in range(1, 13)]
    sequence = [0] + [row['connected_spanning_subgraphs'] for row in records]
    assert sequence[1:3] == [1, 4]
    assert all(sequence[k] == 4 * sequence[k - 1] - 2 * sequence[k - 2] for k in range(2, len(sequence)))
    result = {'source_sha256': {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                               for name in (Path(__file__).name, 'check_ward_interval_transport.py')},
              'cases': records, 'asymptotic_label_bits_per_pair': (2 + log2(2 + sqrt(2))) / 4,
              'scope': 'Exact construction and label counts; no cut-rank or factoring-score improvement asserted.'}
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
