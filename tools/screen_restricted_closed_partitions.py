#!/usr/bin/env python3
"""Screen complete restricted-seed catalogues by adaptive unweighted flags.

Conjugates of one irreducible factor have the same unweighted profiles.
Enumerating the number assigned to the row cut therefore removes only
duplicate assignments, while retaining every possible disjoint partition.
Floating entropy grids are candidate selectors, not exact certificates.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from math import sqrt
from pathlib import Path
import time

from adaptive_mixed_place import adaptive_rate
from search_fused_mixed_weighted import atomic_write


def place_label(factor, index=0):
    if factor is None:
        return None
    factor = tuple(factor)
    return factor[0] if len(factor) == 2 else f'root:{factor}:{index}'


def partitions(record):
    if record['degree'] == 0:
        profiles = record['homogeneous_infinity_orders']
        assert not max(profiles[0]) and not max(profiles[1])
        yield [0], [None], [profiles[0]], [profiles[1]]
        return
    points = [(tuple(point['factor']), point['residue_degree'], point['orders'])
              for point in record['closed_points']]
    points.append((None, 1, record['homogeneous_infinity_orders']))
    options = []
    for factor, multiplicity, profiles in points:
        first, second = bool(max(profiles[0])), bool(max(profiles[1]))
        if first and second:
            options.append(tuple(range(multiplicity+1)))
        elif first:
            options.append((multiplicity,))
        elif second:
            options.append((0,))
        else:
            options.append((None,))
    for allocation in product(*options):
        row, column, row_profiles, column_profiles = [], [], [], []
        for (factor, multiplicity, profiles), count in zip(points, allocation):
            if count is None:
                continue
            row.extend(place_label(factor, i) for i in range(count))
            column.extend(place_label(factor, i) for i in range(count, multiplicity))
            row_profiles.extend([profiles[0]]*count)
            column_profiles.extend([profiles[1]]*(multiplicity-count))
        if row and column:
            yield row, column, row_profiles, column_profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('evidence/finite-fiber-restriction-catalogue.json'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--entropy-steps', type=int, default=4097)
    parser.add_argument('--rate-steps', type=int, default=2049)
    args = parser.parse_args()
    raw = args.input.read_bytes()
    catalogue = json.loads(raw)
    assert catalogue['completed']
    names = (Path(__file__).name, 'adaptive_mixed_place.py', 'search_fused_mixed_weighted.py')
    report = {'scope': __doc__, 'input_sha256': sha256(raw).hexdigest(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'arguments': vars(args) | {'input': str(args.input), 'output': str(args.output)},
              'records': [], 'completed': False, 'completed_partitions': 0}
    started = time.monotonic()
    for index, record in enumerate(catalogue['records']):
        q_x, q_y = len(record['subset_x']), len(record['subset_y'])
        candidates = []
        for row, column, row_profiles, column_profiles in partitions(record):
            value = adaptive_rate(row_profiles, column_profiles, record['degree'], sqrt(q_x*q_y),
                                  args.entropy_steps, args.rate_steps)
            value.update(row_places=row, column_places=column)
            candidates.append(value)
        assert candidates
        result = {key: record[key] for key in ('source_index', 'source_label', 'subset_x', 'subset_y',
                                              'degree', 'tensor_sha256')}
        result.update(catalogue_index=index, paired_dimensions=[q_x, q_y],
                      partition_count=len(candidates), best=min(candidates, key=lambda row: row['rho']),
                      retained_partitions=sorted(candidates, key=lambda row: row['rho'])[:4])
        report['records'].append(result)
        report['completed_partitions'] += len(candidates)
        report['elapsed_seconds'] = time.monotonic()-started
        atomic_write(args.output, report)
        print(json.dumps({'case': index, 'alphabets': [q_x, q_y], 'partitions': len(candidates),
                          'rho': result['best']['rho'], 'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    report.update(completed=True, elapsed_seconds=time.monotonic()-started)
    atomic_write(args.output, report)


if __name__ == '__main__':
    main()
