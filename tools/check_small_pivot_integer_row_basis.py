#!/usr/bin/env python3
"""Check small internal pivots without changing any selected original relation."""

import argparse
from functools import reduce
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
from random import Random
from statistics import median
from time import perf_counter, process_time

from check_sparse_oriented_power_search import read_data
from cofactor_collision_relations import RationalRowBasis
from primitive_integer_row_basis import PrimitiveIntegerRowBasis
from small_pivot_integer_row_basis import SmallPivotIntegerRowBasis

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('small_pivot_integer_row_basis.py', 'check_small_pivot_integer_row_basis.py',
           'primitive_integer_row_basis.py', 'cofactor_collision_relations.py')


def pivot_invariant(basis):
    previous = []
    for pivot, row in basis.pivots.items():
        assert row[pivot] > 0 and reduce(gcd, row.values(), 0) == 1
        assert all(row.get(column, 0) == 0 for column in previous)
        assert pivot not in previous
        previous.append(pivot)


def compare(rows, width):
    expected, actual = RationalRowBasis(width), SmallPivotIntegerRowBasis(width)
    for row in rows:
        assert expected.add(row) == actual.add(row)
        assert expected.rows == actual.rows
        pivot_invariant(actual)
    return actual


def controls():
    rng = Random(69128)
    matrices = decisions = 0
    for width in range(17):
        for case in range(40):
            rows = []
            for index in range(25):
                row = [rng.randint(-10, 10) if rng.randrange(3) else 0 for _ in range(width)]
                if index % 4 == 3:
                    row = [sum(previous[j] for previous in rows[-2:]) for j in range(width)]
                rows.append(row)
            compare(rows, width)
            matrices += 1
            decisions += len(rows)
    trap = [[2, 0], [0, 6], [1, 3], [4, 12]]
    assert compare(trap, 2).rows == trap[:2]
    permuted = compare([[6, 1, 0], [0, 2, 1], [3, 1, 1], [0, 0, 2]], 3)
    assert list(permuted.pivots) != sorted(permuted.pivots)
    for row, width in (([1], 2), ([1.5], 1), (['1'], 1)):
        try:
            SmallPivotIntegerRowBasis(width).add(row)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid row accepted')
    return {'matrices': matrices, 'acceptance_decisions': decisions,
            'every_original_row_subsequence_matches_rational_elimination': True,
            'insertion_order_triangularity_checked_after_every_row': True,
            'original_relation_content_trap_passed': True,
            'nonmonotone_pivot_column_control_passed': True, 'invalid_rows_rejected': 3}


def benchmark(path, repetitions):
    data = read_data(path)
    record = data['result']
    rows = record['unit_reduction']['torsion_rows']
    width = len(record['unit_reduction']['torsion_generator_indices'])
    expected = record['independent_rows']
    methods = {'leftmost_primitive': PrimitiveIntegerRowBasis,
               'small_pivot_primitive': SmallPivotIntegerRowBasis}
    runs = []
    for repetition in range(repetitions):
        names = list(methods) if repetition % 2 == 0 else list(reversed(methods))
        for name in names:
            tick, cpu = perf_counter(), process_time()
            basis = methods[name](width)
            for row in rows:
                basis.add(row)
            elapsed, consumed = perf_counter() - tick, process_time() - cpu
            assert basis.rows == expected
            if name == 'small_pivot_primitive':
                pivot_invariant(basis)
            runs.append({'repetition': repetition, 'method': name,
                         'elapsed_seconds': elapsed, 'cpu_seconds': consumed,
                         'statistics': basis.statistics})
    return {'bits': record['input_bits'], 'source_file': str(path.resolve().relative_to(ROOT)),
            'source_sha256': sha256(path.read_bytes()).hexdigest(),
            'matrix_rows': len(rows), 'matrix_columns': width, 'rank': len(expected),
            'selected_original_rows_sha256': sha256(json.dumps(
                expected, separators=(',', ':')).encode()).hexdigest(),
            'repetitions': repetitions, 'runs': runs,
            'medians': {name: {metric: median(run[metric] for run in runs if run['method'] == name)
                for metric in ('elapsed_seconds', 'cpu_seconds')} for name in methods},
            'same_selected_original_relations_on_every_run': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, action='append')
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--controls-only', action='store_true')
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error('At least one repetition is required.')
    hashes = {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
              for name in SOURCES}
    report = {'controls': controls(), 'benchmarks': [], 'source_sha256': hashes,
              'scope': 'Exact original-row selection with changed internal pivots; no uniform factoring exponent gain.',
              'timing_scope': 'Rank selection only; complete search and certificate replay costs remain separate.'}
    if not args.controls_only:
        paths = args.reference or [ROOT / 'evidence/optimized-oriented-power-96.json.xz']
        report['benchmarks'] = [benchmark(path, args.repetitions) for path in paths]
    for name, digest in hashes.items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({**report, 'benchmarks': [{key: value for key, value in row.items()
          if key != 'runs'} for row in report['benchmarks']]}, indent=2))


if __name__ == '__main__':
    main()
