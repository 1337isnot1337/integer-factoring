#!/usr/bin/env python3
"""Check exact original-row selection and benchmark an actual 96-bit core."""

import argparse
from hashlib import sha256
import json
import lzma
from math import gcd
from functools import reduce
from pathlib import Path
from random import Random
from statistics import median
from time import perf_counter, process_time

from cofactor_collision_relations import RationalRowBasis
from primitive_integer_row_basis import PrimitiveIntegerRowBasis

ROOT = Path(__file__).resolve().parents[1]


def compare(rows, width):
    old, new = RationalRowBasis(width), PrimitiveIntegerRowBasis(width)
    answers = []
    for row in rows:
        first, second = old.add(row), new.add(row)
        assert first == second and old.rows == new.rows
        assert all(reduce(gcd, pivot.values(), 0) == 1 for pivot in new.pivots.values())
        answers.append(first)
    return new, answers


def controls():
    rng, matrices, row_checks = Random(1361907), 0, 0
    for width in range(17):
        for trial in range(40):
            rows = []
            for index in range(2 * width + 9):
                if rows and index % 5 == 0:
                    a, b = rng.sample(range(len(rows)), 2) if len(rows) > 1 else (0, 0)
                    row = [3 * x - 2 * y for x, y in zip(rows[a], rows[b])]
                else:
                    row = [rng.choice((0, 0, 0, -101, -6, -1, 1, 4, 13, 97))
                           for _ in range(width)]
                rows.append(row)
            compare(rows, width)
            matrices += 1
            row_checks += len(rows)
    rows = [[2, 0], [0, 6], [1, 3], [4, 12]]
    basis, answers = compare(rows, 2)
    assert basis.rows == [[2, 0], [0, 6]] and answers == [True, True, False, False]
    assert basis.pivots == {0: {0: 1}, 1: {1: 1}}
    return {'matched_integer_matrices': matrices, 'matched_row_decisions': row_checks,
            'original_nonprimitive_relations_retained': basis.rows,
            'internal_primitive_pivots_not_emitted': True}


def benchmark(path, repetitions):
    source = json.loads(lzma.decompress(path.read_bytes()))
    rows = source['result']['unit_reduction']['torsion_rows']
    width = len(source['result']['unit_reduction']['torsion_generator_indices'])
    functions = {'rational': RationalRowBasis, 'primitive_integer': PrimitiveIntegerRowBasis}
    times = {name: [] for name in functions}
    retained = None
    statistics = None
    for repetition in range(repetitions):
        order = list(functions) if repetition % 2 == 0 else list(reversed(functions))
        for name in order:
            started, cpu_started = perf_counter(), process_time()
            basis = functions[name](width)
            for row in rows:
                basis.add(row)
            times[name].append({'wall_seconds': perf_counter() - started,
                                'cpu_seconds': process_time() - cpu_started})
            if retained is None:
                retained = basis.rows
            assert retained == basis.rows == source['result']['independent_rows']
            if name == 'primitive_integer':
                statistics = basis.statistics
    return {'input_file': str(path.relative_to(ROOT)), 'input_sha256': sha256(path.read_bytes()).hexdigest(),
            'rows': len(rows), 'columns': width, 'independent_original_rows': len(retained),
            'selected_original_rows_sha256': sha256(json.dumps(retained, separators=(',', ':')).encode()).hexdigest(),
            'statistics': statistics,
            'methods': {name: {'runs': value,
                              'median_wall_seconds': median(t['wall_seconds'] for t in value),
                              'median_cpu_seconds': median(t['cpu_seconds'] for t in value)}
                        for name, value in times.items()},
            'scope': 'Original-row selection on a supplied core; no complete factoring timing claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'evidence/optimized-oriented-power-96.json.xz')
    parser.add_argument('--repetitions', type=int, default=5)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error('The repetition count must be positive.')
    control, result = controls(), benchmark(args.input.resolve(), args.repetitions)
    sources = ('primitive_integer_row_basis.py', 'check_primitive_integer_row_basis.py',
               'cofactor_collision_relations.py')
    report = {'claim': 'Primitive temporary rows preserve the exact selected original relations',
              'controls': control, 'benchmark': result,
              'source_sha256': {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in sources}}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'controls': control, 'benchmark': {key: value for key, value in result.items()
                                                       if key != 'methods'},
                     'median_seconds': {name: {'wall': value['median_wall_seconds'],
                                              'cpu': value['median_cpu_seconds']}
                                        for name, value in result['methods'].items()}}, indent=2))


if __name__ == '__main__':
    main()
