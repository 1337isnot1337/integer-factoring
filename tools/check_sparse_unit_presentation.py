#!/usr/bin/env python3
"""Replay sparse unit elimination against integer maps and the reference pass."""
import argparse
from hashlib import sha256
from itertools import product
import json
import lzma
from pathlib import Path
from random import Random
from statistics import median
from time import perf_counter

from signed_residual_graph import eliminate_signed_residual_graph, verify_signed_graph_elimination
from signed_residual_unit_elimination import eliminate_unit_residual_columns
from sparse_unit_presentation import eliminate_sparse_unit_residual_columns

ROOT = Path(__file__).resolve().parents[1]


def canonical(value):
    return json.dumps(json.loads(json.dumps(value)), sort_keys=True, separators=(',', ':'))


def characters(rows, width, modulus):
    return sum(all(sum(a * b for a, b in zip(row, vector)) % modulus == 0 for row in rows)
               for vector in product(range(modulus), repeat=width))


def small_controls():
    rng = Random(20261101073)
    comparisons = pivots = character_checks = 0
    digest = sha256()
    explicit = [
        (1, 5, [[2, 1, 1, 1, 0], [3, -1, 2, 0, 1], [0, 0, 0, 0, 0]]),
        (0, 4, [[1, 1, 1, 0], [-1, -1, -1, 0], [0, 2, 2, 0]]),
        (0, 3, [[2, 0, 0], [0, 2, 0], [1, 1, -1]]),
        (2, 5, [[3, 0, 1, -1, 1], [-2, 1, 0, 1, 1], [1, 1, 1, 0, 2]]),
        (0, 5, []), (2, 4, [[0, 0, 0, 0]]),
    ]
    for trial in range(606):
        if trial < len(explicit):
            base, width, rows = explicit[trial]
        else:
            base, residuals = rng.randrange(4), rng.randrange(8)
            width = base + residuals
            rows = []
            for _ in range(rng.randrange(20)):
                row = [rng.randrange(-5, 6) for _ in range(base)] + [0] * residuals
                if residuals:
                    if rng.randrange(3) == 0:
                        for _ in range(rng.randrange(3)):
                            row[base + rng.randrange(residuals)] += rng.choice((-1, 1))
                    else:
                        for j in range(base, width):
                            row[j] = rng.choice((-3, -2, -1, -1, 0, 0, 0, 1, 1, 2, 3))
                rows.append(row)
        sparse = [{j: value for j, value in enumerate(row) if value} for row in rows]
        graph = eliminate_signed_residual_graph(sparse, base, generator_count=width,
                                                allow_extra_rows=True)
        before = canonical(graph)
        expected = eliminate_unit_residual_columns(sparse, graph)
        actual = eliminate_sparse_unit_residual_columns(sparse, graph)
        assert actual == expected
        assert canonical(graph) == before
        assert eliminate_sparse_unit_residual_columns(rows, json.loads(before)) == expected
        verify_signed_graph_elimination(sparse, json.loads(canonical(actual)))
        comparisons += 2
        pivots += len(actual['unit_elimination_steps'])
        if trial < 120 and width <= 5:
            for modulus in (3, 4):
                original = characters(rows, width, modulus)
                reduced = characters(actual['torsion_rows'], len(actual['torsion_columns']), modulus)
                reduced *= modulus ** len(actual['free_root_columns'])
                assert original == reduced
                character_checks += 1
        digest.update(canonical([rows, actual]).encode())
    return {'presentations': 606, 'complete_reference_equalities': comparisons,
            'unit_pivots': pivots, 'independent_character_counts': character_checks,
            'input_presentations_unchanged': True, 'result_sha256': digest.hexdigest()}


def retained_control(bits, repeats):
    path = ROOT / 'evidence' / f'sparse-oriented-power-{bits}.json.xz'
    document = json.loads(lzma.decompress(path.read_bytes()))
    record = document['result']
    rows = [{int(j): value for j, value in row.items()} for row in record['transported_rows']]
    expected = canonical(record['unit_reduction'])
    samples = {'reference_seconds': [], 'sparse_seconds': []}
    actual = None
    if repeats:
        methods = [('reference_seconds', eliminate_unit_residual_columns),
                   ('sparse_seconds', eliminate_sparse_unit_residual_columns)]
        for iteration in range(repeats):
            for name, method in (methods if iteration % 2 == 0 else methods[::-1]):
                started = perf_counter()
                out = method(rows, record['graph_reduction'], verify=True)
                samples[name].append(perf_counter() - started)
                assert canonical(out) == expected
                actual = out
    else:
        actual = eliminate_sparse_unit_residual_columns(rows, record['graph_reduction'], verify=True)
        assert canonical(actual) == expected
    result = {'bits': bits, 'input_sha256': sha256(path.read_bytes()).hexdigest(),
              'original_rows': len(rows), 'old_generators': actual['old_generator_count'],
              'unit_pivots': len(actual['unit_elimination_steps']),
              'final_rows': len(actual['output_rows']),
              'nonfree_generators': len(actual['torsion_columns']),
              'complete_retained_certificate_equal': True,
              'certificate_sha256': sha256(expected.encode()).hexdigest(),
              'verification': actual['verification']}
    if repeats:
        result['timing'] = {**samples, 'reference_median_seconds': median(samples['reference_seconds']),
                            'sparse_median_seconds': median(samples['sparse_seconds']),
                            'scope': 'Matched final pass only; includes copying and full certificate verification.'}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark-repeats', type=int, default=0)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.benchmark_repeats < 0:
        parser.error('Benchmark repeats must be nonnegative.')
    report = {'scope': 'Exact implementation acceleration of supplied presentations; no search or exponent bound.',
              'small': small_controls(),
              'retained': [retained_control(bits, args.benchmark_repeats) for bits in (64, 80)]}
    sources = ('sparse_unit_presentation.py', 'check_sparse_unit_presentation.py',
               'signed_residual_graph.py', 'signed_residual_unit_elimination.py')
    report['source_sha256'] = {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                               for name in sources}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
