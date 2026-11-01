#!/usr/bin/env python3
"""Independent local-isomorphism, full-map, and compact-certificate controls."""
import argparse
from copy import deepcopy
from hashlib import sha256
from itertools import product
import json
import lzma
from pathlib import Path
from random import Random
from statistics import median
from time import perf_counter, process_time

from signed_residual_graph import eliminate_signed_residual_graph, verify_signed_graph_elimination
from sparse_unit_presentation import eliminate_sparse_unit_residual_columns
from unit_presentation_trace import eliminate_unit_residual_columns_compact, verify_unit_presentation_trace

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ('output_rows', 'new_generator_indices', 'free_root_columns', 'torsion_columns',
          'torsion_generator_indices', 'torsion_rows', 'counts')


def canonical(value):
    return json.loads(json.dumps(value))


def encoded(value):
    return json.dumps(canonical(value), sort_keys=True, separators=(',', ':')).encode()


def compare(full, compact, graph):
    for name in FIELDS:
        assert canonical(full[name]) == canonical(compact[name]), name
    assert [[graph['new_generator_indices'][j], s] for _, j, s in compact['pivots']] == [
        [step['old_generator_index'], step['pivot_coefficient']] for step in full['unit_elimination_steps']]


def character_count(rows, width, modulus):
    return sum(all(sum(a * b for a, b in zip(row, word)) % modulus == 0 for row in rows)
               for word in product(range(modulus), repeat=width))


def small_controls():
    rng = Random(20261101091)
    pivots = characters = 0
    digest = sha256()
    for case in range(432):
        base, residuals = rng.randrange(4), rng.randrange(7)
        width = base + residuals
        rows = [[rng.choice((-3, -2, -1, -1, 0, 0, 0, 1, 1, 2, 3))
                 for _ in range(width)] for _ in range(rng.randrange(16))]
        graph = eliminate_signed_residual_graph(rows, base, generator_count=width, allow_extra_rows=True)
        graph_before = encoded(graph)
        full = eliminate_sparse_unit_residual_columns(rows, graph)
        compact = eliminate_unit_residual_columns_compact(rows, graph)
        compare(full, compact, graph)
        verify_unit_presentation_trace(rows, canonical(graph), canonical(compact))
        assert encoded(graph) == graph_before
        pivots += len(compact['pivots'])
        if case < 100 and width <= 5:
            for modulus in (2, 4):
                expected = character_count(rows, width, modulus)
                actual = character_count(compact['torsion_rows'], len(compact['torsion_columns']), modulus)
                actual *= modulus ** len(compact['free_root_columns'])
                assert expected == actual
                characters += 1
        digest.update(encoded([rows, compact]))
    return {'presentations': 432, 'complete_projection_and_pivot_equalities': 432,
            'local_unit_isomorphisms': pivots, 'independent_character_counts': characters,
            'json_round_trips': 432, 'input_graphs_unchanged': True,
            'result_sha256': digest.hexdigest()}


def rejection_controls():
    rows = [[1, 1, 1, 1], [0, 2, 1, -1]]
    graph = eliminate_signed_residual_graph(rows, 1, allow_extra_rows=True)
    good = eliminate_unit_residual_columns_compact(rows, graph)
    assert good['pivots']
    corrupt = []
    item = deepcopy(good); item['pivots'][0][2] *= -1; corrupt.append(item)
    item = deepcopy(good); item['pivots'][0][1] = 0; corrupt.append(item)
    item = deepcopy(good); item['pivots'][0][0] = 1000; corrupt.append(item)
    item = deepcopy(good); item['pivots'].pop(); corrupt.append(item)
    item = deepcopy(good); item['pivots'].append(item['pivots'][0]); corrupt.append(item)
    item = deepcopy(good); item['input_presentation_sha256'] = '0' * 64; corrupt.append(item)
    item = deepcopy(good); item['new_generator_indices'][0] += 1; corrupt.append(item)
    item = deepcopy(good); item['torsion_rows'].append([1]); corrupt.append(item)
    for item in corrupt:
        try:
            verify_unit_presentation_trace(rows, graph, item)
        except (AssertionError, ValueError, IndexError, KeyError):
            continue
        raise AssertionError('A corrupted unit trace was accepted.')
    return {'corrupted_traces_rejected': len(corrupt)}


def retained_control(bits, repeats):
    path = ROOT / 'evidence' / f'sparse-oriented-power-{bits}.json.xz'
    record = json.loads(lzma.decompress(path.read_bytes()))['result']
    graph, expected = record['graph_reduction'], record['unit_reduction']
    rows = [{int(j): value for j, value in row.items()} for row in record['transported_rows']]
    verify_signed_graph_elimination(rows, graph)
    compact = eliminate_unit_residual_columns_compact(rows, graph)
    compare(expected, compact, graph)
    samples = {'expanded_maps': [], 'compact_trace': []}
    if repeats:
        methods = [('expanded_maps', lambda: eliminate_sparse_unit_residual_columns(rows, graph)),
                   ('compact_trace', lambda: eliminate_unit_residual_columns_compact(
                       rows, graph, verify_graph=False))]
        for iteration in range(repeats):
            for name, method in (methods if iteration % 2 == 0 else methods[::-1]):
                started, cpu = perf_counter(), process_time()
                out = method()
                samples[name].append({'wall_seconds': perf_counter()-started, 'cpu_seconds': process_time()-cpu})
                if name == 'expanded_maps':
                    compare(out, compact, graph)
                else:
                    compare(expected, out, graph)
    result = {'bits': bits, 'input_sha256': sha256(path.read_bytes()).hexdigest(),
              'original_rows': len(rows), 'complete_reference_projection_equal': True,
              'expanded_json_bytes': len(encoded(expected)),
              'compact_json_bytes': len(encoded(compact)),
              'expanded_xz_bytes': len(lzma.compress(encoded(expected))),
              'compact_xz_bytes': len(lzma.compress(encoded(compact))),
              'compact_certificate': compact}
    if repeats:
        result['timing'] = {'samples': samples,
                            'medians': {name: {key: median(sample[key] for sample in values)
                                               for key in ('wall_seconds', 'cpu_seconds')}
                                        for name, values in samples.items()},
                            'scope': 'Search component after its graph is verified; full maps or local trace verified inside timing.'}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--benchmark-repeats', type=int, default=0)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.benchmark_repeats < 0:
        parser.error('Benchmark repeats must be nonnegative.')
    report = {'scope': 'Composed integer presentation certificates, with no new relation-search bound.',
              'small': small_controls(), 'rejections': rejection_controls(),
              'retained': [retained_control(bits, args.benchmark_repeats) for bits in (64, 80)]}
    sources = ('unit_presentation_trace.py', 'check_unit_presentation_trace.py',
               'sparse_unit_presentation.py', 'signed_residual_graph.py')
    report['source_sha256'] = {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                               for name in sources}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({**report, 'retained': [
        {key: value for key, value in item.items() if key != 'compact_certificate'}
        for item in report['retained']]}, indent=2))


if __name__ == '__main__':
    main()
