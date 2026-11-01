#!/usr/bin/env python3
"""Replay incremental atom invariants and matched supplied-prefix costs."""

import argparse
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
from statistics import median
from time import perf_counter, process_time
import tracemalloc

from incremental_oriented_coprime_basis import (
    _ProductTree, IncrementalOrientedBasis, incremental_oriented_coprime_basis,
    verify_incremental_oriented_basis,
)
from oriented_coprime_basis import oriented_coprime_basis, verify_oriented_basis

ROOT = Path(__file__).resolve().parents[1]


def canonical(value):
    return json.loads(json.dumps(value))


def compare(discriminant, keys):
    old = oriented_coprime_basis(discriminant, keys)
    new = incremental_oriented_coprime_basis(discriminant, keys, dense=True)
    assert new['atoms'] == old['atoms'] and new['words'] == old['words']
    verify_oriented_basis(discriminant, keys, new)
    verify_incremental_oriented_basis(canonical(new))
    return new


def tree_controls():
    tree, active, checks = _ProductTree(), {}, 0
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]
    for step in range(48):
        if step % 3 == 2 and active:
            slot = min(active)
            tree.remove(slot)
            del active[slot]
        else:
            used = set(active.values())
            norm = next(p for p in primes if p not in used)
            active[tree.insert(norm)] = norm
        for norm in range(1, 258):
            before = tree.gcd_calls
            actual = sorted(tree.overlaps(norm))
            expected = sorted(slot for slot, value in active.items() if gcd(norm, value) > 1)
            assert actual == expected
            height = tree.capacity.bit_length() - 1
            assert tree.gcd_calls - before <= 1 + 2 * len(expected) * height
            checks += 1
    return {'overlap_queries': checks, 'capacity_doublings': tree.capacity_doublings,
            'leaf_updates': tree.leaf_updates,
            'visited_node_bound': '1 + 2 h log2(capacity)'}


def ideal_controls():
    identities = matched_batches = orientation_batches = 0
    for delta in (-23, -84, -260, -575):
        keys = [(norm, root) for norm in range(1, 100, 2) if gcd(norm, delta) == 1
                for root in range(2 * norm) if (root * root - delta) % (4 * norm) == 0]
        for ordered in (keys, list(reversed(keys)), keys[::2] + keys[1::2]):
            compare(delta, ordered)
            identities += len(ordered)
            matched_batches += 1
        for norm in range(3, 180, 2):
            if gcd(norm, delta) != 1:
                continue
            keys = [(norm, root) for root in range(2 * norm)
                    if (root * root - delta) % (4 * norm) == 0]
            if len(keys) >= 4:
                result = compare(delta, keys)
                assert result['orientation_splits']
                identities += len(keys)
                matched_batches += 1
                orientation_batches += 1
    norms = (9, 27, 243, 169, 1521, 9, 1)
    keys = [(norm, next(root for root in range(2 * norm)
                       if (root * root + 23) % (4 * norm) == 0)) for norm in norms]
    repeated = compare(-23, keys)
    assert [atom['norm'] for atom in repeated['atoms']] == [3, 169]
    identities += len(keys)
    matched_batches += 1
    empty = incremental_oriented_coprime_basis(-23, [])
    assert empty['atoms'] == empty['sparse_words'] == []
    mixed = compare(-23, [(39, 17), (39, 35)])
    assert mixed['words'] == [[-1, -1], [-1, 1]]
    identities += 2
    matched_batches += 1
    invalid = [lambda: IncrementalOrientedBasis(23),
               lambda: IncrementalOrientedBasis(-22),
               lambda: incremental_oriented_coprime_basis(-23, [(2, 1)]),
               lambda: incremental_oriented_coprime_basis(-23, [(23, 23)]),
               lambda: incremental_oriented_coprime_basis(-23, [(3, 3)])]
    for call in invalid:
        try:
            call()
        except ValueError:
            pass
        else:
            raise AssertionError('An invalid request was accepted.')
    corrupt = canonical(mixed)
    corrupt['sparse_words'][0][0][1] *= -1
    try:
        verify_incremental_oriented_basis(corrupt)
    except AssertionError:
        pass
    else:
        raise AssertionError('A wrong orientation certificate was accepted.')
    return {'literal_ideal_identities': identities, 'matched_batches': matched_batches,
            'mixed_orientation_batches': orientation_batches,
            'invalid_inputs_rejected': len(invalid), 'corrupt_certificate_rejected': True,
            'retained_power_atom_norms': [atom['norm'] for atom in repeated['atoms']]}


def prefix_controls(prefix_path):
    source = json.loads(prefix_path.read_text())
    for path, digest in source['source_sha256'].items():
        assert sha256((ROOT / path).read_bytes()).hexdigest() == digest
    cases = []
    for case in source['cases']:
        record = case['result']
        delta, keys = record['refinement']['discriminant'], record['keys']
        state = IncrementalOrientedBasis(delta)
        old_gcd = 0
        checkpoints = []
        for checkpoint in record['checkpoints']:
            count = checkpoint['distinct_residual_keys']
            while len(state.keys) < count:
                state.add(*keys[len(state.keys)])
            current = state.result(dense=True)
            old = oriented_coprime_basis(delta, keys[:count])
            assert current['atoms'] == old['atoms'] and current['words'] == old['words']
            old_gcd += old['decision_gcd_calls']
            checkpoints.append({'queries': checkpoint['queries'], 'keys': count,
                                'incremental_gcd_decisions': current['decision_gcd_calls'],
                                'rebuild_gcd_decisions': old_gcd})
        current = state.result(dense=True)
        assert canonical(current['atoms']) == record['refinement']['atoms']
        assert current['words'] == record['refinement']['words']
        assert old_gcd == record['total_refinement_gcd_decisions']
        verify_oriented_basis(delta, keys, current)
        cases.append({'integer': record['integer'], 'queries': record['norm_queries'],
                      'keys': len(keys), 'atoms': len(current['atoms']),
                      'sparse_entries': current['statistics']['usage_entries'],
                      'dense_entries': len(keys) * len(current['atoms']),
                      'factor_preserved': record['checkpoints'][-1]['new_factor'],
                      'identical_atoms_and_signed_words': True,
                      'checkpoints': checkpoints, 'statistics': current['statistics']})
    return cases, source['cases'][-1]['result']


def matched_benchmark(record, repetitions):
    delta, keys = record['refinement']['discriminant'], record['keys']
    functions = {'rebuild': oriented_coprime_basis,
                 'incremental': incremental_oriented_coprime_basis}
    times = {name: [] for name in functions}
    outputs = {}
    for repetition in range(repetitions):
        order = list(functions) if repetition % 2 == 0 else list(reversed(functions))
        for name in order:
            start, cpu_start = perf_counter(), process_time()
            outputs[name] = functions[name](delta, keys, verify=False)
            times[name].append({'wall_seconds': perf_counter() - start,
                                'cpu_seconds': process_time() - cpu_start})
    results = {}
    for name, function in functions.items():
        tracemalloc.start()
        output = function(delta, keys, verify=False)
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        assert output['atoms'] == outputs[name]['atoms']
        results[name] = {'runs': times[name],
                         'median_wall_seconds': median(t['wall_seconds'] for t in times[name]),
                         'median_cpu_seconds': median(t['cpu_seconds'] for t in times[name]),
                         'tracemalloc_peak_bytes_separate_run': peak,
                         'decision_gcd_calls': output['decision_gcd_calls']}
    assert outputs['rebuild']['atoms'] == outputs['incremental']['atoms']
    sparse = outputs['incremental']['sparse_words']
    assert 'words' not in outputs['incremental']
    for row, word in zip(outputs['rebuild']['words'], sparse):
        assert [(i, e) for i, e in enumerate(row) if e] == word
    return {'integer': record['integer'], 'keys': len(keys), 'repetitions': repetitions,
            'results': results,
            'scope': 'Supplied-basis construction and export only; verification and power generation excluded.',
            'allocation_scope': 'Python allocations from a separate tracemalloc run, not process RSS.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, default=ROOT / 'evidence/oriented-cofactor-prefixes.json')
    parser.add_argument('--repetitions', type=int, default=7)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error('The repetition count must be positive.')
    tree, ideals = tree_controls(), ideal_controls()
    prefixes, fixture = prefix_controls(args.prefix)
    benchmark = matched_benchmark(fixture, args.repetitions)
    sources = ('tools/incremental_oriented_coprime_basis.py',
               'tools/check_incremental_oriented_coprime_basis.py',
               'tools/oriented_coprime_basis.py', 'tools/quadratic_ideal_hnf.py',
               'evidence/oriented-cofactor-prefixes.json')
    report = {'claim': 'Exact sparse incremental oriented ideal refinement',
              'tree_controls': tree, 'ideal_controls': ideals, 'prefixes': prefixes,
              'matched_benchmark': benchmark,
              'source_sha256': {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in sources},
              'scope': 'Supplied-witness implementation improvement; no complete factoring exponent claim.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'tree_controls': tree, 'ideal_controls': ideals,
                      'prefixes': [{k: v for k, v in row.items()
                                    if k not in ('checkpoints', 'statistics')} for row in prefixes],
                      'matched_benchmark': benchmark}, indent=2))


if __name__ == '__main__':
    main()
