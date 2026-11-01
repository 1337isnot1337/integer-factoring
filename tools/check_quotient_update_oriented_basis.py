#!/usr/bin/env python3
"""Compare exact tree contents and large supplied-prefix construction costs."""

import argparse
from hashlib import sha256
import json
import lzma
from math import gcd
from pathlib import Path
from statistics import median
from time import perf_counter, process_time

from incremental_oriented_coprime_basis import (
    _ProductTree, incremental_oriented_coprime_basis, verify_incremental_oriented_basis,
)
from quotient_update_oriented_basis import (
    _QuotientUpdateTree, quotient_update_oriented_coprime_basis,
)

ROOT = Path(__file__).resolve().parents[1]


def canonical(value):
    return json.loads(json.dumps(value))


def payload(record):
    return {key: record[key] for key in ('discriminant', 'keys', 'atoms', 'sparse_words')}


def tree_controls():
    old, new, slots, checks = _ProductTree(), _QuotientUpdateTree(), [], 0
    for step in range(96):
        if step % 3 == 2 and slots:
            slot = slots.pop(0)
            old.remove(slot)
            new.remove(slot)
        else:
            norm = (1 << (5 + step % 60)) + 2 * step + 1
            first, second = old.insert(norm), new.insert(norm)
            assert first == second
            slots.append(first)
        assert old.products == new.products and old.leaves == new.leaves
        assert old.capacity == new.capacity and old.free == new.free
        for norm in range(1, 128):
            assert old.overlaps(norm) == new.overlaps(norm)
            checks += 1
    return {'exact_tree_state_comparisons': 96, 'overlap_comparisons': checks,
            'capacity_doublings': old.capacity_doublings,
            'exact_divisions': new.exact_divisions}


def ideal_controls():
    identities = 0
    for delta in (-23, -84, -260, -575):
        keys = [(a, b) for a in range(1, 100, 2) if gcd(a, delta) == 1
                for b in range(2 * a) if (b * b - delta) % (4 * a) == 0]
        for ordered in (keys, list(reversed(keys))):
            old = incremental_oriented_coprime_basis(delta, ordered)
            new = quotient_update_oriented_coprime_basis(delta, ordered)
            assert payload(old) == payload(new)
            assert old['decision_gcd_calls'] == new['decision_gcd_calls']
            identities += len(keys)
    return {'literal_ideal_identities': identities, 'identical_refinement_batches': 8}


def benchmark(path, repetitions):
    data = json.loads(lzma.decompress(path.read_bytes()))
    supplied = data['result']['refinement']
    delta, keys = supplied['discriminant'], supplied['keys']
    functions = {'ancestor_products': incremental_oriented_coprime_basis,
                 'exact_quotients': quotient_update_oriented_coprime_basis}
    times, results = {name: [] for name in functions}, {}
    for repetition in range(repetitions):
        order = list(functions) if repetition % 2 == 0 else list(reversed(functions))
        for name in order:
            start, cpu_start = perf_counter(), process_time()
            result = functions[name](delta, keys, verify=False)
            times[name].append({'wall_seconds': perf_counter() - start,
                                'cpu_seconds': process_time() - cpu_start})
            assert canonical(payload(result)) == payload(supplied)
            results[name] = result
    verification = verify_incremental_oriented_basis(results['exact_quotients'])
    assert payload(results['ancestor_products']) == payload(results['exact_quotients'])
    digest = sha256(json.dumps(payload(results['exact_quotients']), sort_keys=True,
                              separators=(',', ':')).encode()).hexdigest()
    return {'input_path': str(path.relative_to(ROOT)), 'input_sha256': sha256(path.read_bytes()).hexdigest(),
            'integer': data['result']['integer'], 'keys': len(keys),
            'atoms': len(supplied['atoms']), 'repetitions': repetitions,
            'identical_mathematical_payload_sha256': digest, 'verification': verification,
            'methods': {name: {'runs': times[name],
                              'median_wall_seconds': median(t['wall_seconds'] for t in times[name]),
                              'median_cpu_seconds': median(t['cpu_seconds'] for t in times[name]),
                              'statistics': results[name]['statistics']}
                        for name in functions},
            'scope': 'Supplied-basis construction and sparse export; no complete-search timing claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, default=ROOT / 'evidence/sparse-oriented-power-80.json.xz')
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error('The repetition count must be positive.')
    controls, ideals = tree_controls(), ideal_controls()
    result = benchmark(args.prefix.resolve(), args.repetitions)
    sources = ('quotient_update_oriented_basis.py', 'check_quotient_update_oriented_basis.py',
               'incremental_oriented_coprime_basis.py', 'oriented_coprime_basis.py')
    report = {'claim': 'Exact quotient updates preserve every atom and improve the measured component cost',
              'tree_controls': controls, 'ideal_controls': ideals, 'benchmark': result,
              'source_sha256': {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in sources}}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'tree_controls': controls, 'ideal_controls': ideals,
                     'benchmark': {k: v for k, v in result.items() if k != 'methods'},
                     'medians': {name: {'wall': value['median_wall_seconds'],
                                       'cpu': value['median_cpu_seconds']}
                                 for name, value in result['methods'].items()}}, indent=2))


if __name__ == '__main__':
    main()
