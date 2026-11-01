#!/usr/bin/env python3
"""Check changing coprime forests, prepared overlap corrections and ideal words."""
import argparse
from hashlib import sha256
import json
import lzma
from math import gcd, isqrt, prod
from pathlib import Path
from random import Random
from statistics import median
from time import perf_counter, process_time
from batched_coprime_bucket_basis import PreparedCoprimeBucketForest, BufferedBucketOrientedBasis, batched_bucket_oriented_coprime_basis
from coprime_bucket_forest import bucket_oriented_coprime_basis
from incremental_oriented_coprime_basis import incremental_oriented_coprime_basis, verify_incremental_oriented_basis
from quotient_update_oriented_basis import _QuotientUpdateTree, quotient_update_oriented_coprime_basis

ROOT = Path(__file__).resolve().parents[1]
PAYLOAD = ('discriminant', 'keys', 'atoms', 'sparse_words', 'norm_splits', 'orientation_splits', 'input_bit_bound')


def payload(record):
    return json.loads(json.dumps({key: record[key] for key in PAYLOAD}))


def mutable_tree_controls():
    primes = [p for p in range(2, 350) if all(p % d for d in range(2, isqrt(p) + 1))]
    rng = Random(202611010512)
    queries = transitions = product_checks = 0
    for bucket_size in (2, 8, 64):
        for batch_size in (2, 7, 32):
            forest, reference = PreparedCoprimeBucketForest(bucket_size), _QuotientUpdateTree()
            active = {}
            for prime in primes[:20]:
                norm = prime ** rng.randrange(1, 5)
                slot = reference.insert(norm)
                assert forest.insert(norm) == slot
                active[prime] = (slot, norm)
            for _ in range(12):
                norms = [prod(p ** rng.randrange(1, 7) for p in rng.sample(primes, rng.randrange(1, 5)))
                         for _ in range(batch_size)]
                forest.prepare_batch(norms)
                for norm in norms:
                    expected = sorted(slot for slot, value in active.values() if gcd(norm, value) > 1)
                    assert sorted(forest.overlaps(norm)) == expected == sorted(reference.overlaps(norm))
                    queries += 1
                    if active and rng.randrange(3) != 0:
                        prime = rng.choice(sorted(active))
                        slot, _ = active.pop(prime)
                        reference.remove(slot)
                        forest.remove(slot)
                        transitions += 1
                    missing = [prime for prime in primes if prime not in active]
                    if missing:
                        prime = rng.choice(missing)
                        value = prime ** rng.randrange(1, 7)
                        slot = reference.insert(value)
                        assert forest.insert(value) == slot
                        active[prime] = (slot, value)
                        transitions += 1
                    assert forest.products[1] == reference.products[1] == prod(v for _, v in active.values())
                    for tree in forest.buckets:
                        assert all(tree.products[i] == tree.products[2*i] * tree.products[2*i+1]
                                   for i in range(1, tree.capacity))
                    state = forest.prepared
                    for i, tree in enumerate(forest.buckets):
                        old_live = prod(value for slot, value in state['old_slots'].items()
                                        if slot // bucket_size == i)
                        assert gcd(old_live, state['added'][i]) == 1
                        assert old_live * state['added'][i] == tree.products[1]
                    product_checks += 1
                forest.finish_batch()
    return {'prepared_overlap_queries': queries, 'mutable_leaf_transitions': transitions,
            'all_bucket_node_and_delta_checks': product_checks,
            'bucket_sizes': [2, 8, 64], 'batch_sizes': [2, 7, 32]}


def ideal_controls():
    prefixes = identities = 0
    for delta in (-23, -84, -260, -575):
        keys = [(a, b) for a in range(1, 100, 2) if gcd(a, delta) == 1
                for b in range(2*a) if (b*b-delta) % (4*a) == 0]
        keys += keys[:5]
        for ordered in (keys, list(reversed(keys))):
            state = BufferedBucketOrientedBasis(delta, batch_size=7, bucket_size=8)
            for index, (norm, root) in enumerate(ordered):
                assert state.add(norm, root) == index
                assert state.keys == ordered[:index+1]
                if index % 11 == 0 or index + 1 == len(ordered):
                    actual = state.result()
                    expected = incremental_oriented_coprime_basis(delta, ordered[:index+1])
                    assert payload(actual) == payload(expected)
                    assert not state.pending and state.engine.tree.prepared is None
                    prefixes += 1
                    identities += len(actual['keys'])
    return {'ordered_ideal_batches': 8, 'complete_prefix_equalities': prefixes,
            'literal_ideal_identities': identities, 'unit_norms_and_repeated_keys_included': True}


def fixture(bits, repetitions):
    name = f'evidence/sparse-oriented-power-{bits}.json.xz' if bits == 80 else f'evidence/optimized-oriented-power-{bits}.json.xz'
    path = ROOT / name
    saved = json.loads(lzma.decompress(path.read_bytes()))['result']['refinement']
    methods = [
        ('global_quotient_tree', lambda: quotient_update_oriented_coprime_basis(saved['discriminant'], saved['keys'], verify=False)),
        ('local_buckets', lambda: bucket_oriented_coprime_basis(saved['discriminant'], saved['keys'], bucket_size=512, verify=False)),
        ('prepared_buckets', lambda: batched_bucket_oriented_coprime_basis(saved['discriminant'], saved['keys'], verify=False)),
    ]
    samples, final = {name: [] for name, _ in methods}, {}
    for repetition in range(repetitions):
        shift = repetition % len(methods)
        for name, method in methods[shift:] + methods[:shift]:
            started, cpu = perf_counter(), process_time()
            actual = method()
            samples[name].append({'wall_seconds': perf_counter()-started, 'cpu_seconds': process_time()-cpu})
            assert payload(actual) == payload(saved), name
            final[name] = actual
    verification = verify_incremental_oriented_basis(final['prepared_buckets'])
    digest = sha256(json.dumps(payload(saved), sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return {'bits': bits, 'input_path': str(path.relative_to(ROOT)),
            'input_sha256': sha256(path.read_bytes()).hexdigest(), 'keys': len(saved['keys']),
            'atoms': len(saved['atoms']), 'identical_complete_payload_sha256': digest,
            'verification': verification,
            'methods': {name: {'samples': values,
                               'medians': {key: median(item[key] for item in values) for key in ('wall_seconds', 'cpu_seconds')},
                               'statistics': final[name]['statistics']}
                        for name, values in samples.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repetitions', type=int, default=3)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.repetitions < 1:
        parser.error('The repetition count must be positive.')
    report = {'scope': 'Exact supplied-basis acceleration; no bound on relation-search queries.',
              'mutable_trees': mutable_tree_controls(), 'ideals': ideal_controls(),
              'fixtures': [fixture(bits, args.repetitions) for bits in (80, 96)]}
    sources = ('coprime_bucket_forest.py', 'batched_coprime_bucket_basis.py',
               'check_batched_coprime_bucket_basis.py', 'incremental_oriented_coprime_basis.py',
               'quotient_update_oriented_basis.py', 'oriented_coprime_basis.py')
    report['source_sha256'] = {'tools/'+name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in sources}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({**report, 'fixtures': [
        {**item, 'methods': {name: value['medians'] for name, value in item['methods'].items()}}
        for item in report['fixtures']]}, indent=2))


if __name__ == '__main__':
    main()
