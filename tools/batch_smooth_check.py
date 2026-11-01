#!/usr/bin/env python3
"""Exact controls and opt-in matched benchmarks for batched smooth parts."""

import argparse
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import platform
from random import Random
from time import perf_counter, process_time

from batch_smooth_parts import (product_tree, simultaneous_remainders,
                               smooth_flags, smooth_parts,
                               smooth_sequence_batches, smoothness_exponent)
from check_leading_power_smooth_relations import verify_witnesses, canonical
from check_prime_cofactor_saturation import verify_result
from leading_power_smooth_relations import factor_with_smooth_leading_powers
from prime_cofactor_saturation import factor_with_prime_cofactor_saturation


ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('tools/batch_smooth_parts.py', 'tools/batch_smooth_check.py',
           'tools/leading_power_smooth_relations.py',
           'tools/prime_cofactor_saturation.py')
FIXTURES = {
    92: (4951760154835678088235319297, 1021, 10000),
    124: (21267647932558653899591465697288388633, 4093, 100000),
    144: (17355670471043633865640786569418441084204301, 8191, 200000),
}


def payload(record, kind):
    if kind == 'fixed':
        fields = ('factor', 'method', 'witnesses', 'relation_matrix',
                  'trial_divisions', 'terminal_witness', 'prime_records',
                  'active_primes', 'last_extraction')
    else:
        fields = ('factor', 'method', 'witnesses', 'generators', 'trial_divisions',
                  'sparse_augmented_relation_matrix', 'fixed_relation_sources',
                  'fixed_relation_sha256', 'terminal_witness', 'last_extraction',
                  'saturation_prime_count', 'maximum_saturation_product_bits')
    return {key: record[key] for key in fields}


def payload_digest(record, kind):
    encoded = json.dumps(payload(record, kind), sort_keys=True, separators=(',', ':'))
    return sha256(encoded.encode()).hexdigest()


def compare(reference, batched, kind, batch_size):
    assert payload(reference, kind) == payload(batched, kind)
    assert batched['scanned_norm_queries'] == reference['norm_queries']
    assert batched['norm_queries'] >= reference['norm_queries']
    assert batched['generated_ahead_norms'] == \
        batched['norm_queries'] - reference['norm_queries']
    searches = (len(batched['attempts']) if kind == 'fixed' else
                batched['initial_power_computations'])
    assert batched['generated_ahead_norms'] <= searches * (batch_size - 1)
    assert batched['modular_smoothness_tests'] == batched['norm_queries']
    if kind == 'saturation':
        assert batched['norm_queries'] <= 2 * batched['scanned_norm_queries']
        assert (batched['successive_power_compositions'] -
                reference['successive_power_compositions']) == batched['generated_ahead_norms']
    else:
        assert sum(item['queries'] for item in batched['attempts']) == batched['norm_queries']
        assert sum(item['scanned_queries'] for item in batched['attempts']) == \
            batched['scanned_norm_queries']


def kernel_controls():
    remainder_cases = smooth_part_cases = 0
    stopping_controls = 0
    lists = [(), (1,), (1, 1, 1), tuple(range(1, 100)),
             (12, 18, 24, 36, 72, 12, 1), (2**80, 3**45, 2**80, 5**32)]
    rng = Random(20261101)
    lists += [tuple(rng.randrange(1, 1 << bits) for _ in range(length))
              for bits in (4, 32, 97, 256) for length in (3, 8, 31, 64)]
    for values in lists:
        tree = product_tree(values)
        assert not values or tree[-1] == [prod(values)]
        for dividend in (0, 1, 17, (1 << 1000) - 1):
            assert simultaneous_remainders(dividend, values) == [dividend % a for a in values]
            remainder_cases += len(values)
    for primes in ((), (2,), (2, 3, 5), (3, 7, 13), (2, 5, 11)):
        product = prod(primes)
        values = tuple(range(1, 1201)) + (2**80, 3**45, 5**32, 2**30 * 7**18)
        expected = []
        for value in values:
            residual = value
            for prime in primes:
                while residual % prime == 0:
                    residual //= prime
            expected.append(value // residual)
        assert smooth_parts(values, product) == expected
        assert smooth_parts(values, product**3) == expected
        assert smooth_flags(values, product) == [a == s for a, s in zip(values, expected)]
        smooth_part_cases += len(values)
    invalid = [lambda: product_tree([0]), lambda: simultaneous_remainders(-1, [2]),
               lambda: smooth_parts([2], 0), lambda: smooth_parts([128], 2, exponent=2),
               lambda: smoothness_exponent(0),
               lambda: list(smooth_sequence_batches(1, lambda x: x + 1, 1, 2, 0))]
    for call in invalid:
        try:
            call()
        except ValueError:
            pass
        else:
            raise AssertionError('An invalid arithmetic request was accepted.')
    for count in range(50):
        for batch in (1, 2, 3, 8, 17):
            blocks = list(smooth_sequence_batches(1, lambda x: x + 1, count, 30, batch))
            assert [x for _, block, _ in blocks for x in block] == list(range(1, count + 1))
            assert all(1 <= len(block) <= batch for _, block, _ in blocks)
            assert [s for _, _, parts in blocks for s in parts] == smooth_parts(range(1, count + 1), 30)
            for offset, block, _ in blocks:
                for first_stop in range(offset + 1, offset + len(block) + 1):
                    assert offset + len(block) <= min(first_stop + batch - 1, 2 * first_stop - 1)
                    stopping_controls += 1
    return {'simultaneous_remainders': remainder_cases,
            'smooth_parts_and_flags': smooth_part_cases,
            'invalid_inputs_rejected': len(invalid), 'sequence_controls': 250,
            'geometric_first_stop_bounds': stopping_controls}


def small_controls():
    matched = budget_controls = 0
    cases = ((15, 2, 20), (35, 5, 20), (77, 3, 20), (539, 3, 30),
             (49, 3, 20), (21, 2, 20), (97, 3, 20), (1022117, 7, 1000),
             (1000036000099, 31, 5000), (27221, 5, 100))
    for case in cases:
        saturation_reference = None
        for kind, function, verifier in (
                ('fixed', factor_with_smooth_leading_powers, lambda r: verify_witnesses(canonical(r))),
                ('saturation', factor_with_prime_cofactor_saturation, verify_result)):
            reference = function(*case)
            if kind == 'saturation':
                saturation_reference = reference
            verifier(reference)
            for batch in (2, 8, 32, 128):
                result = function(*case, smooth_batch_size=batch)
                compare(reference, result, kind, batch)
                matched += 1
        for cap in (0, 1, 2, 3, 4, 7, 16, 31, 64, 127):
            result = factor_with_prime_cofactor_saturation(*case, max_queries=cap,
                                                           smooth_batch_size=32)
            verify_result(result)
            assert result['norm_queries'] <= cap
            assert result['witnesses'] == \
                saturation_reference['witnesses'][:len(result['witnesses'])]
            budget_controls += 1
    return {'matched_first_success_runs': matched, 'strict_query_budget_runs': budget_controls}


def benchmark(bits, kind, batches):
    case = FIXTURES[bits]
    reference = None
    output = []
    for batch in (1, *batches):
        start = perf_counter()
        cpu_start = process_time()
        if kind == 'fixed':
            result = factor_with_smooth_leading_powers(*case, smooth_batch_size=batch)
        else:
            result = factor_with_prime_cofactor_saturation(
                *case, smooth_batch_size=batch, verify_candidates=True, include_relations=False)
        elapsed = perf_counter() - start
        cpu_elapsed = process_time() - cpu_start
        if reference is None:
            reference = result
            if kind == 'fixed':
                verify_witnesses(canonical(result))
            else:
                verify_result(result)
        else:
            compare(reference, result, kind, batch)
        record = {'input_bits': bits, 'kind': kind, 'batch_size': batch,
                  'seconds': elapsed, 'cpu_seconds': cpu_elapsed, 'factor': result['factor'],
                  'norm_queries': result['norm_queries'],
                  'scanned_norm_queries': result.get('scanned_norm_queries', result['norm_queries']),
                  'generated_ahead_norms': result.get('generated_ahead_norms', 0),
                  'relation_rank': result['relation_rank'],
                  'prime_product_bits': result.get('maximum_saturation_product_bits',
                                                    result.get('prime_product_bits')),
                  'full_mathematical_payload_sha256': payload_digest(result, kind),
                  'matches_reference_witnesses_and_relations': True}
        output.append(record)
        print(json.dumps(record), flush=True)
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--benchmark', type=int, choices=tuple(FIXTURES))
    parser.add_argument('--kind', choices=('fixed', 'saturation'), default='saturation')
    parser.add_argument('--batches', type=int, nargs='+', default=[32])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if any(value <= 1 for value in args.batches):
        parser.error('Comparison batch sizes must exceed one.')
    report = {'scope': 'Exact arithmetic and finite matched traversals; no uniform window or score improvement.',
              'algorithm_source': 'D. J. Bernstein, How to find smooth parts of integers (2004), Algorithm 2.1',
              'kernel_controls': kernel_controls(), 'selector_controls': small_controls(),
              'source_sha256': {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES},
              'environment': {'python': platform.python_version(), 'machine': platform.machine(),
                              'system': platform.system()}}
    if args.benchmark:
        report['matched_benchmark'] = benchmark(args.benchmark, args.kind, args.batches)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key != 'matched_benchmark'}, indent=2))


if __name__ == '__main__':
    main()
