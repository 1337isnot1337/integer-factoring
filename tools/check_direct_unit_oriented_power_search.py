#!/usr/bin/env python3
"""Replay complete direct-unit searches on retained matched power prefixes."""

import argparse
from hashlib import sha256
import json
import lzma
from pathlib import Path
from time import process_time

from check_compact_witness_oriented_power_search import SOURCES as COMPACT_SOURCES, compact_witnesses
from check_sparse_oriented_power_search import canonical, read_data
from direct_unit_oriented_power_search import (
    _direct, _reference, factor_with_direct_unit_refinement, verify_direct_unit_search,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCES = COMPACT_SOURCES + ('direct_unit_presentation.py',
    'direct_unit_oriented_power_search.py', 'check_direct_unit_oriented_power_search.py')


def source_hashes():
    return {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
            for name in SOURCES}


def binding_controls():
    assert _direct.__code__ is _reference.__code__
    changed = sorted(name for name, value in _reference.__globals__.items()
                     if _direct.__globals__[name] is not value)
    assert changed == ['eliminate_signed_residual_graph', 'eliminate_unit_residual_columns']
    return {'identical_compact_query_code_object': True,
            'changed_bindings': changed, 'original_namespace_unchanged': True}


def prefix_equalities(old, new):
    fields = ('integer', 'input_bits', 'discriminant', 'cutoff', 'max_queries',
        'checkpoint_interval', 'max_extraction_columns', 'active_primes',
        'prime_records', 'norm_queries', 'trial_divisions', 'keys',
        'old_base_rows', 'candidate_collisions', 'transported_rows',
        'terminal_witness', 'factoring_oracle_calls', 'residual_primality_tests')
    for name in fields:
        assert old[name] == new[name], name
    assert compact_witnesses(old['witnesses']) == new['witnesses']
    excluded = {'statistics', 'decision_gcd_calls', 'tree_update_method'}
    assert {k: v for k, v in old['refinement'].items() if k not in excluded} == {
        k: v for k, v in new['refinement'].items() if k not in excluded}
    assert len(old['checkpoints']) == len(new['checkpoints'])
    for before, after in zip(old['checkpoints'], new['checkpoints']):
        for name in ('queries', 'keys', 'atoms', 'old_base_rank', 'old_factor', 'old_stats'):
            assert before[name] == after[name], name
    return {'same_candidate_prefix': True, 'same_complete_oriented_atom_words': True,
            'same_original_integer_relation_matrix': True,
            'same_original_exact_key_comparison': True,
            'final_presentation_allowed_to_differ': True}


def verify_case(case):
    path = ROOT / case['reference_file']
    assert sha256(path.read_bytes()).hexdigest() == case['reference_sha256']
    reference = read_data(path)
    assert prefix_equalities(reference['result'], case['result']) == case['matched_prefix']
    check = verify_direct_unit_search({'input': reference['input'], 'result': case['result']})
    assert check['units']['original_sparse_matrix_checked']
    assert check['units']['full_integer_presentation_preserved']
    return check


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hashes, bindings = source_hashes(), binding_controls()
    if args.replay:
        report = read_data(args.replay)
        assert report['source_sha256'] == hashes
        for case in report['cases']:
            verify_case(case)
    else:
        cases = []
        for bits in (64, 80):
            path = ROOT / f'evidence/sparse-oriented-power-{bits}.json.xz'
            reference = read_data(path)
            old = reference['result']
            tick = process_time()
            result = factor_with_direct_unit_refinement(old['integer'], old['cutoff'],
                old['max_queries'], checkpoint_interval=old['checkpoint_interval'],
                max_extraction_columns=old['max_extraction_columns'])
            cpu = process_time() - tick
            result = canonical(result)
            case = {'bits': bits, 'reference_file': str(path.relative_to(ROOT)),
                    'reference_sha256': sha256(path.read_bytes()).hexdigest(),
                    'result': result, 'search_cpu_seconds': cpu,
                    'matched_prefix': prefix_equalities(old, result)}
            case['independent_verification'] = verify_case(case)
            cases.append(case)
        report = {'claim': 'Direct unit reduction retains complete exact power-prefix evidence',
                  'binding_controls': bindings, 'cases': cases, 'source_sha256': hashes,
                  'scope': 'Finite complete candidate-prefix comparison and integer presentation replay; no uniform exponent claim.',
                  'timing_scope': 'Search CPU includes every checkpoint and its incidence verifier; excludes JSON conversion and saved replay.'}
        if args.output:
            if args.output.suffix != '.xz':
                parser.error('The complete matched result must end in .xz.')
            with lzma.open(args.output, 'wt', encoding='utf-8', preset=6) as stream:
                json.dump(report, stream, sort_keys=True, separators=(',', ':'))
                stream.write('\n')
    assert source_hashes() == hashes
    print(json.dumps({'binding_controls': bindings,
        'cases': [{'bits': case['bits'], 'queries': case['result']['norm_queries'],
                   'factor': case['result']['factor'],
                   'search_cpu_seconds': case['search_cpu_seconds'],
                   'elapsed_seconds': case['result']['elapsed_seconds'],
                   'matched_prefix': case['matched_prefix'],
                   'unit_verification': case['independent_verification']['units']}
                  for case in report['cases']]}, indent=2))


if __name__ == '__main__':
    main()
