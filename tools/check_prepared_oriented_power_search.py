#!/usr/bin/env python3
"""Check the complete prepared search against a retained reference certificate."""

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from time import process_time

from check_compact_oriented_power_search import (
    SOURCES as COMPACT_SOURCES, apply_patch as apply_compact_patch,
    make_patch as make_compact_patch,
)
from check_sparse_oriented_power_search import canonical, read_data
from prepared_oriented_power_search import (
    _prepared, _reference, factor_with_prepared_oriented_refinement, verify_prepared_search,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCES = COMPACT_SOURCES + ('prepared_oriented_power_search.py',
    'check_prepared_oriented_power_search.py', 'small_pivot_integer_row_basis.py',
    'coprime_bucket_forest.py', 'batched_coprime_bucket_basis.py')


def binding_controls():
    assert _prepared.__code__ is _reference.__code__
    assert _prepared.__globals__ is not _reference.__globals__
    changed = sorted(name for name, value in _reference.__globals__.items()
                     if _prepared.__globals__[name] is not value)
    assert changed == ['IncrementalOrientedBasis', 'RationalRowBasis',
                       'eliminate_unit_residual_columns']
    assert _prepared.__kwdefaults__ == _reference.__kwdefaults__
    return {'identical_search_code_object': True,
            'changed_dependency_bindings': changed,
            'original_global_namespace_unchanged': True}


def set_gcd_counters(reference, patch):
    result = deepcopy(reference)
    assert len(result['checkpoints']) == len(patch['checkpoint_gcd_decisions'])
    for row, count in zip(result['checkpoints'], patch['checkpoint_gcd_decisions']):
        row['basis_gcd_decisions'] = count
    result['refinement']['decision_gcd_calls'] = patch['refinement_gcd_decisions']
    return result


def apply_patch(reference, patch):
    assert set(patch) == {'checkpoint_gcd_decisions', 'refinement_gcd_decisions', 'compact_patch'}
    return apply_compact_patch(set_gcd_counters(reference, patch), patch['compact_patch'])


def make_patch(reference, new):
    patch = {'checkpoint_gcd_decisions': [row['basis_gcd_decisions'] for row in new['checkpoints']],
             'refinement_gcd_decisions': new['refinement']['decision_gcd_calls']}
    patch['compact_patch'] = make_compact_patch(set_gcd_counters(reference, patch), new)
    assert apply_patch(reference, patch) == new
    return patch


def verify_case(case):
    path = ROOT / case['reference_file']
    assert sha256(path.read_bytes()).hexdigest() == case['reference_sha256']
    reference = read_data(path)
    data = {**reference, 'result': apply_patch(reference['result'], case['certificate_patch'])}
    verification = verify_prepared_search(data)
    assert verification['units']['complete_original_graph_checked'] is True
    return verification


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--raw-result', type=Path)
    parser.add_argument('--bits', type=int, default=80, choices=(64, 80))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    bindings = binding_controls()
    hashes = {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
              for name in SOURCES}
    if args.replay:
        report = read_data(args.replay)
        for name, digest in report['source_sha256'].items():
            assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
        for case in report['cases']:
            verify_case(case)
    else:
        path = ROOT / f'evidence/sparse-oriented-power-{args.bits}.json.xz'
        reference = read_data(path)
        old = reference['result']
        if args.raw_result:
            raw = read_data(args.raw_result)
            assert raw['input'] == reference['input']
            for name, digest in raw.get('search_source_sha256', {}).items():
                assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
            new, cpu = raw['result'], raw['cpu_seconds']
        else:
            started = process_time()
            new = canonical(factor_with_prepared_oriented_refinement(old['integer'], old['cutoff'],
                old['max_queries'], checkpoint_interval=old['checkpoint_interval'],
                max_extraction_columns=old['max_extraction_columns']))
            cpu = process_time() - started
        case = {'bits': args.bits, 'integer': old['integer'],
                'reference_file': str(path.relative_to(ROOT)),
                'reference_sha256': sha256(path.read_bytes()).hexdigest(),
                'certificate_patch': make_patch(old, new), 'cpu_seconds': cpu,
                'same_queries': new['norm_queries'], 'same_factor': new['factor'],
                'same_complete_mathematical_result': True,
                'reference_elapsed_seconds': old['elapsed_seconds'],
                'prepared_elapsed_seconds': new['elapsed_seconds']}
        case['independent_verification'] = verify_case(case)
        report = {'claim': 'Prepared atom batches and small pivots preserve the complete selected relations',
                  'binding_controls': bindings, 'cases': [case], 'source_sha256': hashes,
                  'scope': 'Complete finite matched search and independent replay; no uniform exponent gain.',
                  'timing_scope': 'All search checkpoints charged; partial atom batches flush in the transport phase.'}
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + '\n')
    for name, digest in hashes.items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    print(json.dumps({'binding_controls': bindings,
        'cases': [{key: case[key] for key in ('bits', 'same_queries', 'same_factor',
            'same_complete_mathematical_result', 'prepared_elapsed_seconds', 'cpu_seconds')}
            for case in report['cases']]}, indent=2))


if __name__ == '__main__':
    main()
