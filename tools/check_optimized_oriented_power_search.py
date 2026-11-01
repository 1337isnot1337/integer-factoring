#!/usr/bin/env python3
"""Check matched optimized runs without duplicating their unchanged certificates."""

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import lzma
from pathlib import Path
from time import process_time

from check_sparse_oriented_power_search import SOURCES as REFERENCE_SOURCES, verify_search
from optimized_sparse_oriented_power_search import (
    _optimized, _reference, factor_with_optimized_oriented_refinement,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCES = REFERENCE_SOURCES + ('optimized_sparse_oriented_power_search.py',
    'check_optimized_oriented_power_search.py', 'quotient_update_oriented_basis.py',
    'sparse_unit_presentation.py')


def canonical(value):
    return json.loads(json.dumps(value))


def binding_controls():
    assert _optimized.__code__ is _reference.__code__
    assert _optimized.__globals__ is not _reference.__globals__
    changed = [name for name, value in _reference.__globals__.items()
               if _optimized.__globals__[name] is not value]
    assert sorted(changed) == ['IncrementalOrientedBasis', 'eliminate_unit_residual_columns']
    assert _optimized.__kwdefaults__ == _reference.__kwdefaults__
    return {'identical_search_code_object': True, 'changed_dependency_bindings': sorted(changed),
            'original_global_namespace_unchanged': True}


def make_patch(reference, new):
    """Permit only measured runtime and implementation-counter differences."""
    assert len(reference['checkpoints']) == len(new['checkpoints'])
    patch = {key: new[key] for key in ('elapsed_seconds', 'phase_seconds', 'component_replacements')}
    patch['checkpoint_elapsed_seconds'] = [row['elapsed_seconds'] for row in new['checkpoints']]
    patch['refinement_statistics'] = new['refinement']['statistics']
    patch['tree_update_method'] = new['refinement']['tree_update_method']
    assert apply_patch(reference, patch) == new
    return patch


def apply_patch(reference, patch):
    """Reconstruct a complete run using only the explicit metadata whitelist."""
    result = deepcopy(reference)
    assert set(patch) == {'elapsed_seconds', 'phase_seconds', 'component_replacements',
                          'checkpoint_elapsed_seconds', 'refinement_statistics', 'tree_update_method'}
    for key in ('elapsed_seconds', 'phase_seconds', 'component_replacements'):
        result[key] = patch[key]
    assert len(result['checkpoints']) == len(patch['checkpoint_elapsed_seconds'])
    for row, elapsed in zip(result['checkpoints'], patch['checkpoint_elapsed_seconds']):
        row['elapsed_seconds'] = elapsed
    result['refinement']['statistics'] = patch['refinement_statistics']
    result['refinement']['tree_update_method'] = patch['tree_update_method']
    return result


def reference_data(bits):
    path = ROOT / f'evidence/sparse-oriented-power-{bits}.json.xz'
    return path, json.loads(lzma.decompress(path.read_bytes()))


def verify_case(case):
    path = ROOT / case['reference_file']
    assert sha256(path.read_bytes()).hexdigest() == case['reference_sha256']
    reference = json.loads(lzma.decompress(path.read_bytes()))
    data = {**reference, 'result': apply_patch(reference['result'], case['metadata_patch'])}
    return verify_search(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--raw-directory', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    bindings = binding_controls()
    if args.replay:
        report = json.loads(args.replay.read_text())
        for name, digest in report['source_sha256'].items():
            assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
        for case in report['cases']:
            verify_case(case)
    else:
        cases = []
        for bits in (64, 80):
            path, reference = reference_data(bits)
            old = reference['result']
            if args.raw_directory:
                raw = json.loads((args.raw_directory / f'optimized-oriented-{bits}-result.json').read_text())
                assert raw['input'] == reference['input']
                new, cpu = raw['result'], raw['cpu_seconds']
                assert raw['unit_source_sha256_at_search'] == sha256(
                    (ROOT / 'tools/sparse_unit_presentation.py').read_bytes()).hexdigest()
            else:
                raw = {}
                started = process_time()
                new = canonical(factor_with_optimized_oriented_refinement(
                    old['integer'], old['cutoff'], old['max_queries'],
                    checkpoint_interval=old['checkpoint_interval'],
                    max_extraction_columns=old['max_extraction_columns']))
                cpu = process_time() - started
            patch = make_patch(old, new)
            case = {'bits': bits, 'integer': old['integer'],
                    'reference_file': str(path.relative_to(ROOT)),
                    'reference_sha256': sha256(path.read_bytes()).hexdigest(),
                    'metadata_patch': patch, 'cpu_seconds': cpu,
                    'same_queries': new['norm_queries'], 'same_factor': new['factor'],
                    'entire_result_equal_except_whitelisted_metadata': True,
                    'baseline_elapsed_seconds': raw.get('paired_medians', {}).get(
                        'baseline', {}).get('wall_seconds', old['elapsed_seconds']),
                    'reference_elapsed_seconds': old['elapsed_seconds'],
                    'optimized_elapsed_seconds': new['elapsed_seconds']}
            if 'paired_runs' in raw:
                case['paired_runs'] = raw['paired_runs']
                case['paired_medians'] = raw['paired_medians']
                case['baseline_timing_origin'] = 'Three paired repeats with alternating order'
            else:
                case['baseline_timing_origin'] = 'Retained complete reference run'
            case['independent_verification'] = verify_case(case)
            cases.append(case)
        report = {'claim': 'Exact component replacements preserve complete search results and reduce measured runtime',
                  'binding_controls': bindings, 'cases': cases,
                  'source_sha256': {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                    for name in SOURCES},
                  'scope': 'Matched finite runs with all checkpoint verification charged; no uniform exponent claim.'}
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'binding_controls': bindings,
                     'cases': [{key: case[key] for key in ('bits', 'same_queries', 'same_factor',
                         'baseline_elapsed_seconds', 'optimized_elapsed_seconds', 'cpu_seconds',
                         'entire_result_equal_except_whitelisted_metadata')} for case in report['cases']]}, indent=2))


if __name__ == '__main__':
    main()
