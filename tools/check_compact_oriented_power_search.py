#!/usr/bin/env python3
"""Replay compact search certificates against retained complete reference runs."""

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from time import process_time

from check_sparse_oriented_power_search import SOURCES as REFERENCE_SOURCES, canonical, read_data
from compact_oriented_power_search import (
    _compact, _reference, factor_with_compact_oriented_refinement, verify_compact_search,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCES = REFERENCE_SOURCES + ('compact_oriented_power_search.py',
    'check_compact_oriented_power_search.py', 'quotient_update_oriented_basis.py',
    'primitive_integer_row_basis.py', 'unit_presentation_trace.py')
FINAL_PRESENTATION_FIELDS = ('counts', 'free_root_columns', 'initial_size',
    'new_generator_indices', 'output_rows', 'torsion_columns',
    'torsion_generator_indices', 'torsion_rows')
PATCH_FIELDS = {'elapsed_seconds', 'phase_seconds', 'component_replacements',
    'checkpoint_elapsed_seconds', 'refinement_statistics', 'tree_update_method',
    'unit_reduction'}


def binding_controls():
    assert _compact.__code__ is _reference.__code__
    assert _compact.__globals__ is not _reference.__globals__
    changed = sorted(name for name, value in _reference.__globals__.items()
                     if _compact.__globals__[name] is not value)
    assert changed == ['IncrementalOrientedBasis', 'RationalRowBasis',
                       'eliminate_unit_residual_columns']
    assert _compact.__kwdefaults__ == _reference.__kwdefaults__
    return {'identical_search_code_object': True, 'changed_dependency_bindings': changed,
            'original_global_namespace_unchanged': True}


def apply_patch(reference, patch):
    """Replace an explicit metadata whitelist and one independently replayed proof."""
    assert set(patch) == PATCH_FIELDS
    result = deepcopy(reference)
    for key in ('elapsed_seconds', 'phase_seconds', 'component_replacements'):
        result[key] = patch[key]
    assert len(result['checkpoints']) == len(patch['checkpoint_elapsed_seconds'])
    for row, elapsed in zip(result['checkpoints'], patch['checkpoint_elapsed_seconds']):
        row['elapsed_seconds'] = elapsed
    result['refinement']['statistics'] = patch['refinement_statistics']
    result['refinement']['tree_update_method'] = patch['tree_update_method']
    trace = patch['unit_reduction']
    assert trace['schema'] == 'integer-unit-substitution-trace-v1'
    for key in FINAL_PRESENTATION_FIELDS:
        assert reference['unit_reduction'][key] == trace[key], key
    result['unit_reduction'] = trace
    return result


def make_patch(reference, new):
    patch = {key: new[key] for key in ('elapsed_seconds', 'phase_seconds',
                                      'component_replacements', 'unit_reduction')}
    patch['checkpoint_elapsed_seconds'] = [row['elapsed_seconds'] for row in new['checkpoints']]
    patch['refinement_statistics'] = new['refinement']['statistics']
    patch['tree_update_method'] = new['refinement']['tree_update_method']
    assert apply_patch(reference, patch) == new
    return patch


def verify_case(case):
    path = ROOT / case['reference_file']
    assert sha256(path.read_bytes()).hexdigest() == case['reference_sha256']
    reference = read_data(path)
    data = {**reference, 'result': apply_patch(reference['result'], case['certificate_patch'])}
    verification = verify_compact_search(data)
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
            new = canonical(factor_with_compact_oriented_refinement(
                old['integer'], old['cutoff'], old['max_queries'],
                checkpoint_interval=old['checkpoint_interval'],
                max_extraction_columns=old['max_extraction_columns']))
            cpu = process_time() - started
        patch = make_patch(old, new)
        case = {'bits': args.bits, 'integer': old['integer'],
                'reference_file': str(path.relative_to(ROOT)),
                'reference_sha256': sha256(path.read_bytes()).hexdigest(),
                'certificate_patch': patch, 'cpu_seconds': cpu,
                'same_queries': new['norm_queries'], 'same_factor': new['factor'],
                'same_complete_mathematical_result': True,
                'same_final_presentation_fields': list(FINAL_PRESENTATION_FIELDS),
                'reference_elapsed_seconds': old['elapsed_seconds'],
                'compact_elapsed_seconds': new['elapsed_seconds'],
                'reference_unit_certificate_json_bytes': len(json.dumps(
                    old['unit_reduction'], sort_keys=True, separators=(',', ':')).encode()),
                'compact_unit_certificate_json_bytes': len(json.dumps(
                    new['unit_reduction'], sort_keys=True, separators=(',', ':')).encode())}
        case['independent_verification'] = verify_case(case)
        report = {'claim': 'Exact component replacements retain all witnesses and final relations while compacting the proof',
                  'binding_controls': bindings, 'cases': [case], 'source_sha256': hashes,
                  'scope': 'Finite complete search and independent certificate replay; no uniform exponent gain.',
                  'timing_scope': 'Single complete run; other research processes were active. Wall and CPU times are separate.'}
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + '\n')
    for name, digest in hashes.items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    print(json.dumps({'binding_controls': bindings,
        'cases': [{key: case[key] for key in ('bits', 'same_queries', 'same_factor',
            'same_complete_mathematical_result', 'compact_elapsed_seconds', 'cpu_seconds',
            'reference_unit_certificate_json_bytes', 'compact_unit_certificate_json_bytes')}
            for case in report['cases']]}, indent=2))


if __name__ == '__main__':
    main()
