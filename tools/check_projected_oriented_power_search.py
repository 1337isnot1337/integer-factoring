#!/usr/bin/env python3
"""Check complete search equality after expanding the sparse intermediate projection."""

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from time import process_time

from check_prepared_oriented_power_search import (
    SOURCES as PREPARED_SOURCES, apply_patch as apply_prepared_patch,
    make_patch as make_prepared_patch,
)
from sparse_projected_signed_graph import materialize_dense_projection
from check_sparse_oriented_power_search import canonical, read_data
from projected_oriented_power_search import (
    _projected, _reference, factor_with_projected_oriented_refinement, verify_projected_search,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCES = PREPARED_SOURCES + ('sparse_projected_signed_graph.py',
    'projected_oriented_power_search.py', 'check_projected_oriented_power_search.py')


def binding_controls():
    assert _projected.__code__ is _reference.__code__
    assert _projected.__globals__ is not _reference.__globals__
    changed = sorted(name for name, value in _reference.__globals__.items()
                     if _projected.__globals__[name] is not value)
    assert changed == ['IncrementalOrientedBasis', 'RationalRowBasis',
        'eliminate_signed_residual_graph', 'eliminate_unit_residual_columns']
    assert _projected.__kwdefaults__ == _reference.__kwdefaults__
    return {'identical_search_code_object': True,
            'changed_dependency_bindings': changed,
            'original_global_namespace_unchanged': True}


def apply_patch(reference, patch):
    assert set(patch) == {'prepared_patch', 'sparse_torsion_rows'}
    result = apply_prepared_patch(reference, patch['prepared_patch'])
    graph = dict(result['graph_reduction'])
    del graph['torsion_rows']
    graph['schema'] = 'sparse-projected-signed-graph-v1'
    graph['sparse_torsion_rows'] = patch['sparse_torsion_rows']
    assert canonical(materialize_dense_projection(graph)) == result['graph_reduction']
    result['graph_reduction'] = graph
    return result


def make_patch(reference, new):
    expanded = deepcopy(new)
    expanded['graph_reduction'] = canonical(materialize_dense_projection(new['graph_reduction']))
    patch = {'prepared_patch': make_prepared_patch(reference, expanded),
             'sparse_torsion_rows': new['graph_reduction']['sparse_torsion_rows']}
    assert apply_patch(reference, patch) == new
    return patch


def verify_case(case):
    path = ROOT / case['reference_file']
    assert sha256(path.read_bytes()).hexdigest() == case['reference_sha256']
    reference = read_data(path)
    data = {**reference, 'result': apply_patch(reference['result'], case['certificate_patch'])}
    verification = verify_projected_search(data)
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
            new = canonical(factor_with_projected_oriented_refinement(old['integer'], old['cutoff'],
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
                'projected_elapsed_seconds': new['elapsed_seconds']}
        case['independent_verification'] = verify_case(case)
        report = {'claim': 'Sparse intermediate projection preserves the complete selected relations',
                  'binding_controls': bindings, 'cases': [case], 'source_sha256': hashes,
                  'scope': 'Complete finite matched search and independent replay; no uniform exponent gain.',
                  'timing_scope': 'All search checkpoints charged; partial atom batches flush in the transport phase.'}
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + '\n')
    for name, digest in hashes.items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    print(json.dumps({'binding_controls': bindings,
        'cases': [{key: case[key] for key in ('bits', 'same_queries', 'same_factor',
            'same_complete_mathematical_result', 'projected_elapsed_seconds', 'cpu_seconds')}
            for case in report['cases']]}, indent=2))


if __name__ == '__main__':
    main()
