#!/usr/bin/env python3
"""Replay the same complete graph maps with a sparse projected row view."""

import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
from types import FunctionType

import check_signed_residual_graph as reference_checks
from check_sparse_oriented_power_search import canonical, read_data
from signed_residual_graph import eliminate_signed_residual_graph
from sparse_projected_signed_graph import (
    eliminate_sparse_projected_signed_graph, materialize_dense_projection,
    verify_sparse_projected_signed_graph,
)
from unit_presentation_trace import eliminate_unit_residual_columns_compact, verify_unit_presentation_trace

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('sparse_projected_signed_graph.py', 'check_sparse_projected_signed_graph.py',
    'signed_residual_graph.py', 'check_signed_residual_graph.py', 'unit_presentation_trace.py')


def bind(function, changes):
    namespace = dict(function.__globals__)
    namespace.update(changes)
    result = FunctionType(function.__code__, namespace, function.__name__,
                          function.__defaults__, function.__closure__)
    result.__kwdefaults__ = dict(function.__kwdefaults__ or {})
    return result


def small_controls():
    equalities = 0

    def equal(*args, **kwargs):
        nonlocal equalities
        original_rows = canonical(args[0])
        expected = eliminate_signed_residual_graph(*args, **kwargs)
        result = eliminate_sparse_projected_signed_graph(*args, **kwargs)
        assert canonical(args[0]) == original_rows
        assert materialize_dense_projection(result) == expected
        verify_sparse_projected_signed_graph(args[0], canonical(result))
        equalities += 1
        return expected

    report = {}
    for name in ('exhaustive_graph_controls', 'mixed_row_controls', 'ideal_controls'):
        check = bind(getattr(reference_checks, name), {'eliminate_signed_residual_graph': equal})
        result = check()
        report[name] = result if name != 'ideal_controls' else {'factor_controls': len(result)}
    report['rejections'] = bind(reference_checks.rejection_controls, {
        'eliminate_signed_residual_graph': eliminate_sparse_projected_signed_graph,
        'verify_signed_graph_elimination': verify_sparse_projected_signed_graph})()
    report['complete_reference_equalities'] = equalities
    return report


def projection_rejections():
    rows = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
    good = eliminate_sparse_projected_signed_graph(rows, 0)
    candidates = []
    bad = deepcopy(good); bad['sparse_torsion_rows'][0][0] += 1; candidates.append(bad)
    bad = deepcopy(good); bad['sparse_torsion_rows'].append({}); candidates.append(bad)
    bad = deepcopy(good); bad['sparse_torsion_rows'][0][9] = 1; candidates.append(bad)
    bad = deepcopy(good); bad['schema'] = 'wrong'; candidates.append(bad)
    for result in candidates:
        try:
            verify_sparse_projected_signed_graph(rows, result)
        except AssertionError:
            pass
        else:
            raise AssertionError('A corrupted sparse projection was accepted')
    return len(candidates)


def wide_control(size=1024):
    rows = [{1: 1, 2: 1, index+3: 1, **({0: index % 7-3} if index % 7 != 3 else {})}
            for index in range(size)]
    count = size+3
    old = eliminate_signed_residual_graph(rows, 1, generator_count=count, allow_extra_rows=True)
    new = eliminate_sparse_projected_signed_graph(rows, 1, generator_count=count, allow_extra_rows=True)
    assert materialize_dense_projection(new) == old
    trace = eliminate_unit_residual_columns_compact(rows, new, verify_graph=False)
    old_trace = eliminate_unit_residual_columns_compact(rows, old, verify_graph=False)
    assert trace == old_trace
    verifier = bind(verify_unit_presentation_trace,
                    {'verify_signed_graph_elimination': verify_sparse_projected_signed_graph})
    verification = verifier(rows, new, trace, verify_graph=True)
    assert verification['complete_original_graph_checked']
    return {'input_rows': size, 'intermediate_columns': len(new['torsion_columns']),
            'dense_projection_cells': len(old['torsion_rows'])*len(old['torsion_columns']),
            'sparse_projection_entries': sum(map(len, new['sparse_torsion_rows'])),
            'complete_graph_and_unit_outputs_match': True,
            'unit_verification': verification}


def fixture_controls():
    result = []
    for bits in (64, 80):
        path = ROOT / f'evidence/sparse-oriented-power-{bits}.json.xz'
        record = read_data(path)['result']
        rows = [{int(j): value for j, value in row.items()} for row in record['transported_rows']]
        old = record['graph_reduction']
        new = eliminate_sparse_projected_signed_graph(rows, old['initial_size'],
            generator_count=old['old_generator_count'], allow_extra_rows=True)
        assert canonical(materialize_dense_projection(new)) == old
        verify_sparse_projected_signed_graph(rows, canonical(new))
        result.append({'bits': bits, 'source_file': str(path.relative_to(ROOT)),
            'source_sha256': sha256(path.read_bytes()).hexdigest(),
            'rows': len(rows), 'dense_projection_cells': len(old['torsion_rows'])*len(old['torsion_columns']),
            'sparse_projection_entries': sum(map(len, new['sparse_torsion_rows'])),
            'complete_graph_certificate_equal': True})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hashes = {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
              for name in SOURCES}
    report = {'claim': 'The complete integer graph certificate is unchanged by a sparse projected view',
              'small_controls': small_controls(), 'projection_corruptions_rejected': projection_rejections(),
              'wide_control': wide_control(), 'fixtures': fixture_controls(), 'source_sha256': hashes,
              'scope': 'Supplied presentation processing; no relation-supply or factoring-exponent claim.'}
    for name, digest in hashes.items():
        assert sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
