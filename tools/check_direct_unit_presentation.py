#!/usr/bin/env python3
"""Check direct unit traces against full row replay and independent characters."""

import argparse
from copy import deepcopy
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random
from time import process_time

from binary_quadratic_forms import ambiguous_factor, inverse_form
from check_sparse_oriented_power_search import canonical, read_data
from cofactor_collision_relations import RationalRowBasis
from direct_unit_presentation import (
    _finish, _identity_presentation, _normalize_rows, _producer, _stream_fingerprint,
    eliminate_direct_unit_presentation, verify_direct_unit_presentation,
)
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from unit_presentation_trace import (
    _fingerprint as reference_fingerprint, _finish as reference_finish,
    eliminate_unit_residual_columns_compact as reference_producer,
    verify_unit_presentation_trace,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('direct_unit_presentation.py', 'check_direct_unit_presentation.py',
    'unit_presentation_trace.py', 'compact_oriented_power_search.py',
    'check_sparse_oriented_power_search.py', 'cofactor_collision_relations.py',
    'relation_two_saturation.py', 'quadratic_ideal_hnf.py', 'binary_quadratic_forms.py')


def source_hashes():
    return {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
            for name in SOURCES}


def binding_controls():
    assert _producer.__code__ is reference_producer.__code__
    assert _finish.__code__ is reference_finish.__code__
    assert [name for name, value in reference_producer.__globals__.items()
            if _producer.__globals__[name] is not value] == ['_finish']
    assert [name for name, value in reference_finish.__globals__.items()
            if _finish.__globals__[name] is not value] == ['_fingerprint']
    return {'identical_unit_producer_code': True,
            'only_fingerprint_serialization_replaced': True,
            'original_namespaces_unchanged': True}


def character_count(rows, width, modulus):
    return sum(all(sum(value * vector[int(j)] for j, value in row.items()) % modulus == 0
                   for row in rows) for vector in product(range(modulus), repeat=width))


def bounded_controls():
    rng = Random(2026110117)
    families = [([], 0, 0), ([{}], 0, 0), ([{}], 1, 4),
                ([{1: 2}], 1, 3), ([{0: 7, 1: -1}], 1, 2),
                ([{1: 2}, {1: 3}], 1, 2),
                ([{0: 2, 1: 1}, {0: -3, 1: -1, 2: 2}], 1, 4)]
    for _ in range(600):
        width = rng.randrange(0, 9)
        base = rng.randrange(width + 1)
        rows = []
        for _ in range(rng.randrange(0, 12)):
            row = {j: rng.randrange(-3, 4) for j in range(width) if rng.randrange(3) == 0}
            rows.append({j: value for j, value in row.items() if value})
        families.append((rows, base, width))
    pivots = substitutions = characters = comparisons = roundtrips = 0
    saved = None
    for rows, base, width in families:
        before = deepcopy(rows)
        identity = _identity_presentation(rows, base, width)
        assert _stream_fingerprint(identity) == reference_fingerprint(identity)
        expected = reference_producer(rows, identity, verify=False, verify_graph=False)
        result = eliminate_direct_unit_presentation(rows, base, width)
        stripped = {k: v for k, v in result.items() if k not in ('input_kind', 'verification')}
        assert stripped == expected
        full_scan = verify_unit_presentation_trace(rows, identity, result, verify_graph=False)
        assert result['verification']['row_substitutions'] == full_scan['row_substitutions']
        assert rows == before
        loaded = canonical(result)
        verify_direct_unit_presentation(canonical(rows), base, width, loaded)
        assert eliminate_direct_unit_presentation(canonical(rows), base, width) == result
        roundtrips += 1
        pivots += len(result['pivots'])
        substitutions += full_scan['row_substitutions']
        comparisons += 1
        if width <= 4:
            for modulus in (2, 3, 4):
                original = character_count(rows, width, modulus)
                presented = character_count(result['output_rows'],
                                            len(result['new_generator_indices']), modulus)
                assert original == presented
                projected = character_count(
                    [{j: value for j, value in enumerate(row) if value}
                     for row in result['torsion_rows']],
                    len(result['torsion_generator_indices']), modulus)
                assert original == projected * modulus ** len(result['free_root_columns'])
                characters += 2
        if saved is None and len(result['pivots']) >= 2 and result['torsion_rows']:
            saved = (rows, base, width, result)
    assert saved is not None
    rows, base, width, good = saved
    rejected = 0
    mutations = (
        lambda r: r.__setitem__('input_kind', 'wrong'),
        lambda r: r.__setitem__('input_presentation_sha256', '0' * 64),
        lambda r: r['pivots'][0].__setitem__(2, -r['pivots'][0][2]),
        lambda r: r['pivots'].append(r['pivots'][0]),
        lambda r: r['pivots'].pop(),
        lambda r: r['new_generator_indices'].reverse(),
        lambda r: r['free_root_columns'].append(0),
        lambda r: r['torsion_rows'][0].__setitem__(0, r['torsion_rows'][0][0] + 1),
    )
    for mutation in mutations:
        bad = deepcopy(good)
        mutation(bad)
        try:
            verify_direct_unit_presentation(rows, base, width, bad)
        except (AssertionError, ValueError, IndexError, KeyError, TypeError):
            rejected += 1
        else:
            raise AssertionError('A corrupted trace was accepted.')
    invalid = [([{0: 1}], 0, 0), ([{1: True}], 0, 2), ([{-1: 2}], 0, 2),
               ([{'0': 1, '00': 2}], 0, 2), ([{0: 1}], 2, 1), ([[]], 0, 1)]
    invalid_rejected = 0
    for rows, base, width in invalid:
        try:
            eliminate_direct_unit_presentation(rows, base, width)
        except (ValueError, TypeError):
            invalid_rejected += 1
        else:
            raise AssertionError('An invalid input was accepted.')
    return {'complete_reference_equalities': comparisons,
            'full_scan_verifier_equalities': comparisons,
            'unit_isomorphisms': pivots, 'row_substitutions': substitutions,
            'independent_character_counts': characters,
            'json_roundtrips': roundtrips, 'corrupted_trace_rejections': rejected,
            'invalid_input_rejections': invalid_rejected,
            'coefficient_two_and_coprime_nonunit_controls': True}


def fixture_control(bits, *, retained_trace=None):
    path = ROOT / f'evidence/sparse-oriented-power-{bits}.json.xz'
    reference = read_data(path)
    record = reference['result']
    rows = record['transported_rows']
    base = len(record['active_primes'])
    atoms = record['refinement']['atoms']
    width = base + len(atoms)
    tick = process_time()
    trace = (eliminate_direct_unit_presentation(rows, base, width)
             if retained_trace is None else retained_trace)
    cpu = process_time() - tick
    replay = verify_direct_unit_presentation(rows, base, width, trace)
    matrix = trace['torsion_rows']
    selector = RationalRowBasis(len(trace['torsion_generator_indices']))
    for row in matrix:
        selector.add(row)
    records = {int(p): value for p, value in record['prime_records'].items()}
    generators = [tuple(records[p][1]) for p in record['active_primes']]
    generators += [tuple(a['form']) for a in atoms]
    selected = [generators[j] for j in trace['torsion_generator_indices']]
    delta = record['discriminant']
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    top, factor, stats = first_useful_relation_top(selected, selector.rows,
        identity, compose_forms_hnf, inverse_form,
        lambda form: ambiguous_factor(form, record['integer']), verify_relations=True)
    assert factor in reference['input']['input_primes']
    return {'bits': bits, 'reference_file': str(path.relative_to(ROOT)),
            'reference_sha256': sha256(path.read_bytes()).hexdigest(),
            'trace': trace, 'producer_with_verification_cpu_seconds': cpu,
            'independent_replay': replay, 'input_rows': len(rows),
            'input_generators': width, 'final_rows': len(matrix),
            'final_columns': len(selected), 'independent_rank': len(selector.rows),
            'factor': factor, 'ambiguous_form': top, 'extraction': stats}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hashes = source_hashes()
    bindings = binding_controls()
    controls = bounded_controls()
    if args.replay:
        report = read_data(args.replay)
        assert report['source_sha256'] == hashes and report['bounded_controls'] == controls
        for case in report['cases']:
            rebuilt = fixture_control(case['bits'], retained_trace=case['trace'])
            for name in rebuilt:
                if name != 'producer_with_verification_cpu_seconds':
                    assert canonical(rebuilt[name]) == case[name], name
    else:
        cases = [fixture_control(bits) for bits in (64, 80)]
        report = {'claim': 'Direct unit traces preserve the full supplied integer presentation',
                  'binding_controls': bindings, 'bounded_controls': controls,
                  'cases': cases, 'source_sha256': hashes,
                  'scope': 'Exact supplied-matrix reduction; no relation-supply or uniform exponent claim.',
                  'timing_scope': 'The component includes its incidence verifier; input loading and final extraction are separate.'}
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + '\n')
    assert source_hashes() == hashes
    print(json.dumps({'binding_controls': bindings, 'bounded_controls': controls,
        'cases': [{k: v for k, v in case.items() if k != 'trace'} for case in report['cases']]}, indent=2))


if __name__ == '__main__':
    main()
