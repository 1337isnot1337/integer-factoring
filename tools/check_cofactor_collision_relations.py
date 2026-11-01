#!/usr/bin/env python3
"""Check exact residual-ideal cancellation and replay retained factoring rows.

The default reconstructs every recorded algebraic witness and reruns the
small searches.  --replay-search reruns every rejected large candidate too.
No finite run is used as a uniform factoring-time theorem.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import gcd, prod
from pathlib import Path

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form, reduce_form, reduced_forms
from cofactor_collision_relations import (
    RationalRowBasis, factor_with_cofactor_collisions, partial_leading_word,
    verify_partial_identity, word_value,
)
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = ROOT / 'evidence/cofactor-collision-relations.json'
SOURCES = ('tools/cofactor_collision_relations.py',
           'tools/check_cofactor_collision_relations.py',
           'tools/leading_power_smooth_relations.py',
           'tools/prime_form_factor_base.py', 'tools/quadratic_ideal_hnf.py',
           'tools/relation_two_saturation.py')
SMALL_CASES = ((15, 2, 10, 25, 100), (21, 2, 15, 25, 100),
               (35, 2, 20, 25, 100), (77, 3, 20, 49, 100),
               (49, 3, 20, 49, 100), (97, 3, 20, 49, 100),
               (539, 3, 30, 49, 100), (959, 2, 30, 25, 100),
               (1022117, 7, 1000, 49, 10000),
               (1000036000099, 31, 10000, 1, 30000),
               (1000036000099, 31, 10000, 31**2, 30000),
               (1000000016000000063, 127, 10000, 1, 30000),
               (1000000016000000063, 127, 10000, 127**2, 30000))


def canonical(value):
    return json.loads(json.dumps(value))


def exact_controls():
    checks = 0
    for integer in range(3, 300, 2):
        delta = -integer if integer % 4 == 3 else -4 * integer
        for cutoff in (2, 5, 13):
            records = {r.prime: (r.root, r.positive)
                       for r in build_prime_form_factor_base(delta, cutoff)['records']}
            if integer % 4 == 1:
                records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
            primes = sorted(records)
            for form in reduced_forms(delta):
                if gcd(form[0], integer) != 1:
                    continue
                split = partial_leading_word(integer, form, primes, records)
                assert verify_partial_identity(integer, form, split, records)
                assert prod(p**e for p, e in split['factors']) == split['smooth_part']
                assert split['smooth_part'] * split['residual'] == form[0]
                checks += 1
    fixtures = []
    for integer, forms in ((959, ((15, 1, 16), (15, -1, 16), (15, 11, 18))),
                            (4969, ((70, 2, 71), (70, -2, 71)))):
        delta = -integer if integer % 4 == 3 else -4 * integer
        records = {r.prime: (r.root, r.positive)
                   for r in build_prime_form_factor_base(delta, 2)['records']}
        if integer % 4 == 1:
            records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
        splits = [partial_leading_word(integer, form, sorted(records), records)
                  for form in forms]
        assert all(verify_partial_identity(integer, form, split, records)
                   for form, split in zip(forms, splits))
        assert splits[0]['key'] == splits[1]['key']
        assert splits[0]['orientation'] == -splits[1]['orientation']
        if integer == 959:
            assert splits[0]['residual'] == splits[2]['residual'] == 15
            assert splits[0]['key'] != splits[2]['key']
            assert splits[0]['residual_form'] != splits[2]['residual_form']
        else:
            assert splits[0]['smooth_part'] == 2 and splits[0]['residual'] == 35
        fixtures.append({'integer': integer, 'forms': forms, 'splits': splits})
    return {'exhaustive_partial_ideal_identities': checks,
            'composite_residual_controls': fixtures,
            'wrong_norm_only_key_rejected': True,
            'ramified_two_and_conjugate_orientations_checked': True}


def verify_record(record):
    integer, delta = record['integer'], record['discriminant']
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    records = {int(p): (value[0], tuple(value[1]))
               for p, value in record['prime_records'].items()}
    primes = record['active_primes']
    indices = {p: i for i, p in enumerate(primes)}
    expected = {r.prime: (r.root, r.positive)
                for r in build_prime_form_factor_base(delta, record['cutoff'])['records']}
    if integer % 4 == 1:
        expected[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    assert records == expected and primes == sorted(expected)
    basis = RationalRowBasis(len(primes))
    for witness, row in zip(record['witnesses'], record['relation_matrix']):
        vectors = []
        for item in (witness['current'], witness['previous']):
            if item is None:
                continue
            target = tuple(item['target'])
            assert target == power(records[item['prime']][1], item['exponent'],
                                   identity, compose_forms_hnf)
            split = partial_leading_word(integer, target, primes, records)
            assert verify_partial_identity(integer, target, split, records)
            for key in ('smooth_part', 'residual', 'factors', 'word', 'key',
                        'orientation', 'residual_form', 'trial_divisions'):
                assert canonical(split[key]) == item[key], key
            assert split['residual'] <= record['max_cofactor']
            vector = [0] * len(primes)
            vector[indices[item['prime']]] = item['exponent']
            for p, exponent in split['word'].items():
                vector[indices[p]] -= exponent
            assert [[j, value] for j, value in enumerate(vector) if value] == item['vector']
            if split['key'] is not None:
                vector = [split['orientation'] * value for value in vector]
            vectors.append(vector)
        if witness['previous'] is None:
            assert witness['current']['residual'] == 1
            reconstructed = vectors[0]
        else:
            assert witness['current']['key'] == witness['previous']['key']
            reconstructed = [a - b for a, b in zip(*vectors)]
        assert reconstructed == row and basis.add(row)
        word = {p: exponent for p, exponent in zip(primes, row) if exponent}
        assert word_value(word, records, identity) == identity
    assert len(record['witnesses']) == len(record['relation_matrix']) == record['relation_rank']
    assert record['independent_collisions'] + record['independent_smooth_rows'] == record['relation_rank']
    assert record['factoring_oracle_calls'] == record['residual_primality_tests'] == 0
    assert record['new_primes_introduced'] == 0
    if record['factor'] is not None:
        assert 1 < record['factor'] < integer and integer % record['factor'] == 0
    return len(record['witnesses'])


def run_case(case):
    integer, cutoff, window, cap, query_cap = case
    return canonical(factor_with_cofactor_collisions(
        integer, cutoff, window, max_cofactor=cap, max_queries=query_cap))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument('--import-large', type=Path, action='append', default=[])
    parser.add_argument('--regenerate', action='store_true')
    parser.add_argument('--replay-search', action='store_true')
    args = parser.parse_args()
    controls = canonical(exact_controls())
    small = [run_case(case) for case in SMALL_CASES]
    source_hashes = {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES}
    if args.regenerate:
        large = [json.loads(path.read_text()) for path in args.import_large]
        evidence = {'status': 'Exact finite selector; uniform admissible-query bound remains open',
                    'source_sha256': source_hashes, 'controls': controls,
                    'small_cases': small, 'large_cases': large,
                    'default_replay_scope': 'Reconstruct all retained ideal and relation witnesses; rerun all small candidates.',
                    'full_replay_flag': '--replay-search', 'factoring_exponent_improved': False}
        args.evidence.write_text(json.dumps(evidence, separators=(',', ':')) + '\n')
    else:
        evidence = json.loads(args.evidence.read_text())
        assert evidence['source_sha256'] == source_hashes
        assert evidence['controls'] == controls and evidence['small_cases'] == small
    relation_checks = sum(verify_record(record) for record in evidence['small_cases'] + evidence['large_cases'])
    full_searches = 0
    if args.replay_search:
        for record in evidence['large_cases']:
            case = (record['integer'], record['cutoff'], record['window'],
                    record['max_cofactor'], record['max_queries'])
            assert run_case(case) == record
            full_searches += 1
    print(json.dumps({'status': 'pass', 'partial_ideal_controls': controls['exhaustive_partial_ideal_identities'],
                      'independent_relation_witnesses': relation_checks,
                      'small_full_searches': len(small), 'large_certificate_replays': len(evidence['large_cases']),
                      'large_full_searches': full_searches,
                      'large_summaries': [{key: record[key] for key in
                                          ('input_bits', 'max_cofactor', 'norm_queries', 'factor',
                                           'relation_rank', 'independent_collisions', 'retained_cofactor_keys')}
                                         for record in evidence['large_cases']]}, indent=2))


if __name__ == '__main__':
    main()
