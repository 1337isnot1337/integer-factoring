#!/usr/bin/env python3
"""Independently replay matched power prefixes and emitted integer rows."""

import argparse
from hashlib import sha256
import json
from math import gcd, isqrt, prod
from pathlib import Path

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form
from cofactor_collision_relations import RationalRowBasis, partial_leading_word, factor_with_cofactor_collisions
from oriented_coprime_basis import verify_oriented_basis, transport_relation_rows
from oriented_cofactor_prefixes import compare_power_prefix, eliminate_atom_columns
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top
from weighted_norm_relations import least_power_above

ROOT = Path(__file__).resolve().parents[1]


def rank(rows, width):
    basis = RationalRowBasis(width)
    for row in rows:
        basis.add(row)
    return len(basis.rows)


def verify_emission(rows, base_size, emitted):
    width = len(rows[0])
    assert len(emitted['rows']) == len(emitted['combinations'])
    for output, combination in zip(emitted['rows'], emitted['combinations']):
        assert all(isinstance(c, int) and 0 <= i < len(rows) for i, c in combination)
        literal = [sum(c * rows[i][j] for i, c in combination) for j in range(width)]
        assert literal == output + [0] * (width - base_size)
    expected = rank(rows, width) - rank([row[base_size:] for row in rows], width - base_size)
    assert rank(emitted['rows'], base_size) == len(emitted['rows']) == expected
    return expected


def verify_prefix(record):
    n = record['integer']
    delta = -n if n % 4 == 3 else -4 * n
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    height = isqrt((-delta) // 3)
    primes = record['active_primes']
    records = {int(p): (value[0], tuple(value[1])) for p, value in record['prime_records'].items()}
    keys, old_parts = [], []
    assert len(record['witnesses']) == record['norm_queries']
    for query, witness in enumerate(record['witnesses']):
        offset, column = divmod(query, len(primes))
        prime = primes[column]
        exponent = least_power_above(prime, height) + offset
        assert witness['prime'] == prime and witness['exponent'] == exponent
        target = power(records[prime][1], exponent, identity, compose_forms_hnf)
        assert tuple(witness['target']) == target and gcd(target[0], n) == 1
        if target != identity and target == inverse_form(target):
            assert ambiguous_factor(target, n) is None
        split = partial_leading_word(n, target, primes, records, prod(primes))
        stored = witness['split']
        assert stored['residual'] == split['residual']
        assert stored['orientation'] == split['orientation']
        assert {int(p): e for p, e in stored['word'].items()} == split['word']
        assert (tuple(stored['key']) if stored['key'] is not None else None) == split['key']
        if split['key'] is not None and split['key'] not in keys:
            keys.append(split['key'])
        row = [0] * len(primes)
        row[column] = exponent
        for p, count in split['word'].items():
            row[primes.index(p)] -= count
        old_parts.append((row, split['key'], split['orientation']))
    assert keys == list(map(tuple, record['keys']))
    rows = []
    for initial, key, sign in old_parts:
        row = initial + [0] * len(keys)
        if key is not None:
            row[len(primes) + keys.index(key)] -= sign
        rows.append(row)
    assert rows == record['old_rows']
    refined = record['refinement']
    verify_oriented_basis(delta, keys, refined)
    transported = transport_relation_rows(rows, len(primes), refined['words'])
    assert transported == record['transported_rows']
    old_rank = verify_emission(rows, len(primes), record['old_emitted'])
    new_rank = verify_emission(transported, len(primes), record['new_emitted'])
    checkpoint = record['checkpoints'][-1]
    assert old_rank == checkpoint['old_base_rank'] and new_rank == checkpoint['new_base_rank']
    assert checkpoint['new_factor'] and n % checkpoint['new_factor'] == 0
    for name in ('old', 'new'):
        _, factor, _ = first_useful_relation_top(
            [records[p][1] for p in primes], record[name + '_emitted']['rows'], identity,
            compose_forms_hnf, inverse_form, lambda form: ambiguous_factor(form, n))
        assert factor == checkpoint[name + '_factor']
    legacy = factor_with_cofactor_collisions(n, record['cutoff'], record['norm_queries'],
        max_queries=record['norm_queries'], verify_candidates=True)
    if checkpoint['old_factor'] is None:
        assert legacy['factor'] is None and legacy['norm_queries'] == record['norm_queries']
        assert legacy['relation_rank'] == old_rank
    else:
        assert legacy['factor'] and n % legacy['factor'] == 0
    return {'literal_power_checks': len(record['witnesses']),
            'old_integer_combinations': old_rank, 'new_integer_combinations': new_rank,
            'legacy_selector': {key: legacy[key] for key in ('factor', 'norm_queries',
                'relation_rank', 'retained_cofactor_keys', 'candidate_relations', 'dependent_relations')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    trap = eliminate_atom_columns([[2, 2], [0, 1]], 1)
    assert trap['rows'] == [[-2]] and trap['combinations'] == [[(0, -1), (1, 2)]]
    cases = []
    for n, y, limit, batch in ((1022117, 7, 128, 8), (1040399, 5, 128, 8),
                               (100160063, 7, 128, 8), (1000036000099, 31, 512, 64)):
        record = compare_power_prefix(n, y, limit, batch_size=batch)
        # Replay the transport representation after a JSON round trip.
        verification = verify_prefix(json.loads(json.dumps(record)))
        cases.append({'result': record, 'verification': verification})
    expected = ((16, 1, 2, 1019), (64, 2, 3, 10007), (320, 0, 3, 1000003))
    for case, (queries, old_rank, new_rank, factor) in zip(cases[1:], expected):
        record = case['result']
        final = record['checkpoints'][-1]
        assert (record['norm_queries'], final['old_base_rank'], final['new_base_rank'],
                final['new_factor']) == (queries, old_rank, new_rank, factor)
        assert final['old_factor'] is None
    sources = ('oriented_cofactor_prefixes.py', 'check_oriented_cofactor_prefixes.py',
               'oriented_coprime_basis.py', 'cofactor_collision_relations.py',
               'leading_power_smooth_relations.py', 'quadratic_ideal_hnf.py',
               'relation_two_saturation.py')
    report = {'claim': 'Exact oriented refinement improves emitted rows on matched selector prefixes',
              'integer_combination_denominator_control': trap, 'cases': cases,
              'source_sha256': {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in sources},
              'scope': 'Finite matched-prefix row and factor gains; no wall-time or uniform exponent claim.'}
    if args.output:
        args.output.write_text(json.dumps(report, sort_keys=True, separators=(',', ':')) + '\n')
    print(json.dumps({'literal_power_checks': sum(c['verification']['literal_power_checks'] for c in cases),
                     'cases': [{'integer': c['result']['integer'],
                                'queries': c['result']['norm_queries'],
                                **{k: v for k, v in c['result']['checkpoints'][-1].items()
                                   if not k.endswith('_stats')}} for c in cases]}, indent=2))


if __name__ == '__main__':
    main()
