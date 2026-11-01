#!/usr/bin/env python3
"""Exact HNF and strict-factor-gain controls for oriented gcd-free atoms."""

import argparse
from hashlib import sha256
import json
from math import gcd, prod
from pathlib import Path

from arithmetic_two_group import power
from binary_quadratic_forms import ambiguous_factor, inverse_form, reduce_form
from cofactor_collision_relations import partial_leading_word, verify_partial_identity, RationalRowBasis
from oriented_coprime_basis import oriented_coprime_basis, transport_relation_rows
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from relation_two_saturation import first_useful_relation_top

ROOT = Path(__file__).resolve().parents[1]


def root_controls():
    identities = mixed_batches = mixed_splits = 0
    for delta in (-23, -84, -260, -575):
        keys = [(a, b) for a in range(1, 100, 2) if gcd(a, delta) == 1
                for b in range(2 * a) if (b * b - delta) % (4 * a) == 0]
        record = oriented_coprime_basis(delta, keys)
        identities += record['verification']['literal_ideal_identities']
        for norm in range(3, 180, 2):
            if gcd(norm, delta) != 1:
                continue
            same_norm = [(norm, b) for b in range(2 * norm)
                         if (b * b - delta) % (4 * norm) == 0]
            if len(same_norm) >= 4:
                record = oriented_coprime_basis(delta, same_norm)
                assert record['orientation_splits'] > 0
                identities += record['verification']['literal_ideal_identities']
                mixed_batches += 1
                mixed_splits += record['orientation_splits']
    powers = (9, 27, 243, 169, 9 * 169)
    power_keys = [(a, next(b for b in range(1, 2 * a, 2) if (b * b + 23) % (4 * a) == 0))
                  for a in powers]
    power_record = oriented_coprime_basis(-23, power_keys)
    assert [a['norm'] for a in power_record['atoms']] == [3, 169]
    assert [list(map(abs, word)) for word in power_record['words']] == \
        [[2, 0], [3, 0], [5, 0], [0, 1], [2, 1]]
    identities += len(power_keys)
    composite = oriented_coprime_basis(-23, [(39, 17)])
    assert composite['atoms'][0]['norm'] == 39
    mixed = oriented_coprime_basis(-23, [(39, 17), (39, 35)])
    assert mixed['orientation_splits'] == 1
    assert [a['norm'] for a in mixed['atoms']] == [3, 13]
    assert mixed['words'] == [[-1, -1], [-1, 1]]
    identities += 3
    invalid = [lambda: oriented_coprime_basis(23, []),
               lambda: oriented_coprime_basis(-22, []),
               lambda: oriented_coprime_basis(-23, [(2, 1)]),
               lambda: oriented_coprime_basis(-23, [(23, 23)]),
               lambda: oriented_coprime_basis(-23, [(3, 3)]),
               lambda: transport_relation_rows([[1]], 1, [[1]])]
    for call in invalid:
        try:
            call()
        except ValueError:
            pass
        else:
            raise AssertionError('An invalid refinement request was accepted.')
    return {'literal_ideal_identities': identities, 'mixed_orientation_batches': mixed_batches,
            'mixed_orientation_splits': mixed_splits, 'invalid_inputs_rejected': len(invalid),
            'composite_atom_retained': composite, 'mixed_root_example': mixed,
            'repeated_power_example': power_record}


def strict_gain():
    integer, cutoff = 1022117, 7
    delta, identity = -4 * integer, (1, 0, integer)
    base = build_prime_form_factor_base(delta, cutoff)
    records = {r.prime: (r.root, r.positive) for r in base['records']}
    records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
    primes, witnesses, keys = sorted(records), [], []
    for prime, exponent in ((7, 9), (7, 10), (3, 21), (7, 45)):
        target = power(records[prime][1], exponent, identity, compose_forms_hnf)
        assert gcd(target[0], integer) == 1 and target != inverse_form(target)
        split = partial_leading_word(integer, target, primes, records, prod(primes))
        assert split['key'] is not None and verify_partial_identity(integer, target, split, records)
        vector = [0] * len(primes)
        vector[primes.index(prime)] = exponent
        for p, value in split['word'].items():
            vector[primes.index(p)] -= value
        keys.append(split['key'])
        witnesses.append({'prime': prime, 'exponent': exponent, 'target': target,
                          'split': split, 'base_vector': vector})
    assert keys == [(209, 98), (187, 78), (437, 60), (391, 78)]
    assert len({key[0] for key in keys}) == len(keys)
    old_rows = []
    for i, witness in enumerate(witnesses):
        row = witness['base_vector'] + [0] * len(witnesses)
        row[len(primes) + i] = -witness['split']['orientation']
        old_rows.append(row)
    old_generators = [records[p][1] for p in primes] + [w['split']['residual_form'] for w in witnesses]
    inspect = lambda form: ambiguous_factor(form, integer)
    old_top, old_factor, old_stats = first_useful_relation_top(
        old_generators, old_rows, identity, compose_forms_hnf, inverse_form, inspect)
    assert old_top is old_factor is None and old_stats['termination'] == 'odd presented torsion'
    refinement = oriented_coprime_basis(delta, keys)
    transported = transport_relation_rows(old_rows, len(primes), refinement['words'])
    generators = [records[p][1] for p in primes] + [a['form'] for a in refinement['atoms']]
    basis = RationalRowBasis(len(generators))
    for row in transported:
        basis.add(row)
    top, factor, stats = first_useful_relation_top(
        generators, basis.rows, identity, compose_forms_hnf, inverse_form, inspect)
    assert factor == 1009 and top == (1009, 0, 1013) and integer % factor == 0
    cycle = [sum(row[j] for row in transported) for j in range(len(generators))]
    assert cycle == [-2, 22, 64, 0, 0, 0, 0]
    cycle_top, cycle_factor, cycle_stats = first_useful_relation_top(
        [records[p][1] for p in primes], [cycle[:len(primes)]], identity,
        compose_forms_hnf, inverse_form, inspect)
    assert cycle_top == top and cycle_factor == factor
    return {'integer': integer, 'cutoff': cutoff, 'active_primes': primes,
            'prime_records': records, 'witnesses': witnesses, 'old_rows': old_rows,
            'old_stats': old_stats, 'old_factor': old_factor,
            'all_old_residual_norms_distinct': True,
            'no_selected_target_is_an_involution_or_has_a_proper_leading_gcd': True,
            'refinement': refinement, 'transported_rows': transported,
            'independent_rows': basis.rows, 'new_stats': stats,
            'base_cycle_relation': cycle[:len(primes)], 'base_half_relation': [-1, 11, 32],
            'cycle_stats': cycle_stats, 'factor': factor, 'ambiguous_form': top,
            'scope': 'Strict gain for supplied exact witnesses; no uniform search bound.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    controls, gain = root_controls(), strict_gain()
    sources = ('oriented_coprime_basis.py', 'check_oriented_coprime_basis.py',
               'quadratic_ideal_hnf.py', 'cofactor_collision_relations.py',
               'relation_two_saturation.py')
    report = {'claim': 'Exact gcd-free oriented ideal refinement and strict torsion-image gain',
              'controls': controls, 'strict_gain': gain,
              'source_sha256': {'tools/' + name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in sources},
              'scope': 'Supplied-witness refinement; no factoring oracle and no complete exponent gain.'}
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'controls': {k: v for k, v in controls.items() if not isinstance(v, dict)},
                      'factor': gain['factor'], 'ambiguous_form': gain['ambiguous_form'],
                      'residual_keys': gain['refinement']['keys'],
                      'base_relation': gain['base_cycle_relation']}, indent=2))


if __name__ == '__main__':
    main()
