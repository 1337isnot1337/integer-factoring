"""Replay the complete fixed-base power scan and verify its exact witnesses."""

import argparse
from hashlib import sha256
import json
from math import gcd, prod
from pathlib import Path

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form, reduce_form, reduced_forms
from leading_power_smooth_relations import (compose_with_prime_form,
                                           factor_with_smooth_leading_powers)
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from weighted_norm_relations import least_power_above


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'evidence/leading-power-smooth-relations.json'
SOURCES = ('tools/leading_power_smooth_relations.py',
           'tools/check_leading_power_smooth_relations.py')
CASES = (
    (15, 2, 20), (35, 5, 20), (77, 3, 20), (539, 3, 30),
    (49, 3, 20), (21, 2, 20), (97, 3, 20), (1022117, 7, 1000),
    (1000036000099, 31, 5000), (1000000016000000063, 127, 10000),
    (4951760154835678088235319297, 1021, 32),
    (4951760154835678088235319297, 1021, 10000),
    (21267647932558653899591465697288388633, 4093, 100000),
)


def canonical(value):
    return json.loads(json.dumps(value))


def finite_controls():
    composition_cases = shared_prime_cases = smoothness_cases = 0
    for integer in range(3, 152, 2):
        delta = -integer if integer % 4 == 3 else -4 * integer
        records = {r.prime: (r.root, r.positive)
                   for r in build_prime_form_factor_base(delta, 19)['records']}
        if integer % 4 == 1:
            records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
        for form in reduced_forms(delta):
            for prime, (root, prime_form) in records.items():
                assert compose_with_prime_form(form, prime, root, prime_form) == \
                    compose_forms_hnf(form, prime_form)
                composition_cases += 1
                shared_prime_cases += form[0] % prime == 0
    for primes in ((), (2,), (2, 3, 5), (3, 7, 13), (2, 5, 11)):
        for height in (7, 31, 127, 500):
            product = prod(primes)
            for n in range(1, height + 1):
                remaining = n
                for prime in primes:
                    while remaining % prime == 0:
                        remaining //= prime
                assert (pow(product, height.bit_length(), n) == 0) == (remaining == 1)
                smoothness_cases += 1
    return {'prime_compositions': composition_cases,
            'shared_prime_fallbacks': shared_prime_cases,
            'modular_smoothness_controls': smoothness_cases}


def verify_witnesses(record):
    delta = record['discriminant']
    identity = (1, delta % 2, (delta % 2 - delta) // 4)
    records = {int(p): value for p, value in record['prime_records'].items()}
    primes = record['active_primes']
    indices = {p: i for i, p in enumerate(primes)}
    target_primes = set()
    for witness, row in zip(record['witnesses'], record['relation_matrix']):
        prime, exponent = witness['prime'], witness['exponent']
        target = tuple(witness['target'])
        assert prime not in target_primes
        target_primes.add(prime)
        assert exponent >= least_power_above(prime, record['height_bound'])
        assert target == power(tuple(records[prime][1]), exponent, identity,
                               compose_forms_hnf)
        assert 1 <= target[0] <= record['height_bound']
        assert gcd(target[0], record['integer']) == 1
        assert all(p in records and e > 0 for p, e in witness['factors'])
        assert len({p for p, _ in witness['factors']}) == len(witness['factors'])
        assert prod(p ** e for p, e in witness['factors']) == target[0]
        value = identity
        expected = [0] * len(primes)
        expected[indices[prime]] = exponent
        for q, e in witness['word'].items():
            q = int(q)
            generator = tuple(records[q][1])
            expected[indices[q]] -= e
            if e < 0:
                generator, e = inverse_form(generator), -e
            value = compose_forms_hnf(value, power(
                generator, e, identity, compose_forms_hnf))
        assert value == target and expected == row
        assert prod(int(q) ** abs(e) for q, e in witness['word'].items()) == target[0]
    assert len(target_primes) == record['relation_rank'] == len(record['relation_matrix'])
    assert record['norm_queries'] == sum(item['queries'] for item in record['attempts'])
    assert record['factoring_oracle_calls'] == record['new_primes_introduced'] == 0
    if record['factor'] is not None:
        assert 1 < record['factor'] < record['integer']
        assert record['integer'] % record['factor'] == 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regenerate', action='store_true')
    args = parser.parse_args()
    controls = finite_controls()
    hashes = {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES}
    computed = [canonical(factor_with_smooth_leading_powers(*case)) for case in CASES]
    reference = canonical(factor_with_smooth_leading_powers(
        4951760154835678088235319297, 1021, 10000, use_primorial=False))
    accelerated = next(r for r in computed
                       if r['integer'] == 4951760154835678088235319297 and r['window'] == 10000)
    for field in ('witnesses', 'attempts', 'relation_matrix', 'factor', 'norm_queries'):
        assert reference[field] == accelerated[field]
    reference_summary = {field: reference[field]
                         for field in ('norm_queries', 'trial_divisions', 'factor')}
    evidence = {'status': 'exact finite searches; uniform short-window bound open',
                'source_sha256': hashes, 'finite_controls': controls,
                'plain_division_reference': reference_summary, 'cases': computed}
    if args.regenerate:
        EVIDENCE.write_text(json.dumps(evidence, separators=(',', ':')) + '\n')
    else:
        assert json.loads(EVIDENCE.read_text()) == evidence
    for record in computed:
        verify_witnesses(record)
    assert accelerated['factor'] == 2147483647
    assert (accelerated['base_size'], accelerated['relation_rank'],
            accelerated['norm_queries'], accelerated['trial_divisions'],
            accelerated['maximum_exponent_tested']) == (85, 82, 218738, 5433, 9892)
    assert reference['trial_divisions'] == 18850094
    assert next(r for r in computed if r['window'] == 32)['factor'] is None
    large = computed[-1]
    assert (large['factor'], large['base_size'], large['relation_rank'],
            large['norm_queries'], large['trial_divisions'],
            large['maximum_leading_bits']) == (2305843009213693951, 283, 274,
                                               3650183, 56571, 63)
    print(json.dumps({'status': 'passed', 'finite_controls': controls,
                      'plain_division_reference': reference_summary,
                      'cases': [{key: r[key] for key in
                                 ('integer', 'factor', 'base_size', 'relation_rank',
                                  'norm_queries', 'trial_divisions', 'method')}
                                for r in computed]}, indent=2))


if __name__ == '__main__':
    main()
