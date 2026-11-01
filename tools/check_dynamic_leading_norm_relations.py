"""Exact certificate-only replay for adaptive reduced leading-norm relations."""

import argparse
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd, isqrt, prod
from pathlib import Path

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form, reduce_form, reduced_forms
from certified_small_factor import CertifiedFactorizer, verify_prime_certificates
from dynamic_leading_norm_relations import (factor_with_dynamic_leading_norms,
                                           leading_norm_word)
from quadratic_ideal_hnf import compose_forms_hnf


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'evidence/dynamic-leading-norm-relations.json'
SOURCES = ('tools/certified_small_factor.py',
           'tools/dynamic_leading_norm_relations.py',
           'tools/check_dynamic_leading_norm_relations.py')
CASES = (
    (35, 5, 10, 4), (15, 2, 10, 4), (77, 3, 20, 4),
    (539, 3, 20, 4), (847, 3, 20, 4), (49, 3, 20, 4),
    (97, 3, 20, 4), (21, 2, 20, 4),
    (1022117, 7, 100, 4), (1000036000099, 11, 100, 4),
    (1000000016000000063, 31, 400, 8),
    (4951760154835678088235319297, 127, 80, 8),
)


def canonical(value):
    return json.loads(json.dumps(value))


def rank(matrix):
    rows = [[Fraction(x) for x in row] for row in matrix]
    current = 0
    for column in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(current, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[current], rows[pivot] = rows[pivot], rows[current]
        scale = rows[current][column]
        rows[current] = [x / scale for x in rows[current]]
        for i in range(current + 1, len(rows)):
            scale = rows[i][column]
            rows[i] = [x - scale * y for x, y in zip(rows[i], rows[current])]
        current += 1
        if current == len(rows):
            break
    return current


def evaluate(word, records, identity):
    value = identity
    for prime, exponent in word.items():
        form = tuple(records[int(prime)][1])
        if exponent < 0:
            form, exponent = inverse_form(form), -exponent
        value = compose_forms_hnf(value, power(
            form, exponent, identity, compose_forms_hnf))
    return value


def small_controls():
    factorizer = CertifiedFactorizer()
    for n in range(1, 1201):
        parts = factorizer.factor(n)
        assert prod(p ** e for p, e in parts) == n
        assert all(not any(p % d == 0 for d in range(2, isqrt(p) + 1))
                   for p, _ in parts)
    verified = verify_prime_certificates(factorizer.certificates)
    malformed = deepcopy(factorizer.certificates)
    malformed[3]['witnesses'] = [(2, 1)]
    try:
        verify_prime_certificates(malformed)
    except ValueError:
        pass
    else:
        raise AssertionError('A malformed Lucas witness was accepted.')
    decoded = ramified_two = split_two = proper_gcds = 0
    for integer in range(3, 302, 2):
        delta = -integer if integer % 4 == 3 else -4 * integer
        parity = delta % 2
        identity = (1, parity, (parity - delta) // 4)
        for form in reduced_forms(delta):
            divisor = gcd(form[0], integer)
            if divisor > 1:
                assert divisor < integer
                proper_gcds += 1
                continue
            records = {}
            if integer % 4 == 1:
                records[2] = (2, reduce_form((2, 2, (integer + 1) // 2)))
            terms, new = leading_norm_word(
                integer, form, factorizer.factor(form[0]), records)
            records.update(new)
            assert evaluate(terms, records, identity) == form
            decoded += 1
            if form[0] % 2 == 0:
                ramified_two += integer % 4 == 1
                split_two += integer % 4 == 3
    return {'trial_checked_integers': 1200,
            'trial_checked_prime_certificates': len(verified),
            'decoded_reduced_forms': decoded,
            'ramified_two_norms': ramified_two,
            'split_two_norms': split_two, 'proper_leading_gcds': proper_gcds}


def replay(record):
    verified = verify_prime_certificates(record['prime_certificates'])
    factorizations = {
        int(n): tuple((int(p), int(e)) for p, e in parts)
        for n, parts in record['norm_factorizations'].items()}
    for n, parts in factorizations.items():
        assert all(p in verified and e > 0 for p, e in parts)
        assert len({p for p, _ in parts}) == len(parts)
        assert prod(p ** e for p, e in parts) == n
    computed = factor_with_dynamic_leading_norms(
        record['integer'], record['cutoff'], record['max_steps'],
        lookahead=record['lookahead'], factor_norm=factorizations.__getitem__)
    expected = {k: v for k, v in record.items()
                if k not in ('prime_certificates', 'factoring_stats')}
    assert canonical(computed) == expected

    records = {int(p): value for p, value in record['prime_records'].items()}
    delta = record['discriminant']
    parity = delta % 2
    identity = (1, parity, (parity - delta) // 4)
    active = record['active_primes']
    assert record['maximum_new_active_prime'] <= record['new_active_prime_bound']
    pivots = []
    for row in record['history']:
        pivot = row['pivot_prime']
        assert pivot not in pivots
        pivots.append(pivot)
        assert power(tuple(records[pivot][1]), row['exponent'], identity,
                     compose_forms_hnf) == tuple(row['target'])
        assert evaluate(row['word'], records, identity) == tuple(row['target'])
        assert prod(int(p) ** abs(e) for p, e in row['word'].items()) == row['norm']
        assert row['norm'] <= record['height_bound'] < pivot ** row['exponent']
    for prime, word in record['expressions'].items():
        assert set(map(int, word)).issubset(active)
        assert evaluate(word, records, identity) == tuple(records[int(prime)][1])
    matrix = [[row.get(str(p), 0) for p in active]
              for row in record['relation_rows']]
    assert rank(matrix) == len(matrix) == record['relation_rank']
    for row in record['relation_rows']:
        assert evaluate(row, records, identity) == identity
    factor = record['factor']
    if factor is not None:
        assert 1 < factor < record['integer'] and record['integer'] % factor == 0
    return {'integer': record['integer'], 'factor': factor,
            'method': record['method'], 'steps': record['steps'],
            'active_generators': len(active), 'retained_rank': len(matrix),
            'norm_queries': record['factor_calls'],
            'unique_norms': len(factorizations), 'verified_primes': len(verified)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regenerate', action='store_true')
    args = parser.parse_args()
    controls = small_controls()
    hashes = {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES}
    if args.regenerate:
        records = [factor_with_dynamic_leading_norms(n, cutoff, steps, lookahead=width)
                   for n, cutoff, steps, width in CASES]
        evidence = canonical({'status': 'exact finite controls; no fast closure bound',
                              'source_sha256': hashes, 'small_controls': controls,
                              'cases': records})
        # Certificate records dominate this fixture. Compact JSON keeps
        # complete replay data without several megabytes of whitespace.
        EVIDENCE.write_text(json.dumps(evidence, separators=(',', ':')) + '\n')
    else:
        evidence = json.loads(EVIDENCE.read_text())
        assert evidence['source_sha256'] == hashes
        assert evidence['small_controls'] == controls
    summaries = [replay(record) for record in evidence['cases']]
    successful = next(r for r in evidence['cases']
                      if r['integer'] == 1000000016000000063)
    assert successful['factor'] == 1000000007
    assert (successful['steps'], successful['active_generators'],
            successful['relation_rank'], successful['unit_eliminations'],
            successful['factor_calls'], successful['maximum_norm_bits']) == (276, 8, 8, 268, 2208, 30)
    unresolved = next(r for r in evidence['cases']
                      if r['integer'] == 4951760154835678088235319297)
    assert unresolved['factor'] is None and unresolved['method'] == 'dynamic budget exhausted'
    print(json.dumps({'small_controls': controls, 'replayed_cases': summaries,
                      'status': 'passed'}, indent=2))


if __name__ == '__main__':
    main()
