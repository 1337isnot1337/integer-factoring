#!/usr/bin/env python3
"""Replay integer-only factoring with deliberately imprimitive smooth norms."""

import argparse
from hashlib import sha256
import json
from math import gcd, prod
from pathlib import Path

from arbitrary_selected_norm_factoring import factor_from_arbitrary_selected_norms
from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form
from prime_form_factor_base import build_prime_form_factor_base, primes_through
from quadratic_ideal_hnf import compose_forms_hnf


def verify(result):
    integer, height = result['integer'], result['height_bound']
    records = build_prime_form_factor_base(-4*integer, result['prime_cutoff'])['records']
    by_prime = {record.prime: record for record in records}
    identity = (1, 0, integer)
    for certificate in result['norm_translations']:
        prime = certificate['generator_prime']
        exponent = certificate['target_exponent']
        form = tuple(certificate['target_form'])
        norm, primitive = certificate['selected_norm'], certificate['primitive_norm']
        content = certificate['integer_content']
        assert prime**(exponent-1) <= height < prime**exponent
        assert form == power(by_prime[prime].positive, exponent, identity, compose_forms_hnf)
        assert norm == content*content*primitive and norm <= height
        assert gcd(norm, 4*integer) == 1
        for vector, expected, expected_content in (
                (certificate['selected_norm_vector'], norm, content),
                (certificate['primitive_vector'], primitive, 1)):
            x, y = vector
            assert gcd(x, y) == expected_content
            assert form[0]*x*x+form[1]*x*y+form[2]*y*y == expected
        assert prod(p**abs(e) for p, e, _ in certificate['prime_terms']) == primitive
        value = identity
        for p, exponent, _ in certificate['prime_terms']:
            generator = by_prime[p].positive
            if exponent < 0:
                generator, exponent = inverse_form(generator), -exponent
            value = compose_forms_hnf(value, power(generator, exponent, identity, compose_forms_hnf))
        assert value == form
    for certificate in result['certificates']:
        row, pivot = certificate['relation_row'], certificate['pivot_column']
        assert records[pivot].prime**abs(row[pivot]) > prod(
            record.prime**abs(exponent) for index, (record, exponent) in enumerate(zip(records, row))
            if index != pivot)
        value = identity
        for record, exponent in zip(records, row):
            generator = record.positive
            if exponent < 0:
                generator, exponent = inverse_form(generator), -exponent
            value = compose_forms_hnf(value, power(generator, exponent, identity, compose_forms_hnf))
        assert value == identity
    if result['factor'] is not None:
        assert 1 < result['factor'] < integer and integer % result['factor'] == 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/arbitrary-selected-norm-factoring.json')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = root/'evidence/height-selected-norm-factoring.json'
    fixtures = json.loads(source.read_text())['records']
    records = []
    for fixture in fixtures:
        integer, cutoff, height = fixture['integer'], fixture['prime_cutoff'], fixture['height_bound']
        inert = next(p for p in primes_through(cutoff) if p > 2 and integer % p
                     and pow((-4*integer) % p, (p-1)//2, p) == p-1)
        selected, padded = [], 0
        for prime, norm in fixture['selected_integers']:
            if norm*inert*inert <= height:
                norm *= inert*inert
                padded += 1
            selected.append((prime, norm))
        assert padded > 0
        result = factor_from_arbitrary_selected_norms(
            integer, cutoff, height, selected, extract_each=False, keep_certificates=True,
        )
        verify(result)
        assert result['factor'] is not None
        assert result['stats']['accepted_relations'] == len(selected)
        assert result['arbitrary_norm_preprocessing']['norms_with_forced_inert_content'] == padded
        record = {
            'integer': integer, 'input_bits': integer.bit_length(),
            'factor': result['factor'], 'prime_cutoff': cutoff, 'height_bound': height,
            'accepted_rows': len(selected), 'padding_inert_prime': inert,
            'selected_norms_requiring_imprimitive_representation': padded,
            'selected_integers': selected,
            'preprocessing_stats': result['arbitrary_norm_preprocessing'],
            'primitive_consumer_stats': result['stats'],
            'extraction_stats': result['last_extraction']['stats'],
            'sample_norm_translation': next(item for item in result['norm_translations']
                                            if item['integer_content'] > 1),
            'all_translations_and_relations_independently_verified': True,
        }
        records.append(record)
        print(json.dumps({key: record[key] for key in (
            'input_bits', 'factor', 'accepted_rows',
            'selected_norms_requiring_imprimitive_representation',
        )}), flush=True)
    control = factor_from_arbitrary_selected_norms(
        103, 31, 7**20, [(7, 9*7**16)], keep_certificates=True,
    )
    verify(control)
    assert control['norm_translations'][0]['integer_content'] == 21
    assert control['norm_translations'][0]['primitive_norm'] == 7**14
    assert control['factor'] is None and control['method'].startswith('no useful top')
    later_gcd = factor_from_arbitrary_selected_norms(
        10807, 31, 21614, [(7, 931), (13, 101)],
        extract_each=False, keep_certificates=True,
    )
    verify(later_gcd)
    assert later_gcd['factor'] == 101 and later_gcd['method'] == 'selected-norm gcd'
    assert later_gcd['arbitrary_norm_preprocessing']['primitive_norms_supplied'] == 1
    assert later_gcd['arbitrary_norm_preprocessing']['proper_candidate_gcds'] == 1
    rejection = factor_from_arbitrary_selected_norms(
        103, 31, 7**20,
        [(7, 7**21), (7, 2), (7, 103), (7, 37), (7, 3),
         (7, 1), (7, 9*7**16), (7, 9*7**11)],
        keep_certificates=True,
    )
    verify(rejection)
    stats = rejection['arbitrary_norm_preprocessing']
    for key in ('over_height_candidates', 'even_candidates', 'noncoprime_candidates',
                'nonsmooth_candidates', 'odd_inert_exponent_candidates',
                'primitive_norms_supplied', 'duplicate_completed_targets'):
        assert stats[key] == 1, key
    assert stats['unrepresented_candidates'] == 2
    malformed = [(2, 31, 10, []), (103, 1, 10, []), (103, 31, 0, []),
                 (103, 31, 7**20, [(2, 1)]), (103, 31, 7**20, [(7, 0)])]
    for call in malformed:
        try:
            factor_from_arbitrary_selected_norms(*call)
        except ValueError:
            pass
        else:
            raise AssertionError('Malformed arbitrary-norm input was accepted.')
    report = {
        'status': 'pass',
        'scope': 'Unsigned smooth-norm processing with forced imprimitive inputs; selection cost remains open.',
        'records': records, 'mixed_content_control': control,
        'proper_gcd_after_valid_relation_control': later_gcd,
        'rejection_and_retry_stats': stats,
        'malformed_input_controls': len(malformed),
        'source_fixture_sha256': sha256(source.read_bytes()).hexdigest(),
        'engine_sha256': sha256((root/'tools/arbitrary_selected_norm_factoring.py').read_bytes()).hexdigest(),
        'decoder_sha256': sha256((root/'tools/arbitrary_norm_representation.py').read_bytes()).hexdigest(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = root/output
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'pass', 'output': str(output)}))


if __name__ == '__main__':
    main()
