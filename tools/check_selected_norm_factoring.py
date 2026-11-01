#!/usr/bin/env python3
"""Replay factor extraction after withholding coordinates and orientations."""

import argparse
import hashlib
import json
from math import gcd, prod
from pathlib import Path

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from selected_norm_factoring import factor_from_selected_norms


def verify_certificates(result):
    integer = result['integer']
    discriminant = -4 * integer
    base = build_prime_form_factor_base(discriminant, result['prime_cutoff'])['records']
    generators = tuple(record.positive for record in base)
    identity = (1, 0, integer)
    for certificate in result['certificates']:
        target = tuple(certificate['target_form'])
        norm = certificate['norm']
        x, y = certificate['primitive_vector']
        assert gcd(x, y) == 1 and gcd(norm, discriminant) == 1
        assert target[0]*x*x + target[1]*x*y + target[2]*y*y == norm
        terms = certificate['prime_terms']
        assert prod(prime**abs(exponent) for prime, exponent, _ in terms) == norm
        assert sum(abs(exponent) for _, exponent, _ in terms) <= result['multiplicity_budget']
        value = identity
        row = certificate['relation_row']
        pivot = certificate['pivot_column']
        assert abs(row[pivot]) > sum(abs(entry) for index, entry in enumerate(row)
                                     if index != pivot)
        for generator, exponent in zip(generators, row):
            if exponent < 0:
                generator, exponent = inverse_form(generator), -exponent
            value = compose_forms_hnf(value, power(generator, exponent, identity,
                                                   compose_forms_hnf))
        assert value == identity
        assert power(generators[pivot], result['diagonal_parameter'], identity,
                     compose_forms_hnf) == target
    if result['factor'] is not None:
        assert 1 < result['factor'] < integer and integer % result['factor'] == 0


def public_summary(result):
    return {
        'integer': result['integer'], 'input_bits': result['integer'].bit_length(),
        'factor': result['factor'], 'method': result['method'],
        'prime_cutoff': result['prime_cutoff'],
        'multiplicity_budget': result['multiplicity_budget'],
        'factor_base_size': len(result['factor_base_primes']),
        'accepted_target_count': len(result['accepted_target_primes']),
        'stats': result['stats'],
        'last_extraction_stats': (result['last_extraction']['stats']
                                  if result['last_extraction'] else None),
        'sample_recovered_certificate': (result['certificates'][0]
                                         if result['certificates'] else None),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/selected-norm-factoring.json')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = json.loads((root / 'evidence/short-norm-factoring-certificates.json').read_text())
    records = []
    for stored in source['records']:
        # These integers are the ONLY witness information passed across
        # the selector interface. No old factorization, coordinate, word,
        # class number, relation determinant, or annihilator is supplied.
        candidates = tuple((witness['generator_prime'], witness['represented_norm'])
                           for witness in stored['witnesses'])
        result = factor_from_selected_norms(
            stored['integer'], stored['smooth_bound'], stored['word_budget'],
            candidates, extract_each=False, keep_certificates=True,
        )
        assert result['factor'] is not None
        assert result['stats']['accepted_relations'] == len(candidates)
        verify_certificates(result)
        records.append(public_summary(result))

    stored = source['records'][0]
    prefix = factor_from_selected_norms(
        stored['integer'], stored['smooth_bound'], stored['word_budget'],
        ((w['generator_prime'], w['represented_norm']) for w in stored['witnesses']),
        extract_each=True, keep_certificates=True,
    )
    assert prefix['factor'] is not None
    assert prefix['stats']['accepted_relations'] < len(stored['witnesses'])
    verify_certificates(prefix)

    # One target is first given four different invalid candidates and then
    # a valid one. The final duplicate cannot add a dependent relation row.
    rejection = factor_from_selected_norms(
        103, 31, 1,
        [(7, 1), (7, 49), (7, 2), (7, 103), (7, 13), (7, 19)],
        keep_certificates=True,
    )
    assert rejection['factor'] is None
    for name in ['unrepresented_candidates', 'over_budget_candidates',
                 'nonsmooth_candidates', 'noncoprime_candidates',
                 'accepted_relations', 'duplicate_completed_targets']:
        assert rejection['stats'][name] == 1, name
    verify_certificates(rejection)

    gcd_case = factor_from_selected_norms(103*107, 3, 1, [(3, 103)])
    trial_case = factor_from_selected_norms(15, 5, 1, [])
    assert gcd_case['factor'] == 103 and gcd_case['method'] == 'selected-norm gcd'
    assert trial_case['factor'] == 3 and trial_case['method'] == 'factor-base trial division'
    invalid_inputs = [
        (2, 31, 1, []), (103, 1, 1, []), (103, 31, -1, []),
        (103, 31, 1, [(2, 1)]), (103, 31, 1, [(7, 0)]),
    ]
    for call in invalid_inputs:
        try:
            factor_from_selected_norms(*call)
        except ValueError:
            pass
        else:
            raise AssertionError('malformed input was accepted')

    report = {
        'status': 'pass',
        'scope': 'Consumes selected integers; no bound for choosing them is claimed.',
        'selected_large_norms_decoded': sum(r['accepted_target_count'] for r in records),
        'large_records': records,
        'partial_prefix_record': public_summary(prefix),
        'rejection_and_retry_control': public_summary(rejection),
        'direct_gcd_controls': 2,
        'malformed_input_controls': len(invalid_inputs),
        'source_fixture': 'short-norm-factoring-certificates.json',
        'engine_sha256': hashlib.sha256((root / 'tools/selected_norm_factoring.py').read_bytes()).hexdigest(),
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({
        'status': report['status'],
        'selected_large_norms_decoded': report['selected_large_norms_decoded'],
        'large_factors': [r['factor'] for r in records],
        'partial_prefix_relations': prefix['stats']['accepted_relations'],
        'full_relation_dimension': len(stored['witnesses']),
        'output': str(output),
    }, indent=2))


if __name__ == '__main__':
    main()
