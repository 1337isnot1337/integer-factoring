#!/usr/bin/env python3
"""Search finite height-only examples, then decode only their integer norms."""

import argparse
from hashlib import sha256
import json
from math import gcd, prod
from pathlib import Path

from arithmetic_two_group import power
from binary_quadratic_forms import inverse_form
from height_selected_norm_factoring import factor_from_height_selected_norms
from prime_form_factor_base import build_prime_form_factor_base
from quadratic_ideal_hnf import compose_forms_hnf
from smooth_norm_sieve import SmoothNormSieve
from weighted_norm_relations import least_power_above


def verify(result):
    integer = result['integer']
    height = result['height_bound']
    records = build_prime_form_factor_base(-4*integer, result['prime_cutoff'])['records']
    identity = (1, 0, integer)
    weighted_only = 0
    for certificate in result['certificates']:
        prime = certificate['generator_prime']
        exponent = certificate['target_exponent']
        norm = certificate['norm']
        pivot = certificate['pivot_column']
        form = tuple(certificate['target_form'])
        x, y = certificate['primitive_vector']
        assert prime**(exponent-1) <= height < prime**exponent
        assert 1 <= norm <= height and gcd(norm, 4*integer) == 1
        assert gcd(x, y) == 1
        assert form[0]*x*x + form[1]*x*y + form[2]*y*y == norm
        assert power(records[pivot].positive, exponent, identity, compose_forms_hnf) == form
        terms = certificate['prime_terms']
        assert prod(p**abs(e) for p, e, _ in terms) == norm
        row = certificate['relation_row']
        left = prime**abs(row[pivot])
        right = prod(record.prime**abs(e) for index, (record, e) in enumerate(zip(records, row))
                     if index != pivot)
        assert left > right
        weighted_only += abs(row[pivot]) <= sum(abs(e) for index, e in enumerate(row) if index != pivot)
        value = identity
        for record, e in zip(records, row):
            generator = record.positive
            if e < 0:
                generator, e = inverse_form(generator), -e
            value = compose_forms_hnf(value, power(generator, e, identity, compose_forms_hnf))
        assert value == identity
    if result['factor'] is not None:
        assert 1 < result['factor'] < integer and integer % result['factor'] == 0
    return weighted_only


def search_case(integer, cutoff, height_bits):
    height = 1 << height_bits
    records = build_prime_form_factor_base(-4*integer, cutoff)['records']
    sieve = SmoothNormSieve(-4*integer, cutoff)
    selected, searched_witnesses = [], []
    search = {
        'attempts': 0, 'coordinate_points_visited': 0, 'modular_strikes': 0,
        'prime_divisions': 0, 'failed_target_primes': [],
        'generator_orders_searched': False, 'complete_class_group_enumerated': False,
    }
    for record in records:
        exponent = least_power_above(record.prime, height)
        target = power(record.positive, exponent, (1, 0, integer), compose_forms_hnf)
        cap = min(height, 4*target[2])
        while True:
            # Ω(n)<=floor(log2 n) is automatic. This parameter is therefore
            # vacuous on the supplied ellipse and does not narrow the search.
            found = sieve.search(target, cap, height.bit_length())
            search['attempts'] += 1
            search['coordinate_points_visited'] += found.coordinate_points_visited
            search['modular_strikes'] += found.modular_strikes
            search['prime_divisions'] += found.prime_divisions
            if found.norm is not None:
                selected.append((record.prime, found.norm))
                searched_witnesses.append({
                    'generator_prime': record.prime, 'target_exponent': exponent,
                    'norm': found.norm, 'search_vector': list(found.vector),
                    'search_cap': cap,
                })
                break
            if cap == height:
                search['failed_target_primes'].append(record.prime)
                break
            cap = min(height, 4*cap)
    result = factor_from_height_selected_norms(
        integer, cutoff, height, selected, extract_each=False, keep_certificates=True,
    )
    weighted_only = verify(result)
    assert result['factor'] is not None
    assert result['stats']['accepted_relations'] == len(selected)
    assert weighted_only > 0
    return {
        'integer': integer, 'input_bits': integer.bit_length(),
        'factor': result['factor'], 'prime_cutoff': cutoff, 'height_bits': height_bits,
        'height_bound': height, 'base_size': len(records),
        'selected_integers': selected, 'searched_witnesses': searched_witnesses,
        'search': search, 'consumer_stats': result['stats'],
        'extraction_stats': result['last_extraction']['stats'],
        'rows_requiring_weighted_dominance': weighted_only,
        'target_exponent_range': [min(c['target_exponent'] for c in result['certificates']),
                                  max(c['target_exponent'] for c in result['certificates'])],
        'sample_recovered_certificate': result['certificates'][0],
        'all_recovered_relations_independently_replayed': True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--small-only', action='store_true')
    parser.add_argument('--output', default='evidence/height-selected-norm-factoring.json')
    args = parser.parse_args()
    cases = [(1000000016000000063, 1021, 42)]
    if not args.small_only:
        cases.append((4951760154835678088235319297, 4093, 64))
    records = []
    for case in cases:
        record = search_case(*case)
        records.append(record)
        print(json.dumps({key: record[key] for key in (
            'input_bits', 'factor', 'base_size', 'rows_requiring_weighted_dominance',
            'target_exponent_range',
        )}), flush=True)

    high_power = factor_from_height_selected_norms(
        103, 31, 7**20, [(7, 7**16)], keep_certificates=True,
    )
    verify(high_power)
    assert high_power['factor'] is None
    assert high_power['stats']['accepted_relations'] == 1
    assert high_power['stats']['maximum_multiplicity'] == 16
    assert high_power['stats']['maximum_distinct_prime_support'] == 1
    rejection = factor_from_height_selected_norms(
        103, 31, 7**20,
        [(7, 7**21), (7, 103), (7, 2), (7, 1), (7, 7**16), (7, 7**11)],
        keep_certificates=True,
    )
    verify(rejection)
    for key in ['over_height_candidates', 'noncoprime_candidates',
                'nonsmooth_candidates', 'unrepresented_candidates',
                'accepted_relations', 'duplicate_completed_targets']:
        assert rejection['stats'][key] == 1, key
    invalid_inputs = [(2, 31, 7**20, []), (103, 1, 7**20, []),
                      (103, 31, 0, []), (103, 31, 7**20, [(2, 1)]),
                      (103, 31, 7**20, [(7, 0)])]
    for call in invalid_inputs:
        try:
            factor_from_height_selected_norms(*call)
        except ValueError:
            pass
        else:
            raise AssertionError('Malformed height-only input was accepted.')
    root = Path(__file__).resolve().parents[1]
    report = {
        'status': 'pass',
        'scope': 'Finite coordinate searches plus integer-only height decoding; no uniform fast selector proved.',
        'records': records,
        'high_multiplicity_control': high_power,
        'rejection_and_retry_stats': rejection['stats'],
        'malformed_input_controls': len(invalid_inputs),
        'engine_sha256': sha256((root / 'tools/height_selected_norm_factoring.py').read_bytes()).hexdigest(),
        'checker_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'status': 'pass', 'output': str(output)}))


if __name__ == '__main__':
    main()
