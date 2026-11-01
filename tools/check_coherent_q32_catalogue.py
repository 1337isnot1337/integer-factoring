#!/usr/bin/env python3
"""Certify all-scale pressure optimality within the binary q32 boundary class.

The default run recomputes every complete cut profile in the 121-map
catalogue. Each product-pressure difference is factored over the integers
and certified nonnegative in either the power basis or the unnormalized
Bernstein basis. This is a fixed finite boundary-class theorem, not an
improved factoring exponent or a statement over extension fields.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import time

from check_coherent_q32_pressure import counting_polynomial, multiply
from coherent_hot_compounds import validate_maps
from coherent_hot_compounds_fast import profiles
from search_coherent_hot_retractions import _atomic_write, q32_catalogue


ROOT = Path(__file__).resolve().parents[1]


def certificate(reference, candidate):
    first = multiply(counting_polynomial(reference['zero_orders_t']),
                     counting_polynomial(reference['infinity_orders_t']))
    second = multiply(counting_polynomial(candidate['zero_orders_t']),
                      counting_polynomial(candidate['infinity_orders_t']))
    difference = [(first[i] if i < len(first) else 0) - (second[i] if i < len(second) else 0)
                  for i in range(max(len(first), len(second)))]
    while difference and difference[-1] == 0:
        difference.pop()
    if not difference:
        return {'method': 'identical_product_polynomials', 'strict_on_open_unit_interval': False}
    shift = next(i for i, value in enumerate(difference) if value)
    quotient, factors = difference[shift:], 0
    while sum(quotient) == 0:
        running, divided = 0, []
        for value in quotient[:-1]:
            running += value
            divided.append(running)
        assert running + quotient[-1] == 0
        quotient = divided
        factors += 1
    rebuilt = quotient
    for _ in range(factors):
        rebuilt = multiply([1, -1], rebuilt)
    assert [0] * shift + rebuilt == difference
    if all(value >= 0 for value in quotient):
        method, positive_coefficients = 'nonnegative_power_basis', quotient
    else:
        degree = len(quotient) - 1
        positive_coefficients = [sum(quotient[j] * comb(degree - j, i - j)
                                     for j in range(i + 1)) for i in range(degree + 1)]
        assert all(value >= 0 for value in positive_coefficients)
        # Q(r)=sum_i b_i*r^i*(1-r)^(degree-i), with no binomial
        # factor absorbed into the displayed b_i. Expand independently.
        expanded = [sum(positive_coefficients[i] * (-1) ** (j - i)
                        * comb(degree - i, j - i) for i in range(j + 1))
                    for j in range(degree + 1)]
        assert expanded == quotient
        method = 'nonnegative_unnormalized_bernstein_basis'
    assert any(value > 0 for value in positive_coefficients)
    return {'method': method, 'strict_on_open_unit_interval': True,
            'r_factor_power': shift, 'one_minus_r_factor_power': factors,
            'quotient_power_coefficients': quotient,
            'nonnegative_basis_coefficients': positive_coefficients}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path,
                        default=ROOT / 'evidence/coherent-hot-q32-retractions.json')
    parser.add_argument('--polynomials-only', action='store_true',
                        help='Check saved exact profiles without a native full-profile replay.')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    saved = json.loads(args.catalogue.read_text())
    assert saved['completed'] and saved['finite_list_exhausted'] and saved['unresolved'] == 0
    assert len(saved['records']) == 121
    winner = next(record for record in saved['records']
                  if record['description'].get('kernel') == 7
                  and record['description'].get('covector') == 15)
    assert winner['index'] == 55
    names = (Path(__file__).name, 'check_coherent_q32_pressure.py',
             'coherent_hot_compounds.py', 'coherent_hot_compounds_fast.py',
             'coherent_hot_compounds_fast.cpp', 'search_coherent_hot_retractions.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'scope': __doc__, 'completed': False, 'winner_index': 55,
              'all_complete_profiles_recomputed': not args.polynomials_only,
              'catalogue_sha256': sha256(args.catalogue.read_bytes()).hexdigest(),
              'started_utc': datetime.now(timezone.utc).isoformat(), 'records': [],
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    started, counts = time.monotonic(), Counter()
    _atomic_write(args.output, report)
    recipes = list(q32_catalogue())
    assert len(recipes) == len(saved['records'])
    for index, (description, embedding, retraction, modes) in enumerate(recipes):
        record = saved['records'][index]
        assert record['index'] == index and record['description'] == description
        assert record['embedding'] == embedding and record['retraction'] == retraction
        _, _, heights, _, checks = validate_maps(embedding, retraction, modes)
        assert sorted(heights) == sorted(winner['heights'])
        counts['deletion_checks'] += checks
        if not args.polynomials_only:
            recomputed = profiles(embedding, retraction, modes)
            for key in ('zero_orders_t', 'infinity_orders_t'):
                assert record[key] == recomputed[key]
                counts['complete_profiles_recomputed'] += 1
                counts['complete_smith_orders_recomputed'] += len(record[key])
        proof = certificate(record, winner)
        counts[proof['method']] += 1
        report['records'].append({'index': index, 'description': description, **proof})
        report['elapsed_seconds'] = time.monotonic() - started
        report['counts'] = dict(counts)
        _atomic_write(args.output, report)
        if index % 20 == 0:
            print(json.dumps({'completed_maps': index + 1, 'elapsed_seconds': report['elapsed_seconds']}),
                  flush=True)
    assert counts['nonnegative_power_basis'] == 82
    assert counts['nonnegative_unnormalized_bernstein_basis'] == 37
    assert counts['identical_product_polynomials'] == 2
    report.update(completed=True, all_121_maps_compared=True,
                  all_scale_catalogue_optimality_proved=True,
                  optimality_interval='0 < r < 1',
                  physical_height_multisets_identical=True,
                  floating_point_pressure_used=False,
                  factoring_exponent_improved=False,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    _atomic_write(args.output, report)
    print(json.dumps({'completed': True, 'counts': dict(counts),
                      'equal_product_indices': [r['index'] for r in report['records']
                                                if r['method'] == 'identical_product_polynomials'],
                      'elapsed_seconds': report['elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()
