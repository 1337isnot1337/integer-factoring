#!/usr/bin/env python3
"""Independent reconstruction and complete-cost replay of the restriction search."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path

from scipy.special import logsumexp

from affine_place_profile_oracle import polynomial
from check_finite_fiber_closed_places import (
    evaluate_matrix, finite_field_tables, numeric_determinant,
)
from check_fused_closed_profiles import from_canonical, multiply, quotient_rank
from f4_closed_place_oracle import F4ClosedPlaceOracle, checked_factor
from f16_identity_operator_profiles import evaluate_one
from fused_six_vertex_crossings import pscale
from physical_split_crossings import flatten
from reflection_laurent_physical_frames import digest


def key(row):
    return row['source_index'], tuple(row['subset_x']), tuple(row['subset_y'])


def paired_types(details):
    result = []
    for name in ('type_x', 'type_y'):
        values = [Fraction(str(v)) for v in details[name]]
        assert min(values) > 0
        total = sum(values)
        result.append([value / total for value in values])
    return result


def height_sides(seed, params):
    dimensions = (seed['q_x'], seed['q_u'], seed['q_u'], seed['q_x'])
    stride = sum(q - 1 for q in dimensions)
    result = []
    for side in range(2):
        values, offset = [], side * stride
        for q in dimensions:
            values.append(Fraction(0))
            values.extend(Fraction(str(v)) for v in params[offset:offset + q - 1])
            offset += q - 1
        result.append(values)
    return result


def replay_cost(seed, best, asymmetric):
    params = best['parameters']
    if asymmetric:
        angle = params[-1]
        params = params[:-1]
        coefficients = (math.cos(angle), math.sin(angle))
        expected = best['support_value']
    else:
        coefficients = (1.0, 1.0)
        expected = best['value']
    heights = height_sides(seed, params)
    scale = Fraction(str(params[-1]))
    types = paired_types(best['details'])
    entropy = -sum(float(v) * math.log(float(v)) for row in types for v in row)
    oracle = F4ClosedPlaceOracle(seed)
    total = seed['degree'] * float(scale)
    profiles = []
    for side, (cut, place) in enumerate(((1, None), (0, (0, 1)))):
        orders = oracle.profile(cut, place, heights[side], scale)
        coefficient = coefficients[side]
        total += coefficient * float(logsumexp([-float(v) / coefficient for v in orders]))
        shift = oracle.type_shift(heights[side], types[0], types[1])
        total += float(shift)
        profiles.append({'cut': cut, 'place': place, 'orders': [str(v) for v in orders],
                         'type_shift': str(shift), 'branch_weight': coefficient})
    observed = total / (entropy if asymmetric else 2 * entropy)
    assert abs(observed - expected) < 3e-8, (observed, expected)
    return {'recorded_cost': expected, 'recomputed_cost': observed,
            'types': [[str(v) for v in row] for row in types], 'profiles': profiles}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    paths = {name: Path('evidence/' + filename) for name, filename in (
        ('input', 'deflated-finite-fiber-input.json'),
        ('search', 'finite-fiber-restriction-search.json'),
        ('catalogue', 'finite-fiber-restriction-catalogue.json'),
        ('selected', 'finite-fiber-restriction-selected-input.json'),
        ('asymmetric', 'finite-fiber-restriction-asymmetric-search.json'),
    )}
    raw = {key: path.read_bytes() for key, path in paths.items()}
    reports = {key: json.loads(value) for key, value in raw.items()}
    search, catalogue = reports['search'], reports['catalogue']
    assert search['completed'] and catalogue['completed'] and reports['asymmetric']['completed']
    assert search['input_sha256'] == catalogue['input_sha256'] == sha256(raw['input']).hexdigest()
    assert reports['asymmetric']['input_sha256'] == sha256(raw['selected']).hexdigest()
    for name in ('search', 'catalogue', 'asymmetric'):
        for source, expected in reports[name]['source_sha256'].items():
            assert sha256(Path(__file__).with_name(source).read_bytes()).hexdigest() == expected
    searched = {key(row): row for row in search['records']}
    algebra = {key(row): row for row in catalogue['records']}
    assert len(searched) == len(algebra) == 200 and set(searched) == set(algebra)
    field, multiplication, inverse, embedding = finite_field_tables()
    coefficient_entries = endpoint_entries = determinant_evaluations = 0
    for identity, record in searched.items():
        source_index, first, second = identity
        original = reports['input']['cases'][source_index]['seed']
        seed = record['seed']
        qx, qy = len(first), len(second)
        common = {int(k): c for k, c in record['scalar_content'].items()}
        endpoint = evaluate_one(common)
        assert endpoint and digest(seed['values']) == record['tensor_sha256']
        assert record['tensor_sha256'] == algebra[identity]['tensor_sha256']
        for index, value in enumerate(seed['values']):
            x = index % qx
            y = index // qx % qy
            u = index // (qx * qy) % qy
            v = index // (qx * qy * qy)
            source = first[x] + 4 * (second[y] + 4 * (second[u] + 4 * first[v]))
            # Independent direct coordinate restriction and reconstruction
            # of the removed scalar; no call to retract() or primitive_seed().
            assert multiply(polynomial(value), common) == pscale(polynomial(original['values'][source]), endpoint)
            assert evaluate_one(polynomial(value)) == int(x == v and y == u)
            coefficient_entries += 1
            endpoint_entries += 1
        n = qx * qy
        for cut in (0, 1):
            matrix = [[polynomial(value) for value in row] for row in flatten(seed, cut)[0]]
            determinant = {int(k): value for k, value in algebra[identity]['determinants'][cut].items()}
            reconstructed = {0: determinant[max(determinant)]}
            for factor in algebra[identity]['factorizations'][cut]:
                checked_factor(tuple(factor['factor']))
                for _ in range(factor['multiplicity']):
                    reconstructed = multiply(reconstructed, from_canonical(factor['factor']))
            assert reconstructed == determinant
            for point in range(n * seed['degree'] + 1):
                evaluated = evaluate_matrix(matrix, point, multiplication, embedding, seed['degree'])
                actual = numeric_determinant(evaluated, multiplication, inverse)
                expected = 0
                for exponent, coefficient in determinant.items():
                    expected ^= multiplication[embedding[coefficient]][field.power(point, exponent)]
                assert actual == expected
                determinant_evaluations += 1
    # Representative quotient-ring checks cover every residue degree
    # present, independently of the bordered-minor profile algorithm.
    quotient_controls = 0
    for degree in range(1, 8):
        choices = [(len(row['subset_x']) * len(row['subset_y']), row, place)
                   for row in catalogue['records'] for place in row['closed_points']
                   if place['residue_degree'] == degree and place['factor'] != [0, 1]]
        _, record, place = min(choices, key=lambda item: item[0])
        seed = searched[key(record)]['seed']
        factor = from_canonical(place['factor'])
        for cut in (0, 1):
            matrix = [[polynomial(value) for value in row] for row in flatten(seed, cut)[0]]
            orders = place['orders'][cut]
            for length in range(1, max(orders) + 2):
                observed = quotient_rank(matrix, factor, length)
                assert observed == degree * sum(max(length - order, 0) for order in orders)
                quotient_controls += 1
    controls = []
    for label, record in sorted(search['winners'].items()):
        controls.append({'mode': 'balanced', 'label': label, 'tensor_sha256': record['tensor_sha256'],
                         **replay_cost(record['seed'], record['best'], False)})
    selected = {case['label']: case for case in reports['selected']['cases']}
    for label, case in selected.items():
        matching = searched[key(case)]
        assert case['tensor_sha256'] == matching['tensor_sha256']
        assert case['initial_params'] == matching['best']['parameters']
        assert case['initial_balanced_rho'] == matching['best']['value']
    for label, record in sorted(reports['asymmetric']['winners'].items()):
        case = selected[label]
        assert record['tensor_sha256'] == case['tensor_sha256']
        controls.append({'mode': 'asymmetric', 'label': label, 'tensor_sha256': case['tensor_sha256'],
                         **replay_cost(case['seed'], record['best'], True)})
    result = {
        'scope': 'Complete restriction reconstruction and determinant catalogues, plus exact rational-order winner replay. Final costs are floating and no global exclusion is claimed.',
        'input_sha256': {name: sha256(value).hexdigest() for name, value in raw.items()},
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in (Path(__file__).name, 'f4_closed_place_oracle.py',
                                       'check_finite_fiber_closed_places.py', 'check_fused_closed_profiles.py')},
        'counts': {'complete_restrictions': 200, 'reconstructed_polynomial_entries': coefficient_entries,
                   'strict_endpoint_entries': endpoint_entries,
                   'independent_determinant_evaluations': determinant_evaluations,
                   'independent_quotient_ranks': quotient_controls,
                   'exact_complete_winner_profiles': 2 * len(controls)},
        'balanced_trial_summary': {'attempted': len(search['records']),
                                   'optimizer_runs_completed': search['completed_trials'],
                                   'started_utc': search['started_utc'], 'finished_utc': search['finished_utc'],
                                   'elapsed_seconds': search['elapsed_seconds'], **search['counters']},
        'asymmetric_trial_summary': {'attempted': len(reports['asymmetric']['records']),
                                     'optimizer_runs_completed': sum(row['completed'] for row in reports['asymmetric']['records']),
                                     'elapsed_seconds': reports['asymmetric']['elapsed_seconds']},
        'controls': controls, 'passed': True,
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'counts': result['counts'], 'passed': True}))


if __name__ == '__main__':
    main()
