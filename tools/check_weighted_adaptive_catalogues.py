#!/usr/bin/env python3
"""Exact weighted jet catalogues, typed ranks, and opposite-flag reconstruction."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from weighted_adaptive_mixed_place import GaugeProfile, jet_costs, catalogue_flag, finite_rank_cost
from check_ternary_basis_mixed import transform_seed, ramified_profile
from check_weighted_multipoint_hermite import typed_power, flatten_polynomials
from check_multipoint_hermite_crossings import local_polynomial, rank, solve
from check_mixed_place_hermite import jet_vector, nullspace, combine, apply_functional


def arithmetic_controls():
    gauges = [GaugeProfile((0, 0, 2, 4), 1, Fraction(0)),
              GaugeProfile((-3, -1, 0, 2), 2, Fraction(-1)),
              GaugeProfile((1, 1, 1, 1), 1, Fraction(0)),
              GaugeProfile((0, 0, 5, 5), 1, Fraction(0))]
    coefficients = 0
    for gauge in gauges:
        for copies in range(1, 5):
            actual = jet_costs(gauge, copies, 12)
            counts = Counter(map(sum, product(gauge.orders, repeat=copies)))
            expected = [sum(count * max(gauge.ramification * j + copies * gauge.mean_shift - total + 1, 0)
                            for total, count in counts.items()) for j in range(13)]
            assert actual == expected
            translated = GaugeProfile(tuple(value + 7 for value in gauge.orders),
                                      gauge.ramification, gauge.mean_shift + 7)
            assert jet_costs(translated, copies, 12) == actual
            coefficients += len(actual)
    catalogue = [[gauges[2], gauges[3]], [gauges[0]]]
    flag = catalogue_flag(catalogue, 2, 8, 1000000)
    costs = [[min(jet_costs(g, 2, 8)[j] for g in group) for j in range(9)] for group in catalogue]
    for length, entry in enumerate(flag, 1):
        best = min(sum(costs[0][:first]) + sum(costs[1][:length - first])
                   for first in range(length + 1))
        assert entry['prefix_cost'] == best
    assert len({jet['gauge_index'] for jet in flag if jet['place_index'] == 0}) > 1
    return {'word_sum_coefficients': coefficients,
            'translated_gauge_coefficients': coefficients,
            'initial_segment_prefixes': len(flag), 'gauge_switch_control': flag}


def exact_flag_basis(rows, columns, row_points, column_points, degree, prime=3):
    assert set(row_points).isdisjoint(column_points)
    first = [jet_vector(degree, row_points[jet['place_index']], jet['order'], prime) for jet in rows]
    second = [jet_vector(degree, column_points[jet['place_index']], jet['order'], prime) for jet in columns]
    assert rank(first, prime) == rank(second, prime) == degree + 1
    vectors = []
    for k in range(degree + 1):
        a, b = first[:k + 1], second[:degree - k + 1]
        equations = [[v[j] for v in a] + [-v[j] % prime for v in b] for j in range(degree + 1)]
        kernel = nullspace(equations, prime)
        assert len(kernel) == 1
        vector = combine(kernel[0][:k + 1], a, prime)
        assert vector == combine(kernel[0][k + 1:], b, prime)
        assert any(vector)
        vectors.append(vector)
    assert rank(vectors, prime) == degree + 1
    weights = solve(list(map(list, zip(*vectors))), [1] * (degree + 1), prime)
    assert combine(weights, vectors, prime) == [1] * (degree + 1)
    return vectors, weights


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    arithmetic = arithmetic_controls()
    basis = [[1, 0, 0], [2, 1, 0], [1, 1, 1]]
    seed = transform_seed(basis, basis)
    half = Fraction(1, 2)
    locations = [(2, 0, [(0, 0, 0, 0), (0, 0, half, half), (0, half, 0, half), (1, -1, 0, 0)]),
                 (None, 0, [(0, 0, 0, 0), (-half, -half, 0, 0), (-half, 0, -half, 0), (0, 0, 1, -1)]),
                 (0, 1, [(0, 0, 0, 0), (-1, -1, 1, 1), (-half, -half, half, half), (1, 0, 0, -1)])]
    local_profiles = []
    for point, cut, slopes_list in locations:
        records = []
        for slopes in slopes_list:
            orders, clearing = ramified_profile(seed, point, cut, slopes)
            records.append({'slopes': [str(Fraction(s)) for s in slopes],
                            'orders': orders, 'ramification': 2, 'clearing': clearing})
        local_profiles.append({'point': point, 'cut': cut, 'gauges': records})
    typed_records = []
    for xcounts, ycounts in [((1, 1, 0), (0, 1, 1)), ((1, 1, 1), (1, 1, 1))]:
        copies = sum(xcounts)
        mean_x = sum(Fraction(a * count, copies) for a, count in enumerate(xcounts))
        mean_y = sum(Fraction(a * count, copies) for a, count in enumerate(ycounts))
        means = (mean_x, mean_y, mean_y, mean_x)
        catalogues = []
        for location in local_profiles:
            gauges = []
            for record in location['gauges']:
                shift = sum(2 * Fraction(s) * average for s, average in zip(record['slopes'], means))
                gauges.append(GaugeProfile(tuple(record['orders']), 2, shift))
            catalogues.append(gauges)
        bound = finite_rank_cost(catalogues[:2], catalogues[2:], copies, 2, xcounts, ycounts)
        table, dimensions = typed_power(seed, 3, copies, 3, xcounts, ycounts)
        full_degree = 2 * copies
        jet_checks, records = 0, []
        for location, gauges in zip(local_profiles, catalogues):
            local = {key: local_polynomial(value, location['point'], full_degree, 3)
                     for key, value in table.items()}
            coefficients = flatten_polynomials(local, dimensions, location['cut'])
            alternatives = [jet_costs(gauge, copies, full_degree) for gauge in gauges]
            for j in range(full_degree + 1):
                actual = rank(coefficients[j], 3)
                predicted = min(bound['cut_dimension'], min(row[j] for row in alternatives))
                assert actual <= predicted
                records.append({'point': location['point'], 'cut': location['cut'], 'jet': j,
                                'actual_rank': actual, 'catalogue_bound': predicted,
                                'individual_gauge_costs': [row[j] for row in alternatives]})
                jet_checks += 1
        vectors, weights = exact_flag_basis(bound['row_flag'], bound['column_flag'], [2, None], [0], full_degree)
        pieces, rank_records = [], []
        for k, (vector, weight) in enumerate(zip(vectors, weights)):
            piece = {key: [weight * apply_functional(vector, value, 3) % 3]
                     for key, value in table.items()}
            ranks = [rank(flatten_polynomials(piece, dimensions, cut)[0], 3) for cut in (0, 1)]
            limits = [bound['row_flag'][k]['prefix_cost'], bound['column_flag'][full_degree - k]['prefix_cost']]
            assert all(actual <= limit for actual, limit in zip(ranks, limits))
            pieces.append(piece)
            rank_records.append({'index': k, 'ranks': ranks, 'bounds': limits})
        endpoints = 0
        for key in product(*(range(d) for d in dimensions)):
            x, y, u, v = key
            expected = int(x == v and y == u)
            assert sum(table.get(key, ())) % 3 == expected
            assert sum(piece.get(key, [0])[0] for piece in pieces) % 3 == expected
            endpoints += 1
        typed_records.append({'copies': copies, 'xcounts': xcounts, 'ycounts': ycounts,
                              'dimensions': dimensions, 'mean_shifts': [[str(g.mean_shift) for g in group] for group in catalogues],
                              'jet_rank_checks': jet_checks, 'jets': records,
                              'adaptive_bound': bound, 'piece_rank_checks': 2 * len(rank_records),
                              'pieces': rank_records, 'endpoint_entries': endpoints})
    # Omitting the physical mean shift falsely predicts a zero finite-place
    # constant coefficient, whose exact uniform-type rank is one.
    orders = tuple(local_profiles[0]['gauges'][1]['orders'])
    assert jet_costs(GaugeProfile(orders, 2, Fraction(2)), 3, 6)[0] == 1
    assert jet_costs(GaugeProfile(orders, 2, Fraction(0)), 3, 6)[0] == 0
    uniform_zero = next(row for row in typed_records[1]['jets']
                        if row['point'] == 2 and row['cut'] == 0 and row['jet'] == 0)
    assert uniform_zero['actual_rank'] == 1
    directory = Path(__file__).parent
    sources = (Path(__file__).name, 'weighted_adaptive_mixed_place.py', 'adaptive_mixed_place.py',
               'check_ternary_basis_mixed.py', 'check_weighted_multipoint_hermite.py',
               'check_mixed_place_hermite.py', 'check_multipoint_hermite_crossings.py')
    result = {'scope': 'Exact weighted catalogue jet bounds and typed opposite-flag reconstruction; no improved asymptotic rho or factoring score claimed.',
              'source_sha256': {name: sha256((directory / name).read_bytes()).hexdigest() for name in sources},
              'arithmetic_controls': arithmetic, 'physical_basis': basis, 'local_profiles': local_profiles,
              'typed_cases': typed_records, 'missing_mean_shift_rejected': True,
              'totals': {'exact_local_profiles': sum(len(row['gauges']) for row in local_profiles),
                         'typed_jet_rank_checks': sum(row['jet_rank_checks'] for row in typed_records),
                         'piece_rank_checks': sum(row['piece_rank_checks'] for row in typed_records),
                         'endpoint_entries': sum(row['endpoint_entries'] for row in typed_records)}}
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['totals']))


if __name__ == '__main__':
    main()
