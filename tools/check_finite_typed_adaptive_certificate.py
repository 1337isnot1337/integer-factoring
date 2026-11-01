#!/usr/bin/env python3
"""Independent finite typed jets, reconstructed pieces, and integer gates."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from check_multipoint_hermite_crossings import local_polynomial, rank
from check_ternary_basis_mixed import transform_seed
from check_weighted_multipoint_hermite import flatten_polynomials, smith_orders, typed_power
from check_weighted_adaptive_catalogues import exact_flag_basis, apply_functional
from finite_typed_adaptive_certificate import (
    amplification_bounds, ceiling_square_root, finite_certificate, gauge_from_heights, read_catalogues,
)
from weighted_adaptive_mixed_place import finite_rank_cost, jet_costs


def exact_profile(seed, dimensions, point, cut, heights, ramification):
    terms = {}
    for key, polynomial in seed.items():
        height = sum(heights[leg][letter] for leg, letter in enumerate(key))
        local = local_polynomial(polynomial, point, 2, 3)
        terms[key] = {ramification * j + height: value for j, value in enumerate(local) if value}
    minimum = min(power for values in terms.values() for power in values)
    clearing = max(0, -minimum)
    table = {}
    for key, values in terms.items():
        if values:
            coefficients = [0] * (max(values) + clearing + 1)
            for power, value in values.items():
                coefficients[power + clearing] = value
            table[key] = coefficients
    return [value - clearing for value in smith_orders(flatten_polynomials(table, dimensions, cut), 3)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root_checks = 0
    for value in range(1001):
        for precision in (0, 3, 16, 64):
            numerator = ceiling_square_root(value, precision)
            target = value << (2 * precision)
            assert numerator * numerator >= target
            assert numerator == 0 or (numerator - 1) ** 2 < target
            root_checks += 1
    sharper = amplification_bounds([(1, 100), (1, 1), (0, 9)], 40)
    assert sharper['rank_sum_gate_passed'] and not sharper['coarse_gate_passed']
    assert sharper['square_root_sum_upper_numerator'] == 11 * (1 << 64)
    assert sharper['zero_piece_bounds'] == 1
    rejections = 0
    for action in (lambda: amplification_bounds([(0, 1)], 4),
                   lambda: amplification_bounds([(1, 1)], 4, (3, 2)),
                   lambda: gauge_from_heights([0] * 6, 2, [[0, .5], [0] * 3, [0] * 3, [0] * 2], (1, 1), (1, 1, 0)),
                   lambda: read_catalogues({'type_counts_x': [1, 1], 'type_counts_y': [1, 1],
                                            'row_places': [{'place': 0, 'gauges': []}],
                                            'column_places': [{'place': 0, 'gauges': []}]})):
        try:
            action()
        except (ValueError, TypeError):
            rejections += 1
        else:
            raise AssertionError('An invalid input was not rejected.')
    basis = [[1, 0, 0], [2, 1, 0], [1, 1, 1]]
    source = transform_seed(basis, basis)
    locations = [(2, 0, (0, 0, 1, 1)), (None, 0, (-1, -1, 0, 0)), (0, 1, (-2, -2, 2, 2))]
    cases = [((1, 1, 0), (0, 1, 1)), ((1, 1, 1), (1, 1, 1)),
             ((2, 1), (1, 1, 1)), ((1, 1, 1), (2, 1))]
    records = []
    totals = {'complete_local_profiles': 0, 'local_channels': 0, 'typed_jet_ranks': 0,
              'reconstructed_piece_ranks': 0, 'endpoint_entries': 0, 'equal_alphabet_agreements': 0}
    for xcounts, ycounts in cases:
        copies = sum(xcounts)
        dimensions = (len(xcounts), len(ycounts), len(ycounts), len(xcounts))
        seed = {key: value for key, value in source.items()
                if all(letter < size for letter, size in zip(key, dimensions))}
        local_records, catalogues = [], []
        for point, cut, slopes in locations:
            group, local = [], []
            for use in (False, True):
                heights = [[(slopes[leg] if use else 0) * letter for letter in range(size)]
                           for leg, size in enumerate(dimensions)]
                orders = exact_profile(seed, dimensions, point, cut, heights, 2)
                gauge = gauge_from_heights(orders, 2, heights, xcounts, ycounts)
                group.append(gauge)
                local.append({'orders': orders, 'heights': heights, 'ramification': 2,
                              'mean_shift': str(gauge.mean_shift)})
                totals['complete_local_profiles'] += 1
                totals['local_channels'] += len(orders)
            catalogues.append(group)
            local_records.append({'point': point, 'cut': cut, 'gauges': local})
        certified = finite_certificate(catalogues[:2], catalogues[2:], copies, 2, xcounts, ycounts)
        if len(xcounts) == len(ycounts):
            previous = finite_rank_cost(catalogues[:2], catalogues[2:], copies, 2, xcounts, ycounts)
            assert previous['row_flag'] == certified['row_flag']
            assert previous['column_flag'] == certified['column_flag']
            assert previous['product'] == certified['amplification']['maximum_rank_product']
            totals['equal_alphabet_agreements'] += 1
        table, typed_dimensions = typed_power(seed, max(dimensions), copies, 3, xcounts, ycounts)
        full_degree = copies * 2
        jets = []
        for (point, cut, _), group in zip(locations, catalogues):
            local = {key: local_polynomial(value, point, full_degree, 3) for key, value in table.items()}
            coefficients = flatten_polynomials(local, typed_dimensions, cut)
            options = [jet_costs(gauge, copies, full_degree) for gauge in group]
            for index in range(full_degree + 1):
                actual = rank(coefficients[index], 3)
                bound = min(certified['cut_dimension'], min(values[index] for values in options))
                assert actual <= bound
                jets.append({'point': point, 'cut': cut, 'order': index, 'rank': actual, 'bound': bound})
                totals['typed_jet_ranks'] += 1
        vectors, weights = exact_flag_basis(certified['row_flag'], certified['column_flag'],
                                             [2, None], [0], full_degree)
        pieces, ranks = [], []
        for index, (vector, weight) in enumerate(zip(vectors, weights)):
            piece = {key: [weight * apply_functional(vector, polynomial, 3) % 3]
                     for key, polynomial in table.items()}
            actual = [rank(flatten_polynomials(piece, typed_dimensions, cut)[0], 3) for cut in (0, 1)]
            assert all(a <= b for a, b in zip(actual, certified['rank_pairs'][index]))
            pieces.append(piece)
            ranks.append({'index': index, 'ranks': actual, 'bounds': certified['rank_pairs'][index]})
            totals['reconstructed_piece_ranks'] += 2
        for key in product(*(range(size) for size in typed_dimensions)):
            expected = int(key[0] == key[3] and key[1] == key[2])
            assert sum(table.get(key, ())) % 3 == expected
            assert sum(piece.get(key, [0])[0] for piece in pieces) % 3 == expected
            totals['endpoint_entries'] += 1
        assert not certified['amplification']['rank_sum_gate_passed']
        records.append({'physical_dimensions': dimensions, 'typed_dimensions': typed_dimensions,
                        'local_profiles': local_records, 'jet_controls': jets, 'piece_controls': ranks,
                        'finite_certificate': certified})
    names = (Path(__file__).name, 'finite_typed_adaptive_certificate.py', 'weighted_adaptive_mixed_place.py',
             'adaptive_mixed_place.py', 'check_ternary_basis_mixed.py', 'check_weighted_multipoint_hermite.py',
             'check_weighted_adaptive_catalogues.py', 'check_multipoint_hermite_crossings.py',
             'check_mixed_place_hermite.py')
    report = {'scope': 'Exact finite typed adaptive arithmetic and real tensor reconstruction controls; no threshold-crossing block is claimed.',
              'status': 'pass', 'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                                for name in names},
              'square_root_upper_bound_controls': root_checks, 'invalid_inputs_rejected': rejections,
              'synthetic_rank_table_arithmetic_control': sharper, 'totals': totals, 'typed_cases': records}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(totals | {'sqrt_controls': root_checks, 'invalid_inputs_rejected': rejections}))


if __name__ == '__main__':
    main()
