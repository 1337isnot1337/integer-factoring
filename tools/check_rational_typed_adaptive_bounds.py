#!/usr/bin/env python3
"""Exact coefficient-majorant, free-prefix, and typed reconstruction checks."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import time

from check_finite_typed_adaptive_certificate import exact_profile
from check_multipoint_hermite_crossings import rank
from check_ternary_basis_mixed import transform_seed
from check_weighted_adaptive_catalogues import exact_flag_basis, apply_functional
from check_weighted_multipoint_hermite import flatten_polynomials, typed_power
from finite_typed_adaptive_certificate import gauge_from_heights
from rational_typed_adaptive_bounds import finite_rational_certificate, rational_jet_costs
from weighted_adaptive_mixed_place import GaugeProfile, jet_costs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    parameters = (Fraction(1, 2), Fraction(2, 3), Fraction(4, 5))
    profiles = [GaugeProfile((0, 0, 2, 4), 1, Fraction(0)),
                GaugeProfile((-3, -1, 0, 2), 2, Fraction(-1)),
                GaugeProfile((1, 1, 1, 1), 1, Fraction(0)),
                GaugeProfile((0, 0, 5, 5), 3, Fraction(1, 2))]
    coefficient_checks = translated_checks = 0
    for profile in profiles:
        for copies in (2, 4, 6, 8):
            actual = jet_costs(profile, copies, 18)
            for parameter in parameters:
                for cap in (31, 10**12):
                    bound = rational_jet_costs(profile, copies, 18, parameter, cap)
                    assert all(min(value, cap) <= estimate <= cap for value, estimate in zip(actual, bound))
                    shift = GaugeProfile(tuple(value + 11 for value in profile.orders), profile.ramification,
                                         profile.mean_shift + 11)
                    assert rational_jet_costs(shift, copies, 18, parameter, cap) == bound
                    coefficient_checks += len(bound)
                    translated_checks += len(bound)
    boundary = rational_jet_costs(GaugeProfile((1, 1), 1, Fraction(0)), 3, 5, Fraction(1, 2), 10000)
    assert boundary[:3] == [0, 0, 0] and boundary[3] == 32
    rejected = 0
    for parameter in (Fraction(0), Fraction(1), Fraction(-1, 2), Fraction(3, 2)):
        try:
            rational_jet_costs(profiles[0], 2, 4, parameter, 100)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError('An invalid rational parameter was accepted.')
    basis = [[1, 0, 0], [2, 1, 0], [1, 1, 1]]
    original = transform_seed(basis, basis)
    locations = [(2, 0, (0, 0, 1, 1)), (None, 0, (-1, -1, 0, 0)), (0, 1, (-2, -2, 2, 2))]
    records, local_count, piece_count, endpoint_count = [], 0, 0, 0
    uniform_catalogues = None
    for xcounts, ycounts in [((1, 1, 0), (0, 1, 1)), ((1, 1, 1), (1, 1, 1)),
                             ((2, 1), (1, 1, 1)), ((1, 1, 1), (2, 1))]:
        copies = sum(xcounts)
        dimensions = (len(xcounts), len(ycounts), len(ycounts), len(xcounts))
        seed = {key: value for key, value in original.items()
                if all(letter < size for letter, size in zip(key, dimensions))}
        catalogues = []
        for point, cut, slopes in locations:
            place = []
            for use in (False, True):
                heights = [[(slopes[leg] if use else 0) * letter for letter in range(size)]
                           for leg, size in enumerate(dimensions)]
                orders = exact_profile(seed, dimensions, point, cut, heights, 2)
                profile = gauge_from_heights(orders, 2, heights, xcounts, ycounts)
                place.append((profile, parameters))
                local_count += 1
            catalogues.append(place)
        if xcounts == ycounts == (1, 1, 1):
            uniform_catalogues = catalogues
        result = finite_rational_certificate(catalogues[:2], catalogues[2:], copies, 2, xcounts, ycounts)
        table, typed_dimensions = typed_power(seed, max(dimensions), copies, 3, xcounts, ycounts)
        vectors, weights = exact_flag_basis(result['row_flag'], result['column_flag'], [2, None], [0], 2 * copies)
        pieces, observed = [], []
        for index, (vector, weight) in enumerate(zip(vectors, weights)):
            piece = {key: [weight * apply_functional(vector, polynomial, 3) % 3]
                     for key, polynomial in table.items()}
            ranks = [rank(flatten_polynomials(piece, typed_dimensions, cut)[0], 3) for cut in (0, 1)]
            assert all(a <= b for a, b in zip(ranks, result['rank_pairs'][index]))
            pieces.append(piece)
            observed.append({'index': index, 'actual_ranks': ranks, 'bounds': result['rank_pairs'][index]})
            piece_count += 2
        for key in product(*(range(size) for size in typed_dimensions)):
            expected = int(key[0] == key[3] and key[1] == key[2])
            assert sum(piece.get(key, [0])[0] for piece in pieces) % 3 == expected
            endpoint_count += 1
        assert not result['amplification']['rank_sum_gate_passed']
        records.append({'type_counts_x': xcounts, 'type_counts_y': ycounts,
                        'typed_dimensions': typed_dimensions, 'actual_pieces': observed,
                        'finite_rational_certificate': result})
    started = time.monotonic()
    large = finite_rational_certificate(uniform_catalogues[:2], uniform_catalogues[2:], 600, 2,
                                        (200, 200, 200), (200, 200, 200))
    assert not large['amplification']['rank_sum_gate_passed']
    rank_table = json.dumps(large['rank_pairs'], separators=(',', ':')).encode()
    large_summary = {'copies': 600, 'tensor_degree': 1200,
                     'external_label_count_bits': [value.bit_length() for value in large['external_label_counts']],
                     'maximum_rank_product_bits': large['amplification']['maximum_rank_product'].bit_length(),
                     'rank_pairs': len(large['rank_pairs']), 'rank_table_sha256': sha256(rank_table).hexdigest(),
                     'rank_sum_gate_passed': False, 'seconds': time.monotonic() - started}
    names = (Path(__file__).name, 'rational_typed_adaptive_bounds.py', 'finite_typed_adaptive_certificate.py',
             'check_finite_typed_adaptive_certificate.py', 'weighted_adaptive_mixed_place.py',
             'adaptive_mixed_place.py', 'check_weighted_adaptive_catalogues.py',
             'check_weighted_multipoint_hermite.py', 'check_ternary_basis_mixed.py',
             'check_multipoint_hermite_crossings.py', 'check_mixed_place_hermite.py')
    result = {'status': 'pass', 'scope': 'Exact rational finite rank upper bounds and actual typed reconstruction; no threshold-crossing block.',
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in names},
              'coefficient_upper_bounds': coefficient_checks, 'translation_invariance_checks': translated_checks,
              'positive_free_prefix_boundary': boundary, 'invalid_parameters_rejected': rejected,
              'complete_local_profiles': local_count, 'actual_piece_cut_ranks': piece_count,
              'endpoint_entries': endpoint_count, 'typed_cases': records, 'large_block_control': large_summary}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('typed_cases', 'source_sha256')}))


if __name__ == '__main__':
    main()
