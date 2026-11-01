#!/usr/bin/env python3
"""Complete exact F4 orbit and all-scale profile dominance for the q5 seed."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import accumulate
import json
from pathlib import Path
import sys

from joint_energy_conjugations import (IDENTITY, N, conjugated_values, field_product,
                                      inverse, model, orbit, unit_characters)
from weighted_smith_oracle_f4 import WeightedSmithOracle

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('joint_energy_conjugations.py', 'check_joint_energy_conjugations.py',
           'search_f4_recentered_crossings.py', 'weighted_smith_oracle_f4.py')


def polynomial(pair):
    a, b = pair
    return sum((((a >> k) & 1) | (((b >> k) & 1) << 1)) << (2*k)
               for k in range(max(a.bit_length(), b.bit_length())))


def cut_matrix(values, crossed):
    result = [[0]*N for _ in range(N)]
    for index, pair in enumerate(values):
        x, y, u, v = ((index//5**i) % 5 for i in range(4))
        row, column = (x+5*u, y+5*v) if crossed else (x+5*y, v+5*u)
        result[row][column] = polynomial(pair)
    return result


def hashes():
    return {name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest() for name in SOURCES}


def operator_coefficients(values):
    result = [[0, 0] for _ in range(19)]
    for index, pair in enumerate(values):
        x, y, u, v = ((index//5**i) % 5 for i in range(4))
        position = N*(x+5*y)+v+5*u
        for k in range(19):
            for b in range(2):
                result[k][b] |= ((pair[b] >> k) & 1) << position
    return list(map(tuple, result))


def weighted_coefficient(values, energies, power):
    result = [0, 0]
    for index, pair in enumerate(values):
        x, y, u, v = ((index//5**i) % 5 for i in range(4))
        row, column = x+5*y, v+5*u
        k = energies[column]-power
        if k >= 0:
            for b in range(2):
                result[b] |= ((pair[b] >> k) & 1) << (N*row+column)
    return tuple(result)


def run(saved=None):
    if saved is not None:
        assert saved['source_sha256'] == hashes()
    data = model()
    characters, full_characters = unit_characters(data)
    representatives, generators = orbit(data)
    # Simultaneous triangularization and independent diagonal characters
    # prove these unit counts over every finite field of characteristic 2.
    leading_units, stabilizer_units = 4**4*3**6, 4**2*3**2
    assert leading_units == 186624 and stabilizer_units == 144
    assert len(representatives) == leading_units//stabilizer_units == 1296
    source = data['values']
    identity_coordinates = data['leading'].coordinates((IDENTITY, 0))
    assert conjugated_values(data, identity_coordinates) == source
    source_coefficients = operator_coefficients(source)
    half = data['half']
    weights = [half[i//5] for i in range(N)]
    baseline, _ = WeightedSmithOracle(cut_matrix(source, True), infinity=True).profile(weights, weights)
    zero, _ = WeightedSmithOracle(cut_matrix(source, False)).profile(data['energies'], [0]*N)
    assert baseline == (-8, -7, -5, -4, -3, -2, -1, 0, 0, 1, 1, 2, 2,
                        3, 3, 4, 4, 4, 5, 5, 6, 7, 7, 8, 12)
    assert zero == tuple(v+16 for v in baseline)
    profiles, records = [], []
    for i, coordinates in enumerate(representatives):
        values = conjugated_values(data, coordinates)
        unit = data['leading'].matrix(coordinates)
        inv = inverse(unit, data['blocks'])
        assert inv is not None
        assert weighted_coefficient(values, data['energies'], -8) == (data['weighted'][-8], 0)
        assert weighted_coefficient(values, data['energies'], -7) == field_product(
            field_product(unit, (data['weighted'][-7], 0)), inv)
        for index, (a, b) in enumerate(values):
            x, y, u, v = ((index//5**j) % 5 for j in range(4))
            assert ((a.bit_count() % 2) | ((b.bit_count() % 2) << 1)) == int(x == v and y == u)
        # Fraction-free polynomial arithmetic; no truncated series or floats.
        orders, _ = WeightedSmithOracle(cut_matrix(values, True), infinity=True).profile(weights, weights)
        assert all(a <= b for a, b in zip(accumulate(orders), accumulate(baseline)))
        if orders not in profiles:
            # A different packed-matrix implementation checks every complete
            # coefficient at one representative of each exact profile class.
            assert operator_coefficients(values) == [field_product(field_product(unit, a), inv)
                                                     for a in source_coefficients]
            profiles.append(orders)
        record = [list(coordinates), profiles.index(orders)]
        records.append(record)
        if saved is not None:
            assert record == saved['records'][i]
            assert orders == tuple(saved['profiles'][record[1]])
        if (i+1) % 128 == 0:
            print(json.dumps({'completed_exact_orbit_profiles': i+1, 'total': 1296}), file=sys.stderr, flush=True)
    assert len(profiles) == 21
    evidence = dict(leading_algebra_dimension=10, full_commutant_dimension=4,
                    leading_diagonal_characters=characters, full_diagonal_characters=full_characters,
                    leading_units=leading_units, stabilizer_units=stabilizer_units,
                    bimodule_dimension=16, orbit_generators=generators,
                    baseline_infinity=baseline, invariant_zero=zero, profiles=profiles, records=records,
                    source_sha256=hashes(),
                    scope='All F4 constant joint-energy conjugations fixing the first weighted infinity coefficient of the q5 reference, at its stated natural physical heights.',
                    factoring_exponent='1/2; unchanged')
    if saved is not None:
        assert json.loads(json.dumps(evidence)) == saved
    report = dict(complete_orbit_size=1296, distinct_exact_profiles=21,
                  fraction_free_profiles_replayed=1298, leading_algebra_dimension=10,
                  full_commutant_dimension=4, compressed_orbit_module_dimension=16,
                  exact_unit_counts=[leading_units, stabilizer_units],
                  conclusion='Every candidate is weakly prefix-dominated by the reference, so cannot improve its cost at any real scale >=1 at the stated heights.',
                  factoring_exponent='1/2; unchanged')
    return evidence, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/joint-energy-conjugations.json')
    args = parser.parse_args()
    saved = None if args.generate else json.loads(args.evidence.read_text())
    evidence, report = run(saved)
    if args.generate:
        args.evidence.write_text(json.dumps(evidence, separators=(',', ':'))+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
