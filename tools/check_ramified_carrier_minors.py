#!/usr/bin/env python3
"""Exact controls for ramification-residue minor domination and its Bell use."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from math import comb, inf
from pathlib import Path
import random

from affine_place_profile_oracle import AffinePlaceOracle
from bell_jet_crt_crossings import construct
from check_fused_multipoint_profiles import determinant, field_rank
from f16_identity_operator_profiles import MUL, add, multiply
from physical_split_crossings import identity_split
from reflection_laurent_physical_frames import digest
from search_entrywise_reflection_jets import bell_values


def matrix_product(first, second):
    result = [[{} for _ in second[0]] for _ in first]
    for i in range(len(first)):
        for k in range(len(second)):
            for j in range(len(second[0])):
                result[i][j] = add(result[i][j], multiply(first[i][k], second[k][j]))
    return result


def integer_certificate():
    numerator = {0: 6, 9: -4, 81: -6, 151: 9, 153: -4, 225: -1}
    quotient = [sum((k - i + 1) * value for i, value in numerator.items() if i <= k)
                for k in range(224)]
    replay = [0] * 226
    for k, value in enumerate(quotient):
        replay[k] += value
        replay[k + 1] -= 2 * value
        replay[k + 2] += value
    assert replay == [numerator.get(i, 0) for i in range(226)]
    transformed = [sum(quotient[k] * comb(223 - k, j - k) for k in range(j + 1))
                   for j in range(224)]
    assert min(transformed) == 6 and all(value > 0 for value in transformed)
    return {'numerator_after_t63': numerator, 'quotient_degree': 223,
            'quotient_coefficients': quotient, 'positive_transform_degree': 223,
            'positive_transform_coefficients': transformed,
            'minimum_positive_coefficient': min(transformed)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(312517)
    records, minor_controls = [], 0
    for trial in range(12):
        size = 3 + trial % 2
        ramification = 2 + trial % 3
        rank_cap = min(ramification - 1, size)
        diagonal = [[{i: 1} if i == j else {} for j in range(size)] for i in range(size)]
        first, second = [], []
        for i in range(size):
            first.append([{0: 1} if i == j else ({1: rng.randrange(1, 4)} if i < j else {})
                          for j in range(size)])
            second.append([{0: 1} if i == j else ({1: rng.randrange(1, 4)} if i > j else {})
                           for j in range(size)])
        unramified = matrix_product(matrix_product(first, diagonal), second)
        source = [[{ramification * exponent: coefficient for exponent, coefficient in value.items()}
                   for value in row] for row in unramified]
        left = [[rng.randrange(4) for _ in range(rank_cap)] for _ in range(size)]
        right = [[rng.randrange(4) for _ in range(size)] for _ in range(rank_cap)]
        update = [[0] * size for _ in range(size)]
        for i in range(size):
            for j in range(size):
                for k in range(rank_cap):
                    update[i][j] ^= MUL[left[i][k]][right[k][j]]
        assert field_rank(update) <= rank_cap
        order = ramification - 1
        eta = {order: 1, order + 1: 2, order + 3: 3}
        changed = [[add(source[i][j], eta, update[i][j]) for j in range(size)] for i in range(size)]
        row_heights = [rng.randrange(-3, 4) for _ in range(size)]
        column_heights = [rng.randrange(-3, 4) for _ in range(size)]
        prefixes = [[], []]
        for width in range(1, size + 1):
            best = [inf, inf]
            for rows in combinations(range(size), width):
                for columns in combinations(range(size), width):
                    old = determinant([[source[i][j] for j in columns] for i in rows])
                    new = determinant([[changed[i][j] for j in columns] for i in rows])
                    old_order, new_order = min(old, default=inf), min(new, default=inf)
                    assert new_order <= old_order
                    shift = sum(row_heights[i] for i in rows) + sum(column_heights[j] for j in columns)
                    best[0] = min(best[0], old_order + shift)
                    best[1] = min(best[1], new_order + shift)
                    minor_controls += 1
            assert best[1] <= best[0] and all(value != inf for value in best)
            for side in (0, 1):
                prefixes[side].append(best[side])
        records.append({'size': size, 'ramification': ramification, 'eta_order': order,
                        'update_rank': field_rank(update), 'row_heights': row_heights,
                        'column_heights': column_heights,
                        'baseline_weighted_prefixes': prefixes[0], 'updated_weighted_prefixes': prefixes[1]})
    counterexample = [[{0: 1}, {1: 1, 2: 1}], [{1: 1, 2: 1}, {2: 1}]]
    assert determinant(counterexample) == {4: 1}
    identity = identity_split(4)
    alternating = [[int(j == (i ^ 1)) for j in range(4)] for i in range(4)]
    row_bell, column_bell = bell_values(alternating, 0), bell_values(identity, 1)
    carrier = []
    for entry in range(256):
        x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
        carrier.append(int(x == v and y == u) ^ MUL[2][row_bell[entry]] ^ MUL[3][column_bell[entry]])
    baseline = construct(alternating, identity, carrier)
    base_oracle = AffinePlaceOracle(baseline)
    expected = (0,) + (8,) * 4 + (16,) * 6 + (24,) * 4 + (32,)
    for cut, point in ((0, (2, 1)), (1, (3, 1))):
        assert base_oracle.profile(cut, point, [0] * 16) == expected
    available = []
    for entry in range(256):
        x, y, u, v = ((entry // 4 ** leg) % 4 for leg in range(4))
        if x != v and y != u:
            available.append(entry)
    weighted = []
    heights_list = [[0] * 16, [rng.randrange(-3, 4) for _ in range(16)]]
    for rank_cap in (1, 3, 7):
        selected = rng.sample(available, rank_cap)
        changed_carrier = list(carrier)
        for entry in selected:
            changed_carrier[entry] ^= rng.randrange(1, 4)
        changed = construct(alternating, identity, changed_carrier)
        oracle = AffinePlaceOracle(changed)
        for heights in heights_list:
            for cut, point in ((0, (2, 1)), (1, (3, 1))):
                old = tuple(round(value) for value in base_oracle.log_profile_with_gradient(cut, point, heights)[3])
                new = tuple(round(value) for value in oracle.log_profile_with_gradient(cut, point, heights)[3])
                assert all(sum(new[:k]) <= sum(old[:k]) for k in range(1, 17))
                weighted.append({'update_rank_bound': rank_cap, 'cut': cut, 'point': list(point),
                                 'heights': heights, 'baseline_orders': old, 'updated_orders': new,
                                 'tensor_sha256': digest(changed['values'])})
    names = (Path(__file__).name, 'bell_jet_crt_crossings.py', 'affine_place_profile_oracle.py',
             'check_fused_multipoint_profiles.py', 'f16_weighted_local_oracle.py',
             'native_matrix_smith_f4.py', 'native_matrix_smith_f4.cpp')
    report = {'status': 'passed',
              'scope': 'Exact residue-minor controls and integer unweighted domination certificate at the two specified Bell places; no exclusion for extra places or full-rank physical weighting.',
              'small_matrix_cases': records, 'all_minor_controls': minor_controls,
              'rank_threshold_counterexample_orders': [2, 4],
              'alternating_baseline_complete_orders': expected,
              'weighted_literal_profile_controls': weighted,
              'integer_positivity_certificate': integer_certificate(),
              'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                for name in names}}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'all_minor_controls': minor_controls,
                      'literal_weighted_profiles': len(weighted), 'positive_integer_coefficients': 224,
                      'minimum_positive_coefficient': 6}))


if __name__ == '__main__':
    main()
