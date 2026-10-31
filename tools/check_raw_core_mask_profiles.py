#!/usr/bin/env python3
"""Literal endpoint and both-cut checks for arbitrary binary/F4 core masks."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import numpy as np

from check_rectangular_f4_repair import fmul
from native_matrix_smith_f4 import MatrixOracle
from raw_core_mask_profiles import CoreMaskPressure, index, normalize_mask, pair_partition_mask
from relabelled_gaussian_sums import gaussian_tensor, flatten, natural_heights


def masked_tensor(interval, mask):
    base = gaussian_tensor(tuple(range(interval)), interval, False)
    full = 1 << interval
    local = lambda state: (state & 1) | (((state >> (interval - 1)) & 1) << 1)
    result = []
    for packed, polynomial in enumerate(base):
        labels = tuple(local((packed >> (leg * interval)) & (full - 1)) for leg in range(4))
        scalar = mask[index(labels)]
        result.append((polynomial if scalar & 1 else 0, polynomial if scalar & 2 else 0))
    return result


def endpoint_controls(masks):
    interval, q = 2, 4
    words = list(product(range(q), repeat=2))
    totals = {word: tuple(sum((state >> i) & 1 for state in word) for i in range(interval))
              for word in words}
    pairs = [(x, v) for x in words for v in words if totals[x] == totals[v]]
    checked = 0
    for mask in masks:
        endpoint = [(lo.bit_count() & 1) + 2 * (hi.bit_count() & 1)
                    for lo, hi in masked_tensor(interval, mask)]
        for x, v in pairs:
            for y, u in pairs:
                value = 1
                for copy in range(2):
                    packed = sum(labels[copy] << (interval * leg)
                                 for leg, labels in enumerate((x, y, u, v)))
                    value = fmul(value, endpoint[packed])
                assert value == int(x == v and y == u)
                checked += 1
    return checked


def direct_profiles(interval, mask, offsets, bulk):
    q = 1 << interval
    local = lambda state: (state & 1) | (((state >> (interval - 1)) & 1) << 1)
    values = masked_tensor(interval, mask)
    heights = natural_heights(tuple(range(interval)), interval).copy()
    for place in range(2):
        for leg in range(4):
            heights[place, leg] += [offsets[place, leg, local(state)] for state in range(q)]
    for state in range(q):
        count = state.bit_count()
        heights[0, 0, state] += bulk[0] * count
        heights[0, 1, state] -= bulk[0] * count
        heights[1, 1, state] += bulk[1] * count
        heights[1, 3, state] += bulk[1] * count
    answer = []
    for place in range(2):
        first, second = ((0, 1), (2, 3)) if place == 0 else ((0, 2), (1, 3))
        rows = [heights[place, first[0], a] + heights[place, first[1], b]
                for a, b in product(range(q), repeat=2)]
        columns = [heights[place, second[0], a] + heights[place, second[1], b]
                   for a, b in product(range(q), repeat=2)]
        answer.append(MatrixOracle(flatten(values, interval, bool(place)), q * q).profile(
            rows, columns, infinity=bool(place))['orders'])
    return answer


def verify():
    rng = np.random.default_rng(528713)
    random_mask = rng.integers(0, 4, 256).tolist()
    for x, y in product(range(4), repeat=2):
        random_mask[index((x, y, y, x))] = 1
    masks = [pair_partition_mask(range(10)), pair_partition_mask((0, 0, 1, 1, 0, 1, 0, 0, 1, 0)),
             normalize_mask(random_mask)]
    endpoints = endpoint_controls(masks)
    records = []
    for interval in (2, 4):
        for kind, mask in enumerate(masks):
            model = CoreMaskPressure(interval, mask)
            for trial in range(2):
                offsets = rng.integers(-8, 9, (2, 4, 4)) / 4
                offsets[:, :, 0] = 0
                bulk = rng.integers(-4, 5, 2) / 4
                expected = direct_profiles(interval, mask, offsets, bulk)
                _, detail = model.evaluate(1.0, offsets, bulk)
                for side in range(2):
                    actual = np.asarray(detail['active_orders'][side])
                    data = model.records[side][0]
                    for order, charge in zip(data['spectator_orders'], data['spectator_charges']):
                        actual = np.r_[actual, actual + order + bulk[side] * charge]
                    assert np.array_equal(np.sort(actual), np.sort(expected[side]))
                    records.append({'interval': interval, 'mask': kind, 'trial': trial,
                                    'place': side, 'rank': len(actual),
                                    'largest_active_charge_block': max(model.oracles[side].block_sizes)})
    root = Path(__file__).resolve().parent.parent
    files = ['tools/blocked_core_charge_pressure.py', 'tools/raw_core_mask_profiles.py',
             'tools/check_raw_core_mask_profiles.py', 'research/raw-core-transition-masks.md']
    return {'status': 'pass', 'typed_endpoint_entries': endpoints,
            'original_weighted_cut_profiles': records,
            'source_sha256': {name: sha256((root / name).read_bytes()).hexdigest() for name in files}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'typed_endpoint_entries': result['typed_endpoint_entries'],
                      'original_weighted_cut_profiles': len(result['original_weighted_cut_profiles'])}))


if __name__ == '__main__':
    main()
