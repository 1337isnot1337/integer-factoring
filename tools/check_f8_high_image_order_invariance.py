#!/usr/bin/env python3
"""Replay the q96 cutoff-47 boundary-image order and pressure control."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps
from coherent_hot_compounds import height
from probe_determinant_inverse_codes import fit

INPUT = Path('evidence/f8-high-image-order-invariance.json')
SOURCES = [
    'probe_f8_hyperoval_retraction.py',
    'coherent_f8_maps.py',
    'coherent_hot_compounds_f8.py',
    'coherent_hot_compounds_f8.cpp',
    'search_coherent_hot_retractions.py',
    'probe_determinant_inverse_codes.py',
]
PIVOTS = [7, 8, 9]


def build_maps(A):
    boundary = sorted(word for word in range(4096)
                      if height(word) == 25 and word.bit_count() == 3)
    low = [word for word in range(4096) if height(word) < 25] + [31]
    nonpivots = [column for column in range(10) if column not in PIVOTS]
    rows = [[0] * 10 for _ in range(3)]
    for index, pivot in enumerate(PIVOTS):
        rows[index][pivot] = 1
    for k, column in enumerate(nonpivots):
        for index in range(3):
            rows[index][column] = A[index][k]
    embedding = [{word: 1} for word in low]
    embedding += [{boundary[column]: rows[index][column]
                   for column in range(10) if rows[index][column]}
                  for index in range(3)]
    retraction = [{word: 1} for word in low]
    retraction += [{boundary[column]: 1} for column in PIVOTS]
    heights, charges = validate_maps(embedding, retraction, 12)
    return embedding, retraction, heights, charges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved['completed'] is True
    assert saved['q'] == 96 and saved['modes'] == 12 and saved['cutoff'] == 47
    assert saved['pivot_columns'] == PIVOTS
    assert saved['source_sha256'] == {
        name: sha256((Path(__file__).resolve().parents[1] / 'tools' / name).read_bytes()).hexdigest()
        for name in SOURCES
    }

    base_orders = None
    result = []
    for image in saved['images']:
        name = image['name']
        embedding, retraction, heights, charges = build_maps(image['A'])
        assert len(embedding) == 96 and len(heights) == 96
        profile = next(row for row in saved['profiles'] if row['name'] == name)
        orders = saved['orders'][name]
        assert len(orders) == 2
        known = fit({'cuts': [{'retained_orders_t': cut} for cut in orders]}, heights)
        augmented = [
            sorted(cut + [saved['cutoff'] + 1] * profile['tails'][index])
            for index, cut in enumerate(orders)
        ]
        upper = fit({'cuts': [{'retained_orders_t': cut} for cut in augmented]}, heights)
        hashes = [sha256(json.dumps(cut).encode()).hexdigest() for cut in orders]
        assert hashes == profile['order_sha256']
        assert abs(known['rho'] - profile['known_lower_rho']) < 1e-14
        assert abs(upper['rho'] - profile['aug_upper_rho']) < 1e-14
        assert known['rho'] > 1 / math.sqrt(2)
        assert known['rho'] <= upper['rho']
        assert hash(tuple(orders[0])) == hash(tuple(orders[1]))
        if base_orders is None:
            base_orders = orders
        else:
            assert orders == base_orders
        result.append({
            'name': name,
            'retained': [len(cut) for cut in orders],
            'tails': profile['tails'],
            'order_sha256': hashes,
            'known_lower_rho': known['rho'],
            'tail_augmented_upper_rho': upper['rho'],
        })

    assert all(row['retained'] == [935, 935] for row in result)
    assert all(row['tails'] == [8281, 8281] for row in result)
    assert all(row['order_sha256'][0] == row['order_sha256'][1] for row in result)
    assert result[0]['order_sha256'] == result[1]['order_sha256']
    assert result[0]['known_lower_rho'] > 1 / math.sqrt(2) + 0.0027

    receipt = {
        'status': 'PASS',
        'q': 96,
        'cutoff': 47,
        'retained_orders_per_cut': 935,
        'unresolved_tails_per_cut': 8281,
        'known_order_lower_rho': result[0]['known_lower_rho'],
        'tail_augmented_upper_rho': result[0]['tail_augmented_upper_rho'],
        'order_hashes_identical': True,
        'common_core_floor': '17/24',
        'complete_profile_recomputed_this_run': False,
        'factoring_exponent_improved': False,
        'scope': saved['scope'],
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
