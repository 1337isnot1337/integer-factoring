#!/usr/bin/env python3
"""Replay the cutoff-47 two-mode drop screen containing mode 11."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps
from probe_determinant_inverse_codes import fit
from probe_f8_hyperoval_retraction import maps

INPUT = Path('evidence/f8-two-mode-drop-screen.json')
BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')
SOURCES = [
    'probe_f8_hyperoval_retraction.py',
    'coherent_f8_maps.py',
    'coherent_hot_compounds_f8.py',
    'coherent_hot_compounds_f8.cpp',
    'search_coherent_hot_retractions.py',
    'probe_determinant_inverse_codes.py',
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved['completed'] is True
    assert saved['q'] == 96 and saved['cutoff'] == 47 and saved['pivots'] == [7, 8, 9]
    assert saved['pairs'] == [[11, mode] for mode in range(11)]
    assert saved['base_profile_sha256'] == sha256(BASE.read_bytes()).hexdigest()
    assert saved['source_sha256'] == {
        name: sha256((Path(__file__).resolve().parents[1] / 'tools' / name).read_bytes()).hexdigest()
        for name in SOURCES
    }

    embedding, retraction, _, _ = maps(saved['pivots'])
    rows = []
    for image in saved['images']:
        pair = image['pair']
        keep = [index for index, vector in enumerate(embedding)
                if not any(((word >> pair[0]) | (word >> pair[1])) & 1
                           for word in set(vector) | set(retraction[index]))]
        assert len(keep) == image['q'] == 96 - image['dropped']
        heights, _ = validate_maps(
            [embedding[index] for index in keep],
            [retraction[index] for index in keep],
            12,
        )
        orders = saved['orders'][str(tuple(pair))]
        known = fit({'cuts': [{'retained_orders_t': cut} for cut in orders]}, heights)
        augmented = [
            sorted(cut + [saved['cutoff'] + 1] * image['tails'][index])
            for index, cut in enumerate(orders)
        ]
        upper = fit({'cuts': [{'retained_orders_t': cut} for cut in augmented]}, heights)
        hashes = [sha256(json.dumps(cut).encode()).hexdigest() for cut in orders]
        assert hashes == image['order_sha256']
        assert abs(known['rho'] - image['known_lower_rho']) < 1e-14
        assert abs(upper['rho'] - image['aug_upper_rho']) < 1e-14
        assert known['rho'] > 1 / math.sqrt(2)
        assert known['rho'] <= upper['rho']
        rows.append({'pair': pair, 'q': image['q'], 'known_lower_rho': known['rho'],
                     'tail_augmented_upper_rho': upper['rho']})

    assert len(rows) == 11
    assert {tuple(row['pair']) for row in rows} == {(11, mode) for mode in range(11)}
    best = min(rows, key=lambda row: row['known_lower_rho'])
    assert best['pair'] == [11, 10]
    assert best['known_lower_rho'] > 0.71005
    receipt = {
        'status': 'PASS',
        'q': 96,
        'cutoff': 47,
        'valid_pairs': len(rows),
        'best_pair': best['pair'],
        'best_q': best['q'],
        'best_known_order_lower_rho': best['known_lower_rho'],
        'best_tail_augmented_upper_rho': best['tail_augmented_upper_rho'],
        'all_lower_bounds_above_critical': True,
        'completed_q94_pressure': 0.7100458191653984,
        'best_pair_above_completed_q94': True,
        'factoring_exponent_improved': False,
        'scope': saved['scope'],
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
