#!/usr/bin/env python3
"""Replay the cutoff-47 single-mode drop census of the q96 F8 map."""
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

INPUT = Path('evidence/f8-single-mode-drop-census.json')
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
    assert saved['base_profile_sha256'] == sha256(BASE.read_bytes()).hexdigest()
    assert saved['source_sha256'] == {
        name: sha256((Path(__file__).resolve().parents[1] / 'tools' / name).read_bytes()).hexdigest()
        for name in SOURCES
    }

    embedding, retraction, _, _ = maps(saved['pivots'])
    rows = []
    for image in saved['images']:
        mode = image['mode']
        keep = [index for index, vector in enumerate(embedding)
                if not any(word >> mode & 1 for word in set(vector) | set(retraction[index]))]
        dropped = [index for index in range(len(embedding)) if index not in keep]
        kept_embedding = [embedding[index] for index in keep]
        kept_retraction = [retraction[index] for index in keep]
        heights, charges = validate_maps(kept_embedding, kept_retraction, 12)
        assert len(kept_embedding) == image['q'] == 96 - image['dropped']
        orders = saved['orders'][str(mode)]
        assert len(orders) == 2
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
        if mode == 11:
            assert sorted(dropped) == [90, 91] and image['q'] == 94
        rows.append({
            'mode': mode,
            'q': image['q'],
            'dropped_indices': dropped,
            'order_sha256': hashes,
            'known_lower_rho': known['rho'],
            'tail_augmented_upper_rho': upper['rho'],
        })

    assert len(rows) == 12
    assert sorted(row['mode'] for row in rows) == list(range(12))
    best = min(rows, key=lambda row: row['known_lower_rho'])
    assert best['mode'] == 11
    assert best['known_lower_rho'] > 1 / math.sqrt(2) + 0.0028
    receipt = {
        'status': 'PASS',
        'q': 96,
        'cutoff': 47,
        'valid_drops': len(rows),
        'best_mode': best['mode'],
        'best_q': best['q'],
        'best_known_order_lower_rho': best['known_lower_rho'],
        'best_tail_augmented_upper_rho': best['tail_augmented_upper_rho'],
        'all_lower_bounds_above_critical': True,
        'complete_profile_recomputed_this_run': False,
        'factoring_exponent_improved': False,
        'scope': saved['scope'],
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
