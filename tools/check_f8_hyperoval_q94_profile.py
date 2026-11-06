#!/usr/bin/env python3
"""Verify the recorded complete q94 F8 hyperoval profile."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps
from search_coherent_hot_retractions import pressure

INPUT = Path('evidence/f8-hyperoval-q94-profile.json')
BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')
SOURCE_NAMES = (
    'coherent_f8_maps.py',
    'coherent_hot_compounds.py',
    'coherent_hot_compounds_f8.py',
    'coherent_hot_compounds_f8.cpp',
    'selected_gaussian_coefficients.cpp',
    'native_rectangular_smith_f4.cpp',
    'search_coherent_hot_retractions.py',
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved['completed'] is True
    assert saved['q'] == 94 and saved['modes'] == 11 and saved['cutoff'] == 135
    assert saved['base_profile_sha256'] == sha256(BASE.read_bytes()).hexdigest()
    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]

    def uses_removed_mode(table):
        return any(word >> 11 & 1 for word in table)

    keep = [index for index, (left, right) in enumerate(zip(embedding, retraction))
            if not uses_removed_mode(left) and not uses_removed_mode(right)]
    assert keep == saved['kept_source_indices']
    embedding = [embedding[index] for index in keep]
    retraction = [retraction[index] for index in keep]
    heights, charges = validate_maps(embedding, retraction, 11)
    assert heights == saved['heights']
    assert charges == saved['charges']
    assert max(heights) == 25 and sum(heights) == 1656
    assert saved['all_orders_resolved'] is True
    assert saved['tail_count'] == [0, 0]
    for orders in (saved['zero_orders_t'], saved['infinity_orders_t']):
        assert len(orders) == 8836
        assert orders == sorted(orders)
        assert min(orders) == 0 and max(orders) <= saved['cutoff']

    profile = {
        'q': saved['q'],
        'modes': saved['modes'],
        'heights': saved['heights'],
        'charges': saved['charges'],
        'zero_orders_t': saved['zero_orders_t'],
        'infinity_orders_t': saved['infinity_orders_t'],
    }
    fresh = pressure(profile)
    assert fresh['completed'] is True
    assert abs(fresh['rho'] - saved['numerical_pressure']['rho']) < 1e-12
    assert abs(fresh['scale_t'] - saved['numerical_pressure']['scale_t']) < 1e-12
    assert fresh['rho'] > 2 ** -0.5
    assert fresh['rho'] > 29 / 41
    for name in SOURCE_NAMES:
        assert sha256((Path('tools') / name).read_bytes()).hexdigest() == saved['source_sha256'][name]

    receipt = {
        'passed': True,
        'q': 94,
        'max_height': 25,
        'orders_per_cut': 8836,
        'tail_count': [0, 0],
        'LP_identity_and_deletion_closure': True,
        'numerical_pressure': fresh,
        'above_inverse_sqrt_two': True,
        'above_height38_floor_value': True,
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'base_profile_sha256': saved['base_profile_sha256'],
        'scope': (
            'Replay of the q94 mode-11-removed map and its recorded complete '
            'finite orders. The map does not retain the full height-below-24 core, so '
            'the floor theorem does not apply; the pressure is numerical and the '
            'complete exponent remains 1/2.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'q': 94, 'rho': fresh['rho'],
                      'above': True}))


if __name__ == '__main__':
    main()
