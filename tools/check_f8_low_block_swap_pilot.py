#!/usr/bin/env python3
"""Verify the recorded cutoff-47 F8 low-block swap pilot."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps
from search_coherent_hot_retractions import pressure

INPUT = Path('evidence/f8-low-block-swap-pilot.json')
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
    assert saved['q'] == 96 and saved['modes'] == 12 and saved['cutoff'] == 47
    assert saved['dropped_word'] == 31 and saved['added_word'] == 46
    assert saved['base_profile_sha256'] == sha256(BASE.read_bytes()).hexdigest()
    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    keep = [index for index, (left, right) in enumerate(zip(embedding, retraction))
            if not (len(left) == 1 and len(right) == 1 and next(iter(left)) == 31)]
    assert keep == saved['kept_source_indices']
    embedding = [embedding[index] for index in keep] + [{46: 1}]
    retraction = [retraction[index] for index in keep] + [{46: 1}]
    heights, charges = validate_maps(embedding, retraction, 12)
    assert heights == saved['heights']
    assert charges == saved['charges']
    assert max(heights) == 26 and sum(heights) == 1704
    assert saved['all_orders_resolved'] is False
    assert saved['tail_count'] == [8281, 8281]
    for orders in (saved['zero_orders_t'], saved['infinity_orders_t']):
        assert len(orders) == 935
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
    assert fresh['rho'] > 0.7098643133625308
    for name in SOURCE_NAMES:
        assert sha256((Path('tools') / name).read_bytes()).hexdigest() == saved['source_sha256'][name]

    receipt = {
        'passed': True,
        'q': 96,
        'max_height': 26,
        'orders_per_cut': 9216,
        'resolved_orders_per_cut': 935,
        'tail_count': [8281, 8281],
        'numerical_pressure': fresh,
        'truncated_pressure_is_lower_bound': True,
        'above_inverse_sqrt_two': True,
        'dominated_by_q97_add46_pilot': True,
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'base_profile_sha256': saved['base_profile_sha256'],
        'scope': (
            'Replay of the cutoff-47 swap pilot dropping word 31 and adding '
            'word 46. The profile is truncated; the pressure is a lower bound, '
            'and it exceeds 1/sqrt(2), so the complete swap map cannot cross. '
            'No complete profile or exponent change is claimed.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'q': 96, 'rho': fresh['rho'],
                      'truncated': True}))


if __name__ == '__main__':
    main()
