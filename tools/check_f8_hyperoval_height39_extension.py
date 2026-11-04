#!/usr/bin/env python3
"""Verify the recorded q99 height-39 F8 hyperoval profile and its scope."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import sys

from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps
from search_coherent_hot_retractions import pressure

INPUT = Path('evidence/f8-hyperoval-height39-extension.json')
BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')
EXTENSION_WORDS = (160, 192, 224)
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
    assert saved['extension_words'] == list(EXTENSION_WORDS)
    assert saved['q'] == 99 and saved['modes'] == 12
    assert saved['base_profile_sha256'] == sha256(BASE.read_bytes()).hexdigest()

    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    for word in EXTENSION_WORDS:
        embedding.append({word: 1})
        retraction.append({word: 1})
    heights, charges = validate_maps(embedding, retraction, 12)
    assert heights == saved['heights']
    assert charges == saved['charges']
    assert heights[-3:] == [26, 28, 39]
    assert saved['all_orders_resolved'] is True
    assert saved['tail_count'] == [0, 0]
    assert len(saved['zero_orders_t']) == 9801
    assert len(saved['infinity_orders_t']) == 9801
    for orders in (saved['zero_orders_t'], saved['infinity_orders_t']):
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
    # A flat numerical minimum need not return identical argmin bits across
    # hosts. Check the pressure at the recorded scale, rather than its identity
    # with a new optimizer's scale.
    scale = saved['numerical_pressure']['scale_t']
    rho = saved['numerical_pressure']['rho']
    assert isfinite(scale) and .01 <= scale <= 2
    target = sum(float(logsumexp([-scale * order for order in orders]))
                 for orders in (profile['zero_orders_t'], profile['infinity_orders_t']))
    height_term = float(logsumexp([-scale * height / rho for height in heights]))
    saved_scale_residual = abs(target - 4 * rho * height_term)
    assert saved_scale_residual < 1e-11
    assert fresh['rho'] > 2 ** -0.5

    for name in SOURCE_NAMES:
        assert sha256((Path('tools') / name).read_bytes()).hexdigest() == saved['source_sha256'][name]

    receipt = {
        'passed': True,
        'q': 99,
        'max_height': 39,
        'orders_per_cut': 9801,
        'tail_count': [0, 0],
        'LP_identity_and_deletion_closure': True,
        'numerical_pressure': fresh,
        'saved_scale_equation_residual': saved_scale_residual,
        'above_inverse_sqrt_two': True,
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'base_profile_sha256': saved['base_profile_sha256'],
        'scope': (
            'Replay of the map identity and numerical pressure from the recorded '
            'complete q99 orders. The full profile generation is not replayed here; '
            'the pressure is numerical and the complete exponent remains 1/2.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'q': 99, 'rho': fresh['rho'], 'above': True}))


if __name__ == '__main__':
    main()
