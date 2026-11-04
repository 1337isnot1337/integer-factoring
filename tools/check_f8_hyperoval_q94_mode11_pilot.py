#!/usr/bin/env python3
"""Replay the q94 mode-11 pilot map and its recorded cutoff-47 pressure."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import isfinite, sqrt
from pathlib import Path
import sys

from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps
from search_coherent_hot_retractions import pressure

INPUT = Path('evidence/f8-hyperoval-q94-mode11-pilot.json')
BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')
DROP = [90, 91]
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
    assert saved['q'] == 94 and saved['modes'] == 11
    assert saved['drop_indices'] == DROP
    assert saved['cutoff'] == 47
    assert saved['base_profile_sha256'] == sha256(BASE.read_bytes()).hexdigest()

    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    bad = [i for i, vector in enumerate(embedding)
           if any(word >> 11 & 1 for word in vector)]
    assert bad == DROP
    keep = [i for i in range(len(embedding)) if i not in bad]
    embedding = [embedding[i] for i in keep]
    retraction = [retraction[i] for i in keep]
    heights, charges = validate_maps(embedding, retraction, 11)
    assert heights == saved['heights']
    assert charges == saved['charges']
    assert max(heights) == 25
    assert len(embedding) == 94

    orders = (saved['zero_orders_t'], saved['infinity_orders_t'])
    assert all(orders_cut == sorted(orders_cut) for orders_cut in orders)
    assert [len(orders_cut) for orders_cut in orders] == [932, 932]
    assert saved['tail_counts'] == [7904, 7904]
    assert saved['all_orders_resolved'] is False

    fresh = pressure({
        'q': saved['q'],
        'modes': saved['modes'],
        'heights': heights,
        'charges': charges,
        'zero_orders_t': orders[0],
        'infinity_orders_t': orders[1],
    })
    assert fresh['completed'] is True
    saved_fit = saved['numerical_pressure']
    assert abs(fresh['rho'] - saved_fit['rho']) < 1e-12
    scale = saved_fit['scale_t']
    rho = saved_fit['rho']
    assert isfinite(scale) and 0.01 <= scale <= 2
    target = sum(float(logsumexp([-scale * order for order in orders_cut]))
                 for orders_cut in orders)
    height_term = float(logsumexp([-scale * height / rho for height in heights]))
    residual = abs(target - 4 * rho * height_term)
    assert residual < 1e-11
    assert rho > 1 / sqrt(2)
    assert rho < 0.7099930431440509

    for name in SOURCE_NAMES:
        assert sha256((Path('tools') / name).read_bytes()).hexdigest() == saved['source_sha256'][name]

    receipt = {
        'passed': True,
        'q': 94,
        'modes': 11,
        'drop_indices': DROP,
        'orders_per_cut': 932,
        'tail_counts': saved['tail_counts'],
        'LP_identity_and_deletion_closure': True,
        'numerical_pressure': fresh,
        'saved_scale_equation_residual': residual,
        'above_inverse_sqrt_two': True,
        'below_q99_pilot': True,
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'scope': (
            'Replay of the q94 mode-11 drop and the recorded cutoff-47 pilot '
            'pressure. The profile is incomplete: 7,904 orders per cut are '
            'unresolved. The pilot value is not a complete profile or an exact '
            'certificate, and the complete exponent remains 1/2.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'q': 94, 'rho': fresh['rho'],
                      'tail': saved['tail_counts'], 'above': True}))


if __name__ == '__main__':
    main()
