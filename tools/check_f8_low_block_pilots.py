#!/usr/bin/env python3
"""Reproduce the F8 low-block cutoff-47 pilots and their crossing exclusion."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps
from search_coherent_hot_retractions import pressure

INPUT = Path('evidence/f8-low-block-pilots.json')
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


def build(name):
    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    if name == 'q95_drop31':
        keep = [i for i, row in enumerate(embedding)
                if not (len(row) == 1 and len(retraction[i]) == 1
                        and next(iter(row)) == 31)]
        modes = 12
    elif name == 'q97_add46':
        keep = list(range(len(embedding)))
        modes = 12
    elif name == 'q97_q94q99':
        keep = [i for i, row in enumerate(embedding)
                if all(not (int(k) >> 11 & 1) for k in set(row) | set(retraction[i]))]
        modes = 11
    else:
        raise ValueError(name)
    embedding = [embedding[i] for i in keep]
    retraction = [retraction[i] for i in keep]
    if name == 'q97_add46':
        embedding.append({46: 1})
        retraction.append({46: 1})
    elif name == 'q97_q94q99':
        for word in (160, 192, 224):
            embedding.append({word: 1})
            retraction.append({word: 1})
    return embedding, retraction, modes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved['completed'] is True
    base = json.loads(BASE.read_text())
    assert sha256(BASE.read_bytes()).hexdigest() == saved['base_profile_sha256']
    results = []
    for pilot in saved['pilots']:
        embedding, retraction, modes = build(pilot['name'])
        heights, charges = validate_maps(embedding, retraction, modes)
        assert heights == pilot['heights'] and charges == pilot['charges']
        assert min(heights) == 0
        assert pilot['q'] == len(heights) and pilot['modes'] == modes
        assert pilot['cutoff'] == 47
        zero = pilot['zero_orders_t']
        infinity = pilot['infinity_orders_t']
        for orders in (zero, infinity):
            assert orders == sorted(orders)
            assert min(orders) >= 0 and max(orders) <= pilot['cutoff']
        assert pilot['tail_counts'] == [pilot['q'] ** 2 - len(zero),
                                        pilot['q'] ** 2 - len(infinity)]
        profile = {
            'q': pilot['q'], 'modes': modes, 'heights': heights, 'charges': charges,
            'zero_orders_t': zero, 'infinity_orders_t': infinity,
        }
        fresh = pressure(profile)
        assert fresh['completed'] is True
        assert abs(fresh['rho'] - pilot['numerical_pressure']['rho']) < 1e-12
        assert abs(fresh['scale_t'] - pilot['numerical_pressure']['scale_t']) < 1e-12
        assert fresh['rho'] > 2 ** -0.5
        results.append({
            'name': pilot['name'], 'q': pilot['q'], 'modes': modes,
            'pilot_orders_per_cut': len(zero), 'tail_counts': pilot['tail_counts'],
            'pilot_rho': fresh['rho'], 'pilot_scale_t': fresh['scale_t'],
            'above_inverse_sqrt_two': True, 'complete_cannot_cross': True,
        })
    receipt = {
        'passed': True,
        'base_profile_sha256': sha256(BASE.read_bytes()).hexdigest(),
        'pilots': results,
        'lemma': 'truncated-profile-lower-bound.md: a resolved subset gives a lower bound on the complete pressure; a pilot above 1/sqrt(2) excludes the complete map from crossing.',
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'source_sha256': {name: sha256((Path('tools') / name).read_bytes()).hexdigest()
                          for name in SOURCE_NAMES},
        'scope': 'Deterministic replay of the recorded cutoff-47 orders and map constructions. No complete profile, certificate, or score change is claimed.',
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'pilots': len(results),
                      'all_above': all(x['pilot_rho'] > 2 ** -0.5 for x in results)}))


if __name__ == '__main__':
    main()
