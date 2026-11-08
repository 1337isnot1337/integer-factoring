#!/usr/bin/env python3
"""Replay the q115 height-gauge F8 candidate crossing at cutoff 47."""
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

INPUT = Path('evidence/f8-q115-height-gauge-crossing.json')
BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')
SOURCES = [
    'coherent_f8_maps.py',
    'coherent_hot_compounds_f8.py',
    'coherent_hot_compounds_f8.cpp',
    'probe_f8_hyperoval_retraction.py',
    'probe_determinant_inverse_codes.py',
]
CRITICAL = 1 / math.sqrt(2)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved['base_profile_sha256'] == sha256(BASE.read_bytes()).hexdigest()
    assert saved['source_sha256'] == {
        name: sha256((Path(__file__).resolve().parents[1] / 'tools' / name).read_bytes()).hexdigest()
        for name in SOURCES
    }

    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    weights = saved['weights']
    assert weights == [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 11]

    support = set()
    for vector in embedding + retraction:
        support.update(vector)
    low = [word for word in range(1, 1 << 12)
           if sum(weights[i] for i in range(12) if word >> i & 1) < 25]
    missing = [word for word in low if word not in support]
    assert missing == saved['new_low_words']
    assert len(missing) == 19
    embedding = embedding + [{word: 1} for word in missing]
    retraction = retraction + [{word: 1} for word in missing]
    assert len(embedding) == saved['q'] == 115

    standard_heights, standard_charges = validate_maps(embedding, retraction, 12)
    assert standard_heights == saved['standard_heights']
    assert standard_charges == saved['standard_charges']
    custom = []
    for left, right in zip(embedding, retraction):
        heights = {sum(weights[i] for i in range(12) if word >> i & 1)
                   for word in set(left) | set(right)}
        assert len(heights) == 1
        custom.append(next(iter(heights)))
    assert custom == saved['custom_heights']
    assert min(custom) == 0 and max(custom) == 25
    support = set().union(*[set(v) for v in embedding + retraction])
    assert all(word in support for word in low)
    assert len(low) == saved['low_core_count'] == 110

    orders = [saved['zero_orders_t'], saved['infinity_orders_t']]
    assert [len(cut) for cut in orders] == saved['retained'] == [935, 935]
    known = fit({'cuts': [{'retained_orders_t': cut} for cut in orders]}, custom)
    augmented = [
        sorted(cut + [saved['cutoff'] + 1] * saved['tails'][index])
        for index, cut in enumerate(orders)
    ]
    upper = fit({'cuts': [{'retained_orders_t': cut} for cut in augmented]}, custom)
    assert abs(known['rho'] - saved['known']['rho']) < 1e-14
    assert abs(upper['rho'] - saved['upper']['rho']) < 1e-14
    assert upper['rho'] < CRITICAL
    assert upper['strictly_below_critical_numerically'] is True
    assert known['rho'] <= upper['rho']

    receipt = {
        'status': 'PASS',
        'q': 115,
        'cutoff': 47,
        'retained_orders_per_cut': 935,
        'unresolved_tails_per_cut': 12290,
        'new_low_words': len(missing),
        'low_core_words': len(low),
        'known_order_lower_rho': known['rho'],
        'tail_augmented_upper_rho': upper['rho'],
        'critical_rho': CRITICAL,
        'upper_margin': CRITICAL - upper['rho'],
        'candidate_crossing_numerically': True,
        'optimizer_globally_certified': False,
        'complete_profile_recomputed_this_run': False,
        'uniform_transfer_proved': False,
        'factoring_exponent_improved': False,
        'scope': saved['scope'],
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'candidate_crossing': True,
                      'known_lower_rho': known['rho'],
                      'tail_augmented_upper_rho': upper['rho'],
                      'margin': CRITICAL - upper['rho']}))


if __name__ == '__main__':
    main()
