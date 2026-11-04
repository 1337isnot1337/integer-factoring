#!/usr/bin/env python3
"""Verify the recorded q45 minimal deletion-closed F8 hyperoval core."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps
from search_coherent_hot_retractions import pressure

INPUT = Path('evidence/f8-hyperoval-minimal-core.json')
BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')
BOUNDARY_WORDS = (52, 76, 82, 97, 138, 145, 262, 265, 517, 1027)
SOURCE_NAMES = (
    'coherent_f8_maps.py',
    'coherent_hot_compounds.py',
    'coherent_hot_compounds_f8.py',
    'coherent_hot_compounds_f8.cpp',
    'selected_gaussian_coefficients.cpp',
    'native_rectangular_smith_f4.cpp',
    'search_coherent_hot_retractions.py',
)


def deletion_closure(words):
    closure = set()
    for word in words:
        closure.add(word)
        level = {word}
        for _ in range(word.bit_count()):
            level = {child for parent in level for bit in range(12)
                     if parent >> bit & 1 for child in (parent ^ (1 << bit),)}
            closure.update(level)
    return closure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved['completed'] is True
    assert saved['q'] == 45 and saved['modes'] == 11 and saved['cutoff'] == 127
    assert saved['base_profile_sha256'] == sha256(BASE.read_bytes()).hexdigest()
    closure = deletion_closure(BOUNDARY_WORDS)
    assert saved['boundary_words'] == list(BOUNDARY_WORDS)
    assert saved['closure_words'] == sorted(closure)
    assert len(closure) == 52

    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    single_indices = []
    for index, (left, right) in enumerate(zip(embedding, retraction)):
        if len(left) == 1 and len(right) == 1:
            word, coefficient = next(iter(left.items()))
            if coefficient == 1 and word in closure:
                single_indices.append(index)
    pair_indices = [index for index, left in enumerate(embedding) if len(left) > 1]
    assert single_indices == saved['single_source_indices']
    assert pair_indices == saved['pair_source_indices'] == [93, 94, 95]
    assert len(single_indices) == 42
    keep = single_indices + pair_indices
    embedding = [embedding[index] for index in keep]
    retraction = [retraction[index] for index in keep]
    heights, charges = validate_maps(embedding, retraction, 11)
    assert heights == saved['heights']
    assert charges == saved['charges']
    assert max(heights) == 25 and sum(heights) == 696
    assert saved['all_orders_resolved'] is True
    assert saved['tail_count'] == [0, 0]
    assert len(saved['zero_orders_t']) == 2025
    assert len(saved['infinity_orders_t']) == 2025
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
    assert abs(fresh['scale_t'] - saved['numerical_pressure']['scale_t']) < 1e-12
    assert fresh['rho'] > 2 ** -0.5
    assert fresh['rho'] > 0.7156

    for name in SOURCE_NAMES:
        assert sha256((Path('tools') / name).read_bytes()).hexdigest() == saved['source_sha256'][name]

    receipt = {
        'passed': True,
        'q': 45,
        'max_height': 25,
        'orders_per_cut': 2025,
        'tail_count': [0, 0],
        'LP_identity_and_deletion_closure': True,
        'numerical_pressure': fresh,
        'above_inverse_sqrt_two': True,
        'above_q96_q99_pressure': True,
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'base_profile_sha256': saved['base_profile_sha256'],
        'scope': (
            'Replay of the minimal deletion-closed q45 map and its recorded '
            'complete finite orders. The pressure is numerical and the complete '
            'exponent remains 1/2; deleting unused q96 core words does not cross.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'q': 45, 'rho': fresh['rho'], 'above': True}))


if __name__ == '__main__':
    main()
