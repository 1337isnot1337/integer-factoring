#!/usr/bin/env python3
"""Build the smallest deletion-closed F8 hyperoval low core (q45)."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from coherent_f8_maps import validate_maps
from coherent_hot_compounds_f8 import partial_profiles
from search_coherent_hot_retractions import pressure

BASE = Path('evidence/f8-hyperoval-q96-full-profile.json')
BOUNDARY_WORDS = (52, 76, 82, 97, 138, 145, 262, 265, 517, 1027)
CUTOFF = 127
MODES = 11
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
        level = {word}
        closure.add(word)
        for _ in range(word.bit_count()):
            level = {child for parent in level for bit in range(12)
                     if parent >> bit & 1 for child in (parent ^ (1 << bit),)}
            closure.update(level)
    return closure


def load_maps():
    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    closure = deletion_closure(BOUNDARY_WORDS)
    single_indices = []
    for index, (left, right) in enumerate(zip(embedding, retraction)):
        if len(left) == 1 and len(right) == 1:
            word, coefficient = next(iter(left.items()))
            if coefficient == 1 and word in closure:
                single_indices.append(index)
    assert len(single_indices) == len(closure) - len(BOUNDARY_WORDS)
    pair_indices = [index for index, left in enumerate(embedding) if len(left) > 1]
    assert pair_indices == [93, 94, 95]
    keep = single_indices + pair_indices
    embedding = [embedding[index] for index in keep]
    retraction = [retraction[index] for index in keep]
    heights, charges = validate_maps(embedding, retraction, MODES)
    return base, closure, single_indices, pair_indices, embedding, retraction, heights, charges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/f8-hyperoval-minimal-core.json'))
    parser.add_argument('--raw-profile', type=Path,
                        help='Reuse an existing raw build instead of recomputing the profile.')
    args = parser.parse_args()

    base, closure, single_indices, pair_indices, embedding, retraction, heights, charges = load_maps()
    started = time.monotonic()
    if args.raw_profile:
        raw = json.loads(args.raw_profile.read_text())
        result = raw['construction']
    else:
        result = partial_profiles(embedding, retraction, MODES, CUTOFF,
                                  progress=lambda row: print(json.dumps(row), flush=True))
    assert result['all_orders_resolved'] is True
    assert all(cut['tail_count'] == 0 for cut in result['cuts'])
    assert [len(cut['retained_orders_t']) for cut in result['cuts']] == [2025, 2025]

    profile = {
        'q': 45,
        'modes': MODES,
        'heights': heights,
        'charges': charges,
        'zero_orders_t': result['cuts'][0]['retained_orders_t'],
        'infinity_orders_t': result['cuts'][1]['retained_orders_t'],
    }
    fit = pressure(profile)
    assert fit['completed'] is True
    assert fit['rho'] > 2 ** -0.5
    source_sha256 = {
        name: sha256((Path('tools') / name).read_bytes()).hexdigest()
        for name in SOURCE_NAMES
    }
    receipt = {
        'completed': True,
        'construction': (
            'The smallest deletion-closed quotient of the q96 F8 hyperoval map: '
            'all subsets of the ten height-25 boundary triples, plus the three '
            'non-single hyperoval embeddings and their retractions.'
        ),
        'boundary_words': list(BOUNDARY_WORDS),
        'closure_words': sorted(closure),
        'single_source_indices': single_indices,
        'pair_source_indices': pair_indices,
        'q': 45,
        'modes': MODES,
        'heights': heights,
        'charges': charges,
        'cutoff': CUTOFF,
        'zero_orders_t': profile['zero_orders_t'],
        'infinity_orders_t': profile['infinity_orders_t'],
        'tail_count': [cut['tail_count'] for cut in result['cuts']],
        'all_orders_resolved': True,
        'numerical_pressure': fit,
        'base_profile': str(BASE),
        'base_profile_sha256': sha256(BASE.read_bytes()).hexdigest(),
        'source_sha256': source_sha256,
        'elapsed_seconds': time.monotonic() - started,
        'scope': (
            'Complete finite natural profile of the smallest deletion-closed '
            'F8 hyperoval low core. The pressure is numerical, not an exact '
            'certificate. It is above 1/sqrt(2) and above the q96/q99 pressure, '
            'so deleting the unused q96 core words does not cross; the complete '
            'factoring exponent is unchanged.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'completed': True, 'q': 45, 'pressure': fit, 'tail_count': [0, 0]}))


if __name__ == '__main__':
    main()
