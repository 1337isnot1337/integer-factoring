#!/usr/bin/env python3
"""Build the complete q99 height-39 extension of the F8 hyperoval profile."""
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
EXTENSION_WORDS = (160, 192, 224)
CUTOFF = 191
SOURCE_NAMES = (
    'coherent_f8_maps.py',
    'coherent_hot_compounds.py',
    'coherent_hot_compounds_f8.py',
    'coherent_hot_compounds_f8.cpp',
    'selected_gaussian_coefficients.cpp',
    'native_rectangular_smith_f4.cpp',
    'search_coherent_hot_retractions.py',
)


def load_maps():
    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    for word in EXTENSION_WORDS:
        embedding.append({word: 1})
        retraction.append({word: 1})
    heights, charges = validate_maps(embedding, retraction, 12)
    assert heights[-3:] == [26, 28, 39]
    return base, embedding, retraction, heights, charges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/f8-hyperoval-height39-extension.json'))
    parser.add_argument('--raw-profile', type=Path,
                        help='Reuse an existing raw build instead of recomputing the profile.')
    args = parser.parse_args()

    base, embedding, retraction, heights, charges = load_maps()
    started = time.monotonic()
    if args.raw_profile:
        raw = json.loads(args.raw_profile.read_text())
        result = raw['construction']
    else:
        result = partial_profiles(embedding, retraction, 12, CUTOFF,
                                  progress=lambda row: print(json.dumps(row), flush=True))
    assert result['all_orders_resolved'] is True
    assert all(cut['tail_count'] == 0 for cut in result['cuts'])
    assert [len(cut['retained_orders_t']) for cut in result['cuts']] == [9801, 9801]

    profile = {
        'q': 99,
        'modes': 12,
        'heights': heights,
        'charges': charges,
        'zero_orders_t': result['cuts'][0]['retained_orders_t'],
        'infinity_orders_t': result['cuts'][1]['retained_orders_t'],
    }
    fit = pressure(profile)
    assert fit['completed'] is True
    source_sha256 = {
        name: sha256((Path('tools') / name).read_bytes()).hexdigest()
        for name in SOURCE_NAMES
    }
    receipt = {
        'completed': True,
        'construction': 'The q96 F8 hyperoval map plus the three single-word coordinates 160, 192, 224 in both embedding and retraction.',
        'extension_words': list(EXTENSION_WORDS),
        'extension_heights': [26, 28, 39],
        'q': 99,
        'modes': 12,
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
            'Complete finite natural profile of the specified q99 height-39 extension; '
            'the pressure is numerical, not an exact certificate. It remains above '
            '1/sqrt(2), and the complete factoring exponent is unchanged.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'completed': True, 'q': 99, 'pressure': fit, 'tail_count': [0, 0]}))


if __name__ == '__main__':
    main()
