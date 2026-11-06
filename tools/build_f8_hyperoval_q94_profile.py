#!/usr/bin/env python3
"""Build the complete q94 F8 hyperoval profile with mode 11 removed."""
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
CUTOFF = 135
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


def load_maps():
    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]

    def uses_removed_mode(table):
        return any(word >> 11 & 1 for word in table)

    keep = [index for index, (left, right) in enumerate(zip(embedding, retraction))
            if not uses_removed_mode(left) and not uses_removed_mode(right)]
    embedding = [embedding[index] for index in keep]
    retraction = [retraction[index] for index in keep]
    heights, charges = validate_maps(embedding, retraction, MODES)
    return base, keep, embedding, retraction, heights, charges


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path('evidence/f8-hyperoval-q94-profile.json'))
    args = parser.parse_args()

    base, keep, embedding, retraction, heights, charges = load_maps()
    assert len(embedding) == 94
    assert max(heights) == 25 and sum(heights) == 1656
    started = time.monotonic()
    result = partial_profiles(embedding, retraction, MODES, CUTOFF,
                              progress=lambda row: print(json.dumps(row), flush=True))
    assert result['all_orders_resolved'] is True
    assert all(cut['tail_count'] == 0 for cut in result['cuts'])
    assert [len(cut['retained_orders_t']) for cut in result['cuts']] == [8836, 8836]
    profile = {
        'q': 94,
        'modes': MODES,
        'heights': heights,
        'charges': charges,
        'zero_orders_t': result['cuts'][0]['retained_orders_t'],
        'infinity_orders_t': result['cuts'][1]['retained_orders_t'],
    }
    fit = pressure(profile)
    assert fit['completed'] is True
    assert fit['rho'] > 2 ** -0.5
    assert fit['rho'] > 29 / 41
    source_sha256 = {name: sha256((Path('tools') / name).read_bytes()).hexdigest()
                     for name in SOURCE_NAMES}
    receipt = {
        'completed': True,
        'construction': (
            'The q96 F8 hyperoval map with mode 11 (height 23) removed from '
            'both the embedding and retraction supports. All remaining '
            'coordinates are retained; the ten canonical height-25 boundary '
            'triples remain because they do not use mode 11.'
        ),
        'base_profile': str(BASE),
        'base_profile_sha256': sha256(BASE.read_bytes()).hexdigest(),
        'kept_source_indices': keep,
        'q': 94,
        'modes': MODES,
        'heights': heights,
        'charges': charges,
        'cutoff': CUTOFF,
        'zero_orders_t': profile['zero_orders_t'],
        'infinity_orders_t': profile['infinity_orders_t'],
        'tail_count': [cut['tail_count'] for cut in result['cuts']],
        'all_orders_resolved': True,
        'numerical_pressure': fit,
        'source_sha256': source_sha256,
        'elapsed_seconds': time.monotonic() - started,
        'scope': (
            'Complete finite natural profile of the q94 map obtained by '
            'dropping mode 11 from q96. The pressure is numerical, not an '
            'exact certificate. The map removes part of the full height-below-24 '
            'core, so the height-38 floor theorem does not apply; its pressure '
            'is nevertheless above 29/41 and 1/sqrt(2). The complete factoring '
            'exponent is unchanged.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'completed': True, 'q': 94, 'pressure': fit,
                      'tail_count': [0, 0]}))


if __name__ == '__main__':
    main()
