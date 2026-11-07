#!/usr/bin/env python3
"""Build the cutoff-47 F8 low-block word-31-for-46 swap pilot."""
from __future__ import annotations

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
OUTPUT = Path('evidence/f8-low-block-swap-pilot.json')
CUTOFF = 47
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
    base = json.loads(BASE.read_text())
    embedding = [{int(k): v for k, v in row.items()} for row in base['embedding']]
    retraction = [{int(k): v for k, v in row.items()} for row in base['retraction']]
    keep = [index for index, (left, right) in enumerate(zip(embedding, retraction))
            if not (len(left) == 1 and len(right) == 1 and next(iter(left)) == 31)]
    embedding = [embedding[index] for index in keep] + [{46: 1}]
    retraction = [retraction[index] for index in keep] + [{46: 1}]
    heights, charges = validate_maps(embedding, retraction, 12)
    started = time.time()
    result = partial_profiles(embedding, retraction, 12, CUTOFF)
    profile = {
        'q': len(embedding),
        'modes': 12,
        'heights': heights,
        'charges': charges,
        'zero_orders_t': result['cuts'][0]['retained_orders_t'],
        'infinity_orders_t': result['cuts'][1]['retained_orders_t'],
    }
    saved = {
        'q': len(embedding),
        'modes': 12,
        'cutoff': CUTOFF,
        'base_profile_sha256': sha256(BASE.read_bytes()).hexdigest(),
        'dropped_word': 31,
        'added_word': 46,
        'kept_source_indices': keep,
        'heights': heights,
        'charges': charges,
        'zero_orders_t': profile['zero_orders_t'],
        'infinity_orders_t': profile['infinity_orders_t'],
        'tail_count': [cut['tail_count'] for cut in result['cuts']],
        'all_orders_resolved': result['all_orders_resolved'],
        'numerical_pressure': pressure(profile),
        'source_sha256': {
            name: sha256((Path('tools') / name).read_bytes()).hexdigest()
            for name in SOURCE_NAMES
        },
    }
    OUTPUT.write_text(json.dumps(saved, indent=2) + '\n')
    print(json.dumps({'q': saved['q'], 'max_height': max(heights),
                      'tail_count': saved['tail_count'],
                      'rho': saved['numerical_pressure']['rho'],
                      'seconds': round(time.time() - started, 1)}))


if __name__ == '__main__':
    main()
