#!/usr/bin/env python3
"""Complete height-39 extension of the q32 height-16 coherent map."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_coherent_q32_pressure import displayed_maps
from coherent_hot_compounds import profiles
from search_coherent_hot_retractions import pressure

EXTENSION = (96, 160, 192, 224)
EXTENSION_HEIGHTS = (24, 26, 28, 39)
THRESHOLD = 1 / math.sqrt(2)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    embedding, retraction = displayed_maps(True)
    assert len(embedding) == len(retraction) == 32
    extended_embedding = list(embedding) + [1 << word for word in EXTENSION]
    extended_retraction = list(retraction) + [1 << word for word in EXTENSION]
    baseline = profiles(embedding, retraction, 12)
    extended = profiles(extended_embedding, extended_retraction, 12)
    assert extended['heights'][-4:] == list(EXTENSION_HEIGHTS)
    assert len(extended['zero_orders_t']) == len(extended['infinity_orders_t']) == 36 * 36
    assert extended['q'] == 36 and max(extended['heights']) == 39

    base_fit = pressure(baseline)
    fit = pressure(extended)
    assert base_fit['completed'] and fit['completed']
    assert fit['rho'] > THRESHOLD
    assert fit['rho'] < base_fit['rho']

    receipt = {
        'completed': True,
        'base': 'q32 binary height-16 coherent winner',
        'base_q': 32,
        'base_maximum_height': 16,
        'extension_words': list(EXTENSION),
        'extension_heights': list(EXTENSION_HEIGHTS),
        'q': 36,
        'modes': 12,
        'maximum_height': 39,
        'order_counts': [36 * 36, 36 * 36],
        'tail_counts': [0, 0],
        'baseline_numerical_pressure': base_fit,
        'numerical_pressure': fit,
        'above_inverse_sqrt_two': fit['rho'] > THRESHOLD,
        'gain_over_baseline': base_fit['rho'] - fit['rho'],
        'scope': (
            'Complete finite profile of the q32 map plus four one-hot height-24/26/28/39 '
            'labels. Numerical pressure only; not an exact certificate, a split, or a '
            'score change.'
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({
        'q': 36,
        'baseline_pressure': base_fit['rho'],
        'pressure': fit['rho'],
        'gain': base_fit['rho'] - fit['rho'],
    }))


if __name__ == '__main__':
    main()
