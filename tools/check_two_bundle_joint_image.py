#!/usr/bin/env python3
"""Compare two-bundle joint images with shared and disjoint input digits.

The proper-coset bundle compresses its row/column code product because
both codecs see the same digits.  This script tests whether that
compression survives across two bundles.  Shared digits give a graph
image; disjoint digits give the full product of the two bundle ranks.
Finite codec diagnostic only; no full pad or factoring claim.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import itertools
import json
from pathlib import Path

from check_factoring_network_structure import Bundle, Fanout


def bundle_pair(bundle, fan, digits):
    row = bundle.row_code(digits)
    messages = tuple(fan.messages(digits)[j] for j in bundle.indices)
    column = bundle.column_code(messages)
    if column is None:
        raise SystemExit('A bundle column code was undefined.')
    return row, column


def same_input_image(h):
    fan = Fanout(h, 2 * h * h)
    left, right = Bundle(fan, 2, 0), Bundle(fan, 2, 1)
    image = set()
    for digits in itertools.product(range(fan.base), repeat=h):
        image.add((bundle_pair(left, fan, digits),
                   bundle_pair(right, fan, digits)))
    return fan, left, right, len(image)


def disjoint_input_image(h):
    left_fan = Fanout(h, 2 * h * h)
    right_fan = Fanout(h, 2 * h * h)
    left, right = Bundle(left_fan, 2, 0), Bundle(right_fan, 2, 0)
    image = set()
    for left_digits in itertools.product(range(left_fan.base), repeat=h):
        for right_digits in itertools.product(range(right_fan.base), repeat=h):
            image.add((bundle_pair(left, left_fan, left_digits),
                       bundle_pair(right, right_fan, right_digits)))
    return left_fan, right_fan, left, right, len(image)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    rows = []
    for h in (1, 2):
        fan, left, right, image = same_input_image(h)
        product = left.rank * right.rank
        if image != fan.base ** h:
            raise SystemExit('Shared-input joint image did not equal the input count.')
        rows.append({
            'h': h,
            'case': 'shared_inputs',
            'inputs_each': fan.base ** h,
            'left_rank': left.rank,
            'right_rank': right.rank,
            'rank_product': product,
            'joint_image': image,
            'compression_factor': product // image,
        })
    h = 1
    left_fan, right_fan, left, right, image = disjoint_input_image(h)
    product = left.rank * right.rank
    if image != product:
        raise SystemExit('Disjoint-input joint image was not the full product.')
    rows.append({
        'h': h,
        'case': 'disjoint_inputs',
        'inputs_each': left_fan.base ** h,
        'left_rank': left.rank,
        'right_rank': right.rank,
        'rank_product': product,
        'joint_image': image,
        'compression_factor': product // image,
    })

    names = ('check_two_bundle_joint_image.py', 'check_factoring_network_structure.py')
    result = {
        'status': 'pass',
        'scope': (
            'Finite two-bundle codec diagnostic. Shared input digits give a '
            'graph image of the same size as the input space; disjoint input '
            'digits give the full rank product at h=1. The full two-cut pad '
            'composition and complete factoring bound are unchanged.'
        ),
        'rows': rows,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass',
                      'rows': [{'h': row['h'], 'case': row['case'],
                                'joint_image': row['joint_image'],
                                'rank_product': row['rank_product'],
                                'compression_factor': row['compression_factor']}
                               for row in rows]}))


if __name__ == '__main__':
    main()
