#!/usr/bin/env python3
"""Compose a shared-input pair with one independent bundle.

The shared pair on one fan has a graph image of one pair per input.
The third bundle uses independent digits, so the joint image is the
product of the shared graph and the third rank.  This is the component
rule needed for a nonuniform pad analysis.  h=1 is enumerated exactly;
h=2 uses the component rule rather than enumerating 2^32 pairs.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import itertools
import json
from pathlib import Path

from check_factoring_network_structure import Bundle, Fanout


def pair(bundle, fan, digits):
    row = bundle.row_code(digits)
    messages = tuple(fan.messages(digits)[j] for j in bundle.indices)
    column = bundle.column_code(messages)
    if column is None:
        raise SystemExit('A bundle column code was undefined.')
    return row, column


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    rows = []
    for h in (1, 2):
        shared_fan = Fanout(h, 2 * h * h)
        left, right = Bundle(shared_fan, 2, 0), Bundle(shared_fan, 2, 1)
        independent_fan = Fanout(h, 2 * h * h)
        third = Bundle(independent_fan, 2, 0)
        shared_image = shared_fan.base ** h
        joint = shared_image * third.rank
        full_product = left.rank * right.rank * third.rank
        if h == 1:
            enumerated = set()
            for digits in itertools.product(range(shared_fan.base), repeat=h):
                for other in itertools.product(range(independent_fan.base), repeat=h):
                    enumerated.add((pair(left, shared_fan, digits),
                                    pair(right, shared_fan, digits),
                                    pair(third, independent_fan, other)))
            if len(enumerated) != joint:
                raise SystemExit('The h=1 component rule did not replay.')
        rows.append({
            'h': h,
            'shared_inputs': shared_fan.base ** h,
            'shared_pair_image': shared_image,
            'independent_rank': third.rank,
            'joint_image': joint,
            'full_rank_product': full_product,
            'compression_factor': full_product // joint,
            'h1_enumerated': h == 1,
        })

    names = ('check_three_bundle_joint_image.py', 'check_factoring_network_structure.py')
    result = {
        'status': 'pass',
        'scope': (
            'Finite three-bundle codec diagnostic. The joint image of a '
            'shared-input pair and one independent bundle is the product of '
            'the shared graph and the independent rank. h=1 is exact by '
            'enumeration; h=2 uses the component rule. The full two-cut pad '
            'and complete factoring bound are unchanged.'
        ),
        'rows': rows,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass',
                      'rows': [{'h': row['h'], 'joint_image': row['joint_image'],
                                'full_rank_product': row['full_rank_product'],
                                'compression_factor': row['compression_factor'],
                                'h1_enumerated': row['h1_enumerated']}
                               for row in rows]}))


if __name__ == '__main__':
    main()
