#!/usr/bin/env python3
"""Measure joint boundary compression in the source Fourier-coset codec.

For a proper Fourier-coset bundle the row and column codecs are both
functions of the same input digits.  The set of joint (row, column)
codes is therefore the image of one map, not the Cartesian product of
the two code spaces.  This script enumerates the source inputs at h=1,2
and compares the joint image size with the product of the code ranks.
It is a finite codec diagnostic, not a full two-cut pad theorem.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import itertools
import json
from pathlib import Path

from check_factoring_network_structure import Bundle, Fanout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    rows = []
    for h in (1, 2):
        fan = Fanout(h, 2 * h * h)
        bundle = Bundle(fan, 2, 0)
        inputs = fan.base ** h
        joint = set()
        row_codes = set()
        column_codes = set()
        missing = 0
        for digits in itertools.product(range(fan.base), repeat=h):
            row = bundle.row_code(digits)
            messages = tuple(fan.messages(digits)[j] for j in bundle.indices)
            column = bundle.column_code(messages)
            if column is None:
                missing += 1
                continue
            joint.add((row, column))
            row_codes.add(row)
            column_codes.add(column)
        if missing:
            raise SystemExit('A bundle column code was undefined.')
        if len(joint) != inputs:
            raise SystemExit('The joint bundle map was not injective.')
        ambient = bundle.rank * bundle.rank
        if len(joint) != bundle.rank:
            raise SystemExit('The joint image did not equal the bundle rank.')
        rows.append({
            'h': h,
            'p': bundle.p,
            'coset': bundle.c,
            'inputs': inputs,
            'bundle_rank': bundle.rank,
            'distinct_row_codes': len(row_codes),
            'distinct_column_codes': len(column_codes),
            'joint_image': len(joint),
            'ambient_code_product': ambient,
            'joint_compression_factor': ambient // len(joint),
        })

    names = ('check_bundle_joint_compression.py', 'check_factoring_network_structure.py')
    result = {
        'status': 'pass',
        'scope': (
            'Finite source-codec diagnostic at h=1,2 for the p=2 Fourier '
            'coset. The joint row/column code image is one input-indexed '
            'table, not the product of the two code spaces. The full '
            'two-cut pad, later gates, and complete factoring bound are '
            'unchanged.'
        ),
        'rows': rows,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass',
                      'rows': [{'h': row['h'], 'bundle_rank': row['bundle_rank'],
                                'joint_image': row['joint_image'],
                                'joint_compression_factor': row['joint_compression_factor']}
                               for row in rows]}))


if __name__ == '__main__':
    main()
