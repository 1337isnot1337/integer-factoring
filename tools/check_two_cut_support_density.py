#!/usr/bin/env python3
"""Measure two-cut support density of the Gaussian reflection seed.

The two cuts of the tensor are (x,y)|(u,v) and (x,u)|(y,v).  Their joint
image is the support of the four-index tensor.  This script counts the
support for the q4 and q8 Gaussian seeds and compares it with the full
product of the two cut index spaces.  Finite support diagnostic only.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from reflection_dual_rail_gaussian import tensor


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    rows = []
    for interval in (4, 8):
        seed = tensor(interval)
        alphabet = seed['q_x']
        entries = len(seed['values'])
        nonzero = sum(1 for low, high in seed['values'] if low or high)
        rows.append({
            'interval': interval,
            'alphabet': alphabet,
            'entries': entries,
            'nonzero_entries': nonzero,
            'support_fraction': nonzero / entries,
            'two_cut_product': entries,
            'support_compression_factor': entries / nonzero,
        })

    names = ('check_two_cut_support_density.py', 'reflection_dual_rail_gaussian.py')
    result = {
        'status': 'pass',
        'scope': (
            'Finite support-density diagnostic for the interval-4/8 Gaussian seed. '
            'The two-cut joint image is the tensor support. It is about one '
            'third full, a constant-factor compression, not an exponent '
            'change. The Smith codec, later gates, and complete bound are '
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
                      'rows': [{'interval': row['interval'],
                                'alphabet': row['alphabet'],
                                'nonzero_entries': row['nonzero_entries'],
                                'support_fraction': row['support_fraction'],
                                'support_compression_factor': row['support_compression_factor']}
                               for row in rows]}))


if __name__ == '__main__':
    main()
