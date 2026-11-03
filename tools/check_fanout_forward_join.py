#!/usr/bin/env python3
"""Measure the source fan-out relation by forward evaluation.

The source edge tensor has one row per input digit vector and one column
per message vector.  Its entry is one exactly when the message vector is
the forward evaluation of the input.  Enumerating the Cartesian product
of inputs and messages is unnecessary: evaluating the function visits
one message per input.  This finite diagnostic counts the forward
messages and their multiplicities for the source fan-out at h=1,2.  It
does not compress the later evaluator pad or prove a factoring rate.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import itertools
import json
from pathlib import Path

from check_factoring_network_structure import Fanout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    rows = []
    for h in (1, 2):
        fan = Fanout(h, 2 * h * h)
        inputs = fan.base ** h
        cartesian = inputs * fan.modulus ** fan.a
        messages = Counter()
        for digits in itertools.product(range(fan.base), repeat=h):
            messages[fan.messages(digits)] += 1
        distinct = len(messages)
        if distinct != inputs:
            raise SystemExit('The full source fan-out was not injective.')
        rows.append({
            'h': h,
            'digit_bits': 2 * h * h,
            'base': fan.base,
            'wires': fan.a,
            'modulus': fan.modulus,
            'inputs': inputs,
            'distinct_message_vectors': distinct,
            'forward_work': inputs,
            'cartesian_message_pairs': cartesian,
            'forward_reduction_factor': cartesian // inputs,
            'max_messages_per_input': 1,
        })

    names = ('check_fanout_forward_join.py', 'check_factoring_network_structure.py')
    result = {
        'status': 'pass',
        'scope': (
            'Finite forward-evaluation diagnostic for the source fan-out at '
            'h=1,2. The full recipient map is injective, so one message is '
            'visited per input instead of the Cartesian input-message set. '
            'The later evaluator pad, full contracted tensor, and complete '
            'factoring bound are unchanged.'
        ),
        'rows': rows,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass',
                      'rows': [{'h': row['h'], 'inputs': row['inputs'],
                                'distinct_message_vectors': row['distinct_message_vectors'],
                                'forward_reduction_factor': row['forward_reduction_factor']}
                               for row in rows]}))


if __name__ == '__main__':
    main()
