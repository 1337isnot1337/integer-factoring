#!/usr/bin/env python3
"""Compare Cartesian enumeration and a hash join on the identity minor.

The source's finite identity-minor check enumerates every ordered pair of
messages and compares them.  The same sum is an equality join, so a hash
table computes it in one pass over the messages plus one pass over the
output.  This script reproduces both computations exactly for s=0,1 and
counts the work.  It is a finite source-definition diagnostic, not a
certificate for the full tensor or a factoring rate.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path


def messages_for(a, b, h, A, modulus):
    values = 1 << (2 * h * h)
    powers = [pow(16, A * i * b, modulus) for i in range(A)]
    base = pow(16, a * b, modulus)
    output = []
    for v in range(values):
        digits = [v] + [0] * (h - 1)
        total = sum(digits[i] * powers[i] for i in range(h))
        output.append(base * total % modulus)
    return output


def brute_join(messages):
    count = 0
    for left in messages:
        for right in messages:
            count += left == right
    return count


def hash_join(messages):
    counts = Counter(messages)
    return sum(multiplicity * multiplicity for multiplicity in counts.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()

    rows = []
    for s in (0, 1):
        h = 1 << s
        A = 2 * h
        d = 2 * h * h
        modulus = (1 << (4 * d)) + 1
        values = 1 << d
        brute_pairs = 0
        brute_answer = 0
        hash_answer = 0
        messages_built = 0
        output_size = 0
        for a in range(A):
            for b in range(A):
                messages = messages_for(a, b, h, A, modulus)
                messages_built += len(messages)
                current = brute_join(messages)
                joined = hash_join(messages)
                if current != joined:
                    raise SystemExit('Hash join differs from Cartesian enumeration.')
                brute_answer += current
                hash_answer += joined
                brute_pairs += len(messages) ** 2
                output_size += current
        if brute_answer != hash_answer:
            raise SystemExit('Aggregate answers differ.')
        if brute_answer != A * A * values:
            raise SystemExit('The identity minor did not produce the expected permutation count.')
        rows.append({
            's': s,
            'h': h,
            'A': A,
            'd': d,
            'modulus': modulus,
            'wires': A * A,
            'values_per_wire': values,
            'brute_pairs': brute_pairs,
            'messages_built': messages_built,
            'hash_output_pairs': output_size,
            'hash_join_work': messages_built + output_size,
            'pair_reduction_factor': brute_pairs / (messages_built + output_size),
            'answer': brute_answer,
        })

    if rows[1]['modulus'] != 641 * 6700417:
        raise SystemExit('The s=1 composite modulus changed.')

    names = ('check_hash_join_identity_minor.py',)
    result = {
        'status': 'pass',
        'scope': (
            'Finite identity-minor source diagnostic. The equality join is '
            'computed exactly in linear message-plus-output work on these '
            'cases. This does not assert that the full edge tensors are '
            'functional, does not compress the evaluator pad, and does not '
            'change the factoring exponent.'
        ),
        'rows': rows,
        'source_sha256': {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in names},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'pass',
                      'rows': [{'s': row['s'],
                                'brute_pairs': row['brute_pairs'],
                                'hash_join_work': row['hash_join_work'],
                                'pair_reduction_factor': row['pair_reduction_factor']}
                               for row in rows]}))


if __name__ == '__main__':
    main()
