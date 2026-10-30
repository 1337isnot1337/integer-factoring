#!/usr/bin/env python3
"""Both-cut Ward bounds using the one-sided charge gauge."""

import argparse
import hashlib
import json
from math import isqrt, log2
from pathlib import Path

from check_modular_crossing import rank_two
from screen_ward_code_bounds import evaluate
from ward_rectangular_bounds import rectangular_bounds, ward_score


def channel_free_score(pairs, length, forms, labels, total_gram_rank=None):
    old = ward_score(pairs, length, forms, labels)
    previous = {(a, b): (r, c) for a, b, r, c in old['profile']}
    column = {}
    largest = old['largest_explicit_block']
    for j in range(forms + 1):
        for shift in range(forms - j + 1):
            for state_shift in range(-shift, shift + 1, 2):
                grades = (length - 2 * j - shift, length - state_shift,
                          length, length - 2 * forms)
                profile, size = rectangular_bounds(pairs, *grades)
                largest = max(largest, size)
                for (b, a), rank in profile.items():
                    target = a, b + forms + j
                    if sum(target) <= length:
                        column[target] = column.get(target, 0) + rank
    order = forms if total_gram_rank is None else max(forms, 2 * forms - total_gram_rank)
    total, result = 0, []
    for (a, b), (row, old_column) in sorted(previous.items()):
        if a + b < order:
            continue
        col = min(old_column, column.get((a, b), 0))
        if row or col:
            value = row * col
            root = isqrt(value)
            total += root + int(root * root < value)
            result.append([a, b, row, col])
    return {'pairs': pairs, 'length': length, 'forms': forms, 'labels': labels,
            'joint_order_used': order, 'total_form_Gram_rank': total_gram_rank,
            'integer_rank_sum_upper_bound': total,
            'previous_integer_rank_sum_upper_bound': old['integer_rank_sum_upper_bound'],
            'rho_of_upper_bound': log2(total) / (2 * log2(labels)),
            'largest_explicit_block': largest, 'profile': result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append', default=[], help='block_size:forms:length')
    parser.add_argument('--half', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cases = args.case or ['4:2:8', '4:3:10', '8:2:12', '8:3:18']
    with args.output.open('w') as output:
        for case in cases:
            block_size, forms, length = map(int, case.split(':'))
            base = evaluate(block_size, forms, length, full_code=not args.half)
            # Disjoint indicator forms have total Gram equal to (b mod2)Id.
            gram_rank = forms * (block_size % 2)
            result = channel_free_score(block_size * forms, length, forms, base['labels'], gram_rank)
            result.update(block_size=block_size, orientation_half=args.half)
            result['scope'] = 'Certified upper bound for the stated admissible-syndrome code; no obstruction inferred from failure.'
            result['source_sha256'] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                       for name in (Path(__file__).name, 'ward_rectangular_bounds.py')}
            output.write(json.dumps(result, sort_keys=True) + '\n')
            output.flush()
            print(json.dumps({key: value for key, value in result.items()
                              if key not in ('profile', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
