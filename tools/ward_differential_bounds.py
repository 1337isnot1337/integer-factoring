#!/usr/bin/env python3
"""Global-form Ward bounds after replacing the Clifford operator by a derivative."""

import argparse
import hashlib
import json
from math import isqrt, log2
from pathlib import Path

from screen_ward_code_bounds import radical_dimensions
from ward_rectangular_bounds import rectangular_bounds, ward_score


def differential_score(pairs, length, labels):
    specifications = [(length, length - 2, length, length),
                      (length - 2, length - 2, length, length),
                      (length, length, length, length - 2),
                      (length - 2, length, length, length - 2)]
    profiles = []
    largest = 0
    for grades in specifications:
        profile, size = rectangular_bounds(pairs, *grades)
        profiles.append(profile)
        largest = max(largest, size)
    old = ward_score(pairs, length, 1, labels)
    old_profile = {(a, b): (r, c) for a, b, r, c in old['profile']}
    answer, total = [], 0
    for a in range(length + 1):
        for b in range(length + 1 - a):
            parity = (length + a + b - 1) % 2
            row = parity * profiles[0].get((a, b - 1), 0) + profiles[1].get((a, b - 2), 0)
            col = parity * profiles[2].get((b - 1, a), 0) + profiles[3].get((b - 2, a), 0)
            old_row, old_col = old_profile.get((a, b), (0, 0))
            row, col = min(row, old_row), min(col, old_col)
            if row or col:
                product = row * col
                root = isqrt(product)
                total += root + int(root * root < product)
                answer.append([a, b, row, col])
    return {'pairs': pairs, 'length': length, 'labels': labels,
            'integer_rank_sum_upper_bound': total, 'profile': answer,
            'previous_integer_rank_sum_upper_bound': old['integer_rank_sum_upper_bound'],
            'rho_of_upper_bound': log2(total) / (2 * log2(labels)),
            'largest_explicit_block': max(largest, old['largest_explicit_block'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append', default=[])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cases = args.case or ['4:4', '16:10', '24:14', '32:18']
    with args.output.open('w') as output:
        for case in cases:
            pairs, length = map(int, case.split(':'))
            result = differential_score(pairs, length, radical_dimensions(pairs)[length])
            result['scope'] = 'Certified coefficient-cost upper bound; no lower-bound inference from failure.'
            result['source_sha256'] = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                       for name in (Path(__file__).name, 'ward_rectangular_bounds.py')}
            output.write(json.dumps(result, sort_keys=True) + '\n')
            output.flush()
            print(json.dumps({k: v for k, v in result.items() if k not in ('profile', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    main()
