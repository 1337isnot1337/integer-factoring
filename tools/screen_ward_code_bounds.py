#!/usr/bin/env python3
"""Evaluate certified Ward coefficient bounds for global or disjoint-form codes."""

import argparse
import hashlib
import json
from math import comb
from pathlib import Path

from ward_rectangular_bounds import ward_score


def radical_dimensions(pairs):
    result = []
    for length in range(2 * pairs + 1):
        total = 0
        for singles in range(pairs + 1):
            if (length - singles) % 2 or (pairs - singles) % 2:
                continue
            even, grade = pairs - singles, (length - singles) // 2
            if even >= 2 and 0 <= grade - 1 <= even - 2:
                total += comb(pairs, singles) * 2 ** singles * comb(even - 2, grade - 1)
        result.append(total)
    return result


def multiply(first, second):
    result = [0] * (len(first) + len(second) - 1)
    for i, x in enumerate(first):
        for j, y in enumerate(second):
            result[i + j] += x * y
    return result


def evaluate(block_size, forms, length, full_code=False):
    dimensions, fixed = [1], [1]
    local = radical_dimensions(block_size)
    fixed_local = [comb(block_size - 2, h - 1) if 1 <= h <= block_size - 1 else 0
                   for h in range(block_size + 1)] if block_size % 2 == 0 else [0]
    for _ in range(forms):
        dimensions = multiply(dimensions, local)
        fixed = multiply(fixed, fixed_local)
    zero = fixed[length // 2] if length % 2 == 0 and length // 2 < len(fixed) else 0
    labels = dimensions[length] if full_code else (dimensions[length] - zero) // 2
    if labels <= 1:
        raise ValueError('at least two orientation-half labels required')
    # Indicator rows on disjoint blocks: an independent ell-subset of weight
    # columns chooses ell distinct blocks, one site per block, and two modes.
    channels = [2 ** ell * comb(forms, ell) * block_size ** ell for ell in range(forms + 1)]
    result = ward_score(block_size * forms, length, forms, labels, channels)
    result.update(block_size=block_size, channel_counts=channels, orientation_half=not full_code)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append', default=[], help='block_size:forms:length')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--full-code', action='store_true', help='retain the full homogeneous product image')
    args = parser.parse_args()
    cases = args.case or ['16:1:10', '24:1:14', '32:1:18', '4:2:8', '8:3:18']
    source = Path(__file__).with_name('ward_rectangular_bounds.py')
    with args.output.open('w') as output:
        for case in cases:
            result = evaluate(*map(int, case.split(':')), full_code=args.full_code)
            result['bound_engine_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
            result['screen_source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            result['scope'] = 'Computed upper bound; failure to certify is not an actual-rank obstruction.'
            output.write(json.dumps(result, sort_keys=True) + '\n')
            output.flush()
            print(json.dumps({k: v for k, v in result.items() if k != 'profile'}), flush=True)


if __name__ == '__main__':
    main()
