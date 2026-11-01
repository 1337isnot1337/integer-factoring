#!/usr/bin/env python3
"""Exact all-size coefficient certificate for scalar Gaussian color clones."""
from __future__ import annotations
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, factorial, prod
from pathlib import Path


def transformed_coefficients(colors):
    polynomial = [0] * 13
    for r in range(5):
        polynomial[3*r] += comb(colors, r)
    for r in range(4):
        polynomial[4*r] -= comb(3, r) * colors**r
    assert polynomial[:3] == [0, 0, 0] and polynomial[3] == colors
    quotient = polynomial[3:]
    return [sum(quotient[i] * comb(9-i, j-i) for i in range(j+1)) for j in range(10)]


@lru_cache(None)
def partitions(total, minimum=1):
    if not total:
        return ((),)
    return tuple((first,) + tail for first in range(minimum, total+1)
                 for tail in partitions(total-first, first))


def blocks(colors):
    records = []
    for total in range(5):
        for partition in partitions(total):
            length = len(partition)
            if length > colors:
                continue
            multiplicity = prod(range(colors-length+1, colors+1))
            multiplicity //= prod(factorial(n) for n in Counter(partition).values())
            size = factorial(4) // (factorial(4-total) * prod(factorial(n) for n in partition))
            records.append({'partition': partition, 'multiplicity': multiplicity, 'block_size': size})
    assert sum(x['multiplicity'] * x['block_size'] for x in records) == (colors+1)**4
    assert all(x['block_size'] <= 24 for x in records)
    for r in range(min(colors, 4)+1):
        item = next(x for x in records if x['partition'] == (1,)*r)
        assert item['multiplicity'] == comb(colors, r)
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    samples = [transformed_coefficients(colors) for colors in range(33, 38)]
    differences = []
    for coefficient in range(10):
        values = [sample[coefficient] for sample in samples]
        row = []
        while values:
            row.append(values[0])
            values = [b-a for a, b in zip(values, values[1:])]
        assert row[0] > 0 and all(value >= 0 for value in row)
        differences.append(row)
    expected = [[33,1,0,0,0], [198,6,0,0,0], [396,12,0,0,0], [528,33,1,0,0],
                [1782,156,6,0,0], [1881,210,9,0,0], [176,300,29,1,0],
                [3102,825,78,3,0], [5973,963,81,3,0], [7634,2447,352,28,1]]
    assert differences == expected
    checks = 0
    for colors in (33,34,38,63,255,1023,65535,10**9):
        direct = transformed_coefficients(colors)
        newton = [sum(value * comb(colors-33, degree) for degree, value in enumerate(row)
                      if degree <= colors-33) for row in differences]
        assert direct == newton and min(direct) > 0
        checks += len(direct)
    all_blocks = {str(colors): blocks(colors) for colors in (1,2,3,4,7,33,255,10**9)}
    assert 3 * 3**2 > 4**2
    result = {'status': 'passed', 'minimum_colors': 33, 'minimum_physical_alphabet': 34,
              'strict_profile_floor': '3/4', 'quotient_degree': 9,
              'newton_expansion_base': 33, 'positive_coefficient_newton_table': differences,
              'additional_coefficient_controls': checks, 'color_blocks': all_blocks,
              'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'scope': 'Every homogeneous scalar Gaussian kernel replicated over independent colors; '
                       'vacuum/single-particle letters; common physical type and common-scale unweighted profiles.',
              'factoring_exponent': '1/2; unchanged'}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status':'passed','minimum_colors':33,'strict_profile_floor':'3/4',
                      'exact_newton_coefficients':50,'additional_controls':checks}, sort_keys=True))


if __name__ == '__main__':
    main()
