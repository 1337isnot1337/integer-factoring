#!/usr/bin/env python3
"""Normalize bounded-source kernel slabs into low-degree crossing directions.

For a nonzero binary polynomial slab K and g=gcd(K's entries), the slab
(1+w)*K/g is a polynomial endpoint-zero direction.  This is an enlargement
of the available strict-crossing constructions; it need not stay in the
original bounded-source operator image.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


def divide(polynomial, divisor):
    if divisor <= 0:
        raise ValueError('The polynomial divisor must be nonzero.')
    quotient = 0
    while polynomial and polynomial.bit_length() >= divisor.bit_length():
        power = polynomial.bit_length()-divisor.bit_length()
        quotient ^= 1 << power
        polynomial ^= divisor << power
    return quotient, polynomial


def polynomial_gcd(first, second):
    while second:
        first, second = second, divide(first, second)[1]
    return first


def construct(source, degree=18):
    if source['dimensions'] != [5]*4:
        raise ValueError('This constructor uses the five-symbol bounded-source slab layout.')
    primitives = []
    excluded = []
    for index, slab in enumerate(source['kernel_slabs']):
        common = 0
        for entry in slab:
            common = polynomial_gcd(common, entry)
        if not common:
            continue
        normalized = []
        for entry in slab:
            quotient, remainder = divide(entry, common)
            assert remainder == 0
            normalized.append(quotient ^ (quotient << 1))
        assert all(value.bit_count() % 2 == 0 for value in normalized)
        source_degree = max(value.bit_length()-1 for value in slab)
        normalized_degree = max(value.bit_length()-1 for value in normalized)
        item = {'source_slab': index, 'source_degree': source_degree, 'entry_gcd': common,
                'entry_gcd_degree': common.bit_length()-1, 'normalized_degree': normalized_degree,
                'normalized_slab': normalized}
        (primitives if normalized_degree <= degree else excluded).append(item)
    descriptors = []
    for primitive_index, item in enumerate(primitives):
        for logical in range(5):
            for shift in range(degree-item['normalized_degree']+1):
                descriptors.append((primitive_index, logical, shift))
    return {'constructor': 'saturated_bounded_dual_kernel_directions', 'dimensions': [5]*4,
            'degree': degree, 'primitives': primitives, 'excluded_by_degree': excluded,
            'direction_descriptors': descriptors, 'direction_count_before_rank_reduction': len(descriptors),
            'scope': 'Every listed polynomial direction vanishes at the strict endpoint; division leaves the original operator image unrestricted.'}


def literal_directions(record):
    result = []
    for primitive_index, logical, shift in record['direction_descriptors']:
        slab = record['primitives'][primitive_index]['normalized_slab']
        direction = [(0, 0)]*625
        for j, value in enumerate(slab):
            direction[logical+5*j] = (value << shift, 0)
        result.append(tuple(direction))
    return tuple(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--degree', type=int, default=18)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = construct(json.loads(args.source.read_text()), args.degree)
    result.update(source_record=str(args.source), source_record_sha256=sha256(args.source.read_bytes()).hexdigest(),
                  constructor_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'primitives': len(result['primitives']),
                      'normalized_degrees': [item['normalized_degree'] for item in result['primitives']],
                      'directions_before_rank_reduction': result['direction_count_before_rank_reduction']}))


if __name__ == '__main__':
    main()
