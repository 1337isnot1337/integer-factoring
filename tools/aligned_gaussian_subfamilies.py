#!/usr/bin/env python3
"""Reproducible raw-column and independent-layer Gaussian search exports."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from saturated_dual_kernel_directions import literal_directions


def raw_family(source, path, source_hash):
    primitives = []
    for row in source['raw_slabs']:
        if row['direct_degree18_direction']:
            slab = [value >> 8 for value in row['raw_slab_in_z']]
            primitives.append({'complement_mask': row['complement_mask'], 'normalized_slab': slab,
                               'normalized_degree': max(value.bit_length()-1 for value in slab)})
    return {'constructor': 'aligned_gaussian_raw_complement_columns', 'dimensions': [5]*4, 'degree': 18,
            'primitives': primitives, 'direction_descriptors': [(i, x, 0) for i in range(len(primitives)) for x in range(5)],
            'source_record': str(path), 'source_record_sha256': source_hash,
            'scope': 'Complementary Gaussian columns with only the same z^-8 monomial normalization as the reference; no gcd deflation or extra shifts.'}


def independent_layers(source, path, source_hash):
    directions = literal_directions(source)
    encoded = [sum(low << (19*entry) for entry, (low, high) in enumerate(direction)) for direction in directions]
    remaining = list(range(len(directions)))
    records = []
    while remaining:
        pivots = {}
        retained = []
        omitted = []
        for index in remaining:
            value = encoded[index]
            while value:
                pivot = value.bit_length()-1
                if pivot not in pivots:
                    pivots[pivot] = value
                    retained.append(index)
                    break
                value ^= pivots[pivot]
            else:
                omitted.append(index)
        if not retained:
            assert all(encoded[index] == 0 for index in omitted)
            break
        records.append({'constructor': 'aligned_gaussian_independent_generator_layer', 'dimensions': [5]*4,
                        'degree': 18, 'layer': len(records), 'primitives': source['primitives'],
                        'direction_descriptors': [source['direction_descriptors'][index] for index in retained],
                        'original_direction_indices': retained, 'source_record': str(path),
                        'source_record_sha256': source_hash,
                        'scope': 'An independent subset of original generator directions; all layers together cover every nonzero original direction, preserving its literal polynomial shape.'})
        remaining = omitted
    assert sorted(index for record in records for index in record['original_direction_indices']) == [i for i, value in enumerate(encoded) if value]
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--raw-output', type=Path)
    parser.add_argument('--layer-prefix', type=Path)
    args = parser.parse_args()
    source = json.loads(args.source.read_text())
    digest = sha256(args.source.read_bytes()).hexdigest()
    if args.raw_output:
        record = raw_family(source, args.source, digest)
        args.raw_output.write_text(json.dumps(record, indent=2)+'\n')
        print(json.dumps({'raw_directions': len(record['direction_descriptors'])}))
    if args.layer_prefix:
        records = independent_layers(source, args.source, digest)
        for index, record in enumerate(records):
            Path(str(args.layer_prefix)+str(index)+'.json').write_text(json.dumps(record, indent=2)+'\n')
        print(json.dumps({'independent_layer_sizes': [len(record['direction_descriptors']) for record in records]}))


if __name__ == '__main__':
    main()
