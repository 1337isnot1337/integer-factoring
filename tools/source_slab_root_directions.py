"""Literal Frobenius-root compression of saturated source directions.

This operation creates new endpoint-zero polynomial directions. It does
not assert that a root remains in the original source rational module.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from source_slab_module_saturation import direction_record


def compress_exponents(value, factor):
    if factor < 1:
        raise ValueError('The exponent compression factor must be positive.')
    result = 0
    while value:
        bit = value & -value
        value ^= bit
        exponent = bit.bit_length() - 1
        if exponent % factor:
            raise ValueError('A nonzero exponent is not divisible by the factor.')
        result ^= 1 << (exponent // factor)
    return result


def inflate_exponents(value, factor):
    result = 0
    while value:
        bit = value & -value
        value ^= bit
        result ^= 1 << ((bit.bit_length() - 1) * factor)
    return result


def root_direction_record(source, factor=4, degree=18):
    rows = source['minimal_basis']['generators']
    compressed = [[compress_exponents(value, factor) for value in row] for row in rows]
    result = direction_record({'dimensions': source['dimensions'],
                               'kernel_slabs': compressed}, degree)
    result.update(constructor='source_slab_literal_root_directions',
                  exponent_compression_factor=factor,
                  compressed_source_generators=compressed,
                  original_generator_degrees=source['minimal_basis']['degrees'],
                  scope='Exact new endpoint-zero polynomial family from literal exponent compression; no membership in the original source module is asserted.')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--factor', type=int, default=4)
    parser.add_argument('--degree', type=int, default=18)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    record = root_direction_record(json.loads(args.source.read_text()), args.factor, args.degree)
    record['source_record_sha256'] = sha256(args.source.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'minimal_degrees': record['minimal_basis']['degrees'],
                      'endpoint_degrees': [item['normalized_degree'] for item in record['primitives']],
                      'directions': len(record['direction_descriptors'])}))


if __name__ == '__main__':
    main()
