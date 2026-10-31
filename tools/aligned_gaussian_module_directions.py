"""Complete polynomial saturation of the actual q5 Gaussian complements."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

from source_slab_module_saturation import direction_record


def aligned_module_record(source, degree=18):
    slabs = [item['raw_slab_in_z'] for item in source['raw_slabs']]
    result = direction_record({'dimensions': source['dimensions'], 'kernel_slabs': slabs}, degree)
    result.update(constructor='aligned_gaussian_complement_module_directions',
                  reference_digest=source['reference_digest'],
                  reference_identity=source['reference_identity'],
                  complementary_masks=[item['complement_mask'] for item in source['raw_slabs']],
                  scope='Complete polynomial rational-span saturation of literal complementary slabs of the retained five-letter Gaussian; endpoint factors are restored exactly.')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--degree', type=int, default=18)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = aligned_module_record(json.loads(args.source.read_text()), args.degree)
    result['source_record_sha256'] = sha256(args.source.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'rank': result['rational_row_certificate']['rank'],
                      'minimal_degrees': result['minimal_basis']['degrees'],
                      'directions': len(result['direction_descriptors']),
                      'saturation': result['saturation_certificate']}))


if __name__ == '__main__':
    main()
