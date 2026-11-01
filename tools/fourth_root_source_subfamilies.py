#!/usr/bin/env python3
"""Export the seven new nonconstant fourth-root source primitives for screening."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


def export(source, source_path, digest, unshifted):
    assert source['dimensions'] == [5]*4 and source['degree'] == 18
    assert [row['normalized_degree'] for row in source['primitives']] == [1, 1, 1, 3, 3, 3, 4, 5, 5, 5]
    descriptors = [row for row in source['direction_descriptors']
                   if row[0] >= 3 and (not unshifted or row[2] == 0)]
    assert len(descriptors) == (35 if unshifted else 525)
    return {'constructor': 'nonconstant_fourth_root_source_primitive_subset',
            'dimensions': [5]*4, 'degree': 18, 'primitives': source['primitives'],
            'direction_descriptors': descriptors, 'source_record': str(source_path),
            'source_record_sha256': digest, 'primitive_indices': list(range(3, 10)),
            'unshifted_only': unshifted,
            'scope': 'The seven nonconstant endpoint-zero fourth-root primitives; the three previously screened constant supports are omitted.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('--unshifted-output', type=Path)
    parser.add_argument('--shifted-output', type=Path)
    args = parser.parse_args()
    source = json.loads(args.source.read_text())
    digest = sha256(args.source.read_bytes()).hexdigest()
    for path, unshifted in ((args.unshifted_output, True), (args.shifted_output, False)):
        if path:
            result = export(source, args.source, digest, unshifted)
            path.write_text(json.dumps(result, indent=2)+'\n')
            print(json.dumps({'output': str(path), 'directions': len(result['direction_descriptors'])}))


if __name__ == '__main__':
    main()
