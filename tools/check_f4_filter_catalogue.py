#!/usr/bin/env python3
"""Replay every representative in a losslessly compressed F4 catalogue."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from compact_f4_filter_screen import compact, expand_rows
from f4_repaired_identity_filters import compiled_monomial_profiles


def verify(path):
    record = json.loads(path.read_text())
    roundtrip = compact(record)
    assert roundtrip['exact_catalogue_sha256'] == record['exact_catalogue_sha256']
    rows, cache = expand_rows(record), {}
    for row in rows:
        exponents = tuple(row['exponents'])
        if exponents not in cache:
            cache[exponents] = compiled_monomial_profiles(
                record['n'], record['ramification'], record['modes'], exponents, record['geometry'])
        actual = cache[exponents](tuple(row['coefficients']))
        assert [(block['charge'], list(block['orders'])) for block in actual] == [
            (block['charge'], block['orders']) for block in row['zero_charge_blocks_at_J1']]
    return {'status': 'pass', 'source': str(path), 'representative_profiles_replayed': len(rows),
            'exponent_tables': len(cache),
            'total_recorded_coefficient_multiplicity': sum(row['coefficient_multiplicity'] for row in rows),
            'exact_catalogue_sha256': record['exact_catalogue_sha256']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('record', type=Path, nargs='?',
                        default=Path('evidence/f4-repaired-three-mode-initial-n32-screen.json'))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify(args.record)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
