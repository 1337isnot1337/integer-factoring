#!/usr/bin/env python3
"""Lossless compact catalogue for the exact data in an F4 filter screen.

Repeated per-block metadata is shared, and each sorted Smith profile is
run-length encoded. Every exponent/color representative and coefficient
multiplicity is retained, together with its floating search observations.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path


ROW_FIELDS = ('exponents', 'coefficients', 'coefficient_multiplicity',
              'zero_charge_blocks_rle', 'critical_grid', 'relative_grid', 'refined')


def expand_rows(record):
    if record.get('catalogue_format') != 'f4-filter-screen-rle-v1':
        return record['rows']
    rows = []
    for values in record['catalogue']:
        row = dict(zip(ROW_FIELDS, values))
        blocks = row.pop('zero_charge_blocks_rle')
        row['zero_charge_blocks_at_J1'] = [
            {'charge': charge, 'orders': [order for order, count in runs for _ in range(count)]}
            for charge, runs in blocks]
        if row['refined'] is None:
            del row['refined']
        rows.append(row)
    return rows


def compact(record):
    rows = expand_rows(record)
    result = {key: value for key, value in record.items() if key not in ('rows', 'catalogue')}
    result['catalogue_format'] = 'f4-filter-screen-rle-v1'
    result['catalogue_row_fields'] = ROW_FIELDS
    result['catalogue'] = []
    exact = []
    for row in rows:
        runs = [[block['charge'], sorted(Counter(block['orders']).items())]
                for block in row['zero_charge_blocks_at_J1']]
        result['catalogue'].append([
            row['exponents'], row['coefficients'], row['coefficient_multiplicity'], runs,
            row['critical_grid'], row['relative_grid'], row.get('refined')])
        exact.append(result['catalogue'][-1][:4])
    result['exact_catalogue_sha256'] = sha256(json.dumps(exact, separators=(',', ':')).encode()).hexdigest()
    assert sum(row['coefficient_multiplicity'] for row in rows) == result['completed_coefficient_cases']
    restored = expand_rows(json.loads(json.dumps(result)))
    for old, new in zip(rows, restored):
        for key in ('exponents', 'coefficients', 'coefficient_multiplicity',
                    'critical_grid', 'relative_grid'):
            assert old[key] == new[key]
        assert [(block['charge'], list(block['orders'])) for block in old['zero_charge_blocks_at_J1']] == [
            (block['charge'], block['orders']) for block in new['zero_charge_blocks_at_J1']]
        assert old.get('refined') == new.get('refined')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.input.read_text())
    result = compact(source)
    encoded = json.dumps(result, separators=(',', ':')) + '\n'
    args.output.write_text(encoded)
    print(json.dumps({'distinct_profile_classes': len(result['catalogue']),
                      'coefficient_cases': result['completed_coefficient_cases'],
                      'output_bytes': len(encoded.encode()),
                      'exact_catalogue_sha256': result['exact_catalogue_sha256']}, indent=2))


if __name__ == '__main__':
    main()
