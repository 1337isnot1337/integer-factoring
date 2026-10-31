#!/usr/bin/env python3
"""Losslessly compact repeated exact profiles in a polynomial lattice screen."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path


ROW_FIELDS = ('degree_bounds', 'row_targets', 'profile_index', 'endpoint_values',
              'endpoint_image_dimension', 'critical_bound_grid', 'relative_bound_grid',
              'witness_polynomial_sha256')


def compact(source):
    source = json.loads(json.dumps(source))
    result = {key: value for key, value in source.items() if key != 'rows'}
    profiles, lookup, rows = [], {}, []
    for row in source['rows']:
        profile = [row['actual_degrees'],
                   [[block['charge'], block['orders']] for block in row['exact_zero_blocks_at_J1']],
                   row['actual_critical_grid'], row['actual_relative_grid']]
        key = json.dumps(profile, separators=(',', ':'))
        if key not in lookup:
            lookup[key] = len(profiles)
            profiles.append(profile)
        record = dict(row)
        record['profile_index'] = lookup[key]
        if set(record) - set(ROW_FIELDS) - {'actual_degrees', 'exact_zero_blocks_at_J1',
                                          'actual_critical_grid', 'actual_relative_grid'}:
            raise ValueError('This lossless format must be extended before compacting extra row fields.')
        rows.append([record[field] for field in ROW_FIELDS])
    result.update({'catalogue_format': 'joint-lattice-screen-rle-v1',
                   'profile_fields': ['actual_degrees', 'zero_blocks', 'actual_critical_grid', 'actual_relative_grid'],
                   'row_fields': list(ROW_FIELDS), 'profile_table': profiles, 'rows': rows,
                   'distinct_actual_profile_pairs': len(profiles),
                   'exact_expanded_rows_sha256': sha256(json.dumps(source['rows'], separators=(',', ':')).encode()).hexdigest()})
    assert expand_rows(result) == source['rows']
    return result


def expand_rows(record):
    if record.get('catalogue_format') != 'joint-lattice-screen-rle-v1':
        return record['rows']
    result = []
    for values in record['rows']:
        row = dict(zip(record['row_fields'], values))
        profile = record['profile_table'][row.pop('profile_index')]
        row['actual_degrees'] = profile[0]
        row['exact_zero_blocks_at_J1'] = [{'charge': charge, 'orders': orders} for charge, orders in profile[1]]
        row['actual_critical_grid'] = profile[2]
        row['actual_relative_grid'] = profile[3]
        # Preserve the original generator's field insertion order, so its
        # canonical row digest can be independently reproduced.
        result.append({key: row[key] for key in (
            'degree_bounds', 'row_targets', 'critical_bound_grid', 'relative_bound_grid',
            'actual_degrees', 'endpoint_values', 'endpoint_image_dimension',
            'actual_critical_grid', 'actual_relative_grid', 'exact_zero_blocks_at_J1',
            'witness_polynomial_sha256')})
    assert sha256(json.dumps(result, separators=(',', ':')).encode()).hexdigest() == record['exact_expanded_rows_sha256']
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    value = compact(json.loads(args.source.read_text()))
    args.output.write_text(json.dumps(value, separators=(',', ':')) + '\n')
    print(json.dumps({'rows': len(value['rows']), 'profiles': len(value['profile_table']),
                      'bytes': args.output.stat().st_size,
                      'exact_expanded_rows_sha256': value['exact_expanded_rows_sha256']}))


if __name__ == '__main__':
    main()
