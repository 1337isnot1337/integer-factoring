#!/usr/bin/env python3
"""Finite F4 diagonal-filter screen on a proved mixed-end three-mode core."""
from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
from time import monotonic

from compact_f4_filter_screen import compact
from f4_repaired_identity_filters import compiled_monomial_profiles
from search_f4_repaired_identity_filters import candidate_exponents, pressure_grid


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=32)
    parser.add_argument('--ramification', type=int, default=1)
    parser.add_argument('--geometry', choices=('split:1:2', 'split:2:1'), default='split:1:2')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--raw-output', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    n, ramification, geometry = args.n, args.ramification, args.geometry
    tables = candidate_exponents(n, ramification, 3)
    if args.limit:
        tables = tables[:args.limit]
    colors = tuple((1,) + values for values in product((1, 2, 3), repeat=7))
    scales = sorted({value / n for value in (.1, .3, 1, 2, 3, 5, 8)} | {.125, .25, .5, 1., 2., 4.})
    charges = sorted({1 - 2 * n, -n, -3, -1, 0, 1, 2, 3, n / 2, n, 2 * n})
    evaluate = pressure_grid(n, ramification, geometry, scales, charges, modes=3)
    rows, seen, per_table = [], {}, []
    count, start = 0, monotonic()
    for index, exponents in enumerate(tables):
        before = len(rows)
        profile = compiled_monomial_profiles(n, ramification, 3, exponents, geometry)
        for coefficients in colors:
            blocks = profile(coefficients)
            count += 1
            key = exponents, tuple((block['charge'], block['orders']) for block in blocks)
            if key in seen:
                rows[seen[key]]['coefficient_multiplicity'] += 1
                continue
            seen[key] = len(rows)
            critical, relative = evaluate(exponents, blocks)
            row = {'exponents': exponents, 'coefficients': coefficients,
                   'zero_charge_blocks_at_J1': blocks, 'coefficient_multiplicity': 1,
                   'critical_grid': critical, 'relative_grid': relative}
            rows.append(row)
            if critical[0] < -1e-9 or relative[0] < -1e-9:
                print('CANDIDATE_GAIN', row, flush=True)
        per_table.append({'exponents': exponents, 'coefficient_cases': len(colors),
                          'distinct_profile_classes': len(rows) - before})
        if (index + 1) % 5 == 0:
            print('PROGRESS', index + 1, len(tables), count, len(rows), monotonic() - start, flush=True)
    result = {
        'status': 'Finite floating screen of exact split-core profiles; no certified pressure witness.',
        'characteristic': 2, 'coefficient_field': 'F4', 'modes': 3, 'occupation_states': 8,
        'n': n, 'ramification': ramification, 'geometry': geometry,
        'rho': '1/sqrt(2)', 'evaluation_scale': 't=P*log(s)',
        'exponent_vectors': len(tables), 'nonzero_coefficient_vectors_per_exponent': len(colors),
        'completed_coefficient_cases': count, 'distinct_profile_count': len(rows),
        'scale_grid': scales, 'charge_grid': charges, 'gain_detection_tolerance': 1e-9,
        'minimum_sampled_critical_pressure': min(row['critical_grid'][0] for row in rows),
        'minimum_sampled_relative_log_cost': min(row['relative_grid'][0] for row in rows),
        'same_gauge_gain_below_tolerance': any(row['relative_grid'][0] < -1e-9 for row in rows),
        'negative_critical_pressure_found': any(row['critical_grid'][0] < -1e-9 for row in rows),
        'elapsed_seconds': monotonic() - start, 'per_exponent_table': per_table, 'rows': rows,
    }
    plain = json.loads(json.dumps(result))
    if args.raw_output:
        args.raw_output.write_text(json.dumps(plain, separators=(',', ':')) + '\n')
    portable = compact(plain)
    encoded = json.dumps(portable, separators=(',', ':')) + '\n'
    if len(encoded.encode()) > 2_000_000:
        keep = set(sorted(range(len(rows)), key=lambda i: rows[i]['critical_grid'][0])[:100])
        keep.update(sorted(range(len(rows)), key=lambda i: rows[i]['relative_grid'][0])[:100])
        for exponents in tables:
            group = [i for i, row in enumerate(rows) if row['exponents'] == exponents]
            keep.add(min(group, key=lambda i: rows[i]['critical_grid'][0]))
        catalogue = portable.pop('catalogue')
        portable['catalogue_format'] = 'f4-filter-screen-selected-v1'
        portable['record_scope'] = 'Selected pressure representatives, full aggregate/table counts, and full exact-catalogue digest.'
        portable['representative_multiplicities_are_not_total_counts'] = True
        portable['selected_catalogue'] = [catalogue[i] for i in sorted(keep)]
        portable['selected_profile_classes'] = len(keep)
        encoded = json.dumps(portable, separators=(',', ':')) + '\n'
    args.output.write_text(encoded)
    print('DONE', count, len(rows), len(encoded.encode()), result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
