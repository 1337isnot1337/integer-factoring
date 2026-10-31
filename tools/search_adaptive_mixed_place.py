#!/usr/bin/env python3
"""Search adaptive opposite-flag rates and optional exact finite prefixes.

Only the finite rank costs are exact. Floating entropy grids select candidates
for the separate rational interval certificate generator.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import log
from pathlib import Path
import time

from adaptive_mixed_place import adaptive_rate, finite_rank_cost
from search_mixed_place_profiles import partitions


def save(path, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(result, indent=2) + '\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('profiles', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--alphabet', type=int)
    parser.add_argument('--entropy-steps', type=int, default=8193)
    parser.add_argument('--rate-steps', type=int, default=4097)
    parser.add_argument('--copies', type=int, nargs='*', default=[])
    args = parser.parse_args()
    raw = args.profiles.read_bytes()
    records = json.loads(raw)
    if isinstance(records, dict):
        records = records['cases']
    if args.alphabet:
        records = [r for r in records if r['alphabet'] == args.alphabet]
    output = []
    started = time.monotonic()
    for record in records:
        candidates = []
        choices = list(partitions(record))
        for row, col in choices:
            result = adaptive_rate([p['cuts'][0] for p in row], [p['cuts'][1] for p in col],
                                   record['degree'], record['alphabet'],
                                   args.entropy_steps, args.rate_steps)
            result['row_places'] = [p['point'] for p in row]
            result['column_places'] = [p['point'] for p in col]
            candidates.append(result)
        index = min(range(len(candidates)), key=lambda i: candidates[i]['rho'])
        best, (row, col) = candidates[index], choices[index]
        entry = {key: record[key] for key in ('m', 'theta', 'alphabet', 'degree') if key in record}
        entry.update({'source_profiles_sha256': sha256(raw).hexdigest(),
                      'partitions': len(candidates), 'best': best, 'finite_checks': []})
        output.append(entry)
        save(args.output, output)
        print(json.dumps({'alphabet': record['alphabet'], 'best_numeric_rho': best['rho'],
                          'partitions': len(candidates)}), flush=True)
        for copies in args.copies:
            result = finite_rank_cost([p['cuts'][0] for p in row], [p['cuts'][1] for p in col],
                                      copies, record['degree'])
            product = result['product']
            finite = {'copies': copies, 'critical_index': result['index'],
                      'exact_max_product_hex': hex(product),
                      'normalized_log_product': log(product) / (4 * copies * log(record['alphabet'])),
                      'scope': 'exact greedy-flag bound for the displayed selected partition'}
            entry['finite_checks'].append(finite)
            entry['elapsed_seconds'] = time.monotonic() - started
            save(args.output, output)
            print(json.dumps({key: finite[key] for key in ('copies', 'critical_index',
                                                          'normalized_log_product')}), flush=True)


if __name__ == '__main__':
    main()
