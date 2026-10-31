#!/usr/bin/env python3
"""Finite raw extreme-core pair-mask pressure search with exact profiles."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

from raw_core_mask_profiles import CoreMaskPressure, additive_pair_partitions, pair_partition_mask
from strict_core_charge_profiles import optimize_bulk


def candidates(include_binary=True):
    found = additive_pair_partitions()
    ordered = [((0,) * 10, {'family': 'Gaussian'}),
               (tuple(range(10)), {'family': 'all unordered pairs separated'})]
    ordered += [(key, {'family': 'finite additive charges', **value}) for key, value in found.items()]
    if include_binary:
        ordered += [((0,) + tuple((mask >> i) & 1 for i in range(9)),
                     {'family': 'binary partition'}) for mask in range(1, 512)]
    seen = set()
    for partition, detail in ordered:
        labels = {}
        partition = tuple(labels.setdefault(value, len(labels)) for value in partition)
        if partition not in seen:
            seen.add(partition)
            yield partition, detail


def run(interval, output, include_binary=True, limit=None):
    started = time.monotonic()
    records, best, best_nonconstant = [], None, None
    digest = sha256()
    all_candidates = list(candidates(include_binary))
    selected = all_candidates if limit is None else all_candidates[:limit]
    for partition, label in selected:
        model = CoreMaskPressure(interval, pair_partition_mask(partition))
        profiles = []
        for record, oracle in zip(model.records, model.oracles):
            data, metadata, left, right = record
            result = oracle.profile()
            if not metadata['exact_polynomial_matrix'] and max(result['orders']) >= metadata['known_below']:
                raise ArithmeticError('The exact profile exceeds the certified coefficient window.')
            profiles.append({'place': data['place'], 'active_orders': result['orders'],
                'active_charges': [left['row_charges'][i] for i in result['rows']],
                'spectator_orders': data['spectator_orders'],
                'spectator_charges': data['spectator_charges'],
                'active_rank': len(result['orders']), 'profile_precision': metadata['profile_precision'],
                'block_sizes': oracle.block_sizes})
        digest.update(bytes(partition))
        digest.update(json.dumps(profiles, sort_keys=True, separators=(',', ':')).encode())
        record = {'partition': partition, **label, **optimize_bulk(interval, profiles)}
        records.append(record)
        if best is None or record['rho'] < best['rho']:
            best = {**record, 'profiles': profiles}
        if max(partition) and (best_nonconstant is None or record['rho'] < best_nonconstant['rho']):
            best_nonconstant = {**record, 'profiles': profiles}
        output.write_text(json.dumps({'status': 'Exact active profiles; exploratory floating pressure roots.',
            'interval': interval, 'physical_core_indices': [0, interval - 1],
            'seconds': time.monotonic() - started, 'candidate_count': len(all_candidates),
            'selected_count': len(selected), 'completed_count': len(records),
            'complete_selected_screen': len(records) == len(selected),
            'profile_digest_sha256': digest.hexdigest(), 'best': best,
            'best_nonconstant': best_nonconstant, 'records': records}, indent=2) + '\n')
        print(json.dumps({'candidate': len(records), **record,
                          'best_nonconstant_rho': None if best_nonconstant is None else best_nonconstant['rho'],
                          'seconds': time.monotonic() - started}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=8)
    parser.add_argument('--additive-only', action='store_true')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.interval, args.output, not args.additive_only, args.limit)


if __name__ == '__main__':
    main()
