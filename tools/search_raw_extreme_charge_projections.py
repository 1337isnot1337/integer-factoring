#!/usr/bin/env python3
"""All affine charge classes on the raw Gaussian's two extreme modes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

from strict_core_charge_profiles import affine_charge_representatives, optimize_bulk
from strict_core_charge_gauges import CoreChargePressure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--interval', type=int, default=8)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    records = []
    started = time.monotonic()
    for charges in affine_charge_representatives(2):
        model = CoreChargePressure(args.interval, charges, family='raw_extremes')
        profiles = []
        for record, oracle in zip(model.records, model.oracles):
            data, metadata, left, right = record
            result = oracle.profile()
            if not metadata['exact_polynomial_matrix'] and max(result['orders']) >= metadata['known_below']:
                raise ArithmeticError('The complete profile exceeds the certified coefficient window.')
            profiles.append({'place': data['place'], 'active_orders': result['orders'],
                'active_charges': [left['row_charges'][i] for i in result['rows']],
                'spectator_orders': data['spectator_orders'],
                'spectator_charges': data['spectator_charges'],
                'active_rank': len(result['orders']), 'profile_precision': metadata['profile_precision']})
        record = {'charges': list(charges), 'profiles': profiles, **optimize_bulk(args.interval, profiles)}
        records.append(record)
        args.output.write_text(json.dumps({'status': 'exact active profiles; exploratory floating bulk-gauge roots',
            'family': 'raw Gaussian with typed endpoint and extreme-core charge mask',
            'interval': args.interval, 'physical_core_indices': [0, args.interval - 1],
            'seconds': time.monotonic() - started, 'records': records}, indent=2) + '\n')
        print(json.dumps({key: value for key, value in record.items() if key != 'profiles'}), flush=True)


if __name__ == '__main__':
    main()
