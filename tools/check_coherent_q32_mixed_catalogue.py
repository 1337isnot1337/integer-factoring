#!/usr/bin/env python3
"""Exact all-scale optimum across all 496 binary height-boundary maps.

Every saved profile is complete. The default verifies exact map recipes,
source hashes, and integer polynomial certificates. The optional native
replay recomputes the complete profiles rather than reusing the receipts.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import time

from check_coherent_q32_catalogue import certificate
from coherent_hot_mixed import profiles, validate_maps
from search_coherent_mixed_boundary import boundary_maps, catalogue_keys


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path,
                        default=ROOT / 'evidence/coherent-q32-mixed-boundary.json')
    parser.add_argument('--recompute-profiles', action='store_true')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    saved = json.loads(args.catalogue.read_text())
    assert saved['completed'] and saved['finite_list_exhausted'] and saved['unresolved'] == 0
    assert len(saved['records']) == saved['total_maps'] == 496
    for name, digest in saved['source_sha256'].items():
        assert sha256((ROOT / 'tools' / name).read_bytes()).hexdigest() == digest, name
    winner = next(record for record in saved['records']
                  if record['description']['kernel'] == 7 and record['description']['covector'] == 15)
    assert winner['completed']
    keys = catalogue_keys()
    assert len(keys) == len(set(keys)) == 496
    counts, records = Counter(), []
    start = time.monotonic()
    for key, record in zip(keys, saved['records']):
        description, embedding, retraction = boundary_maps(*key)
        assert record['completed'] and record['description'] == description
        assert record['embedding'] == embedding and record['retraction'] == retraction
        _, _, heights, checks = validate_maps(embedding, retraction, 8)
        assert list(heights) == record['heights']
        assert sorted(heights) == sorted(winner['heights'])
        counts['deletion_checks'] += checks
        for name in ('zero_orders_t', 'infinity_orders_t'):
            assert len(record[name]) == 1024 and record[name] == sorted(record[name])
            assert all(isinstance(order, int) and order >= 0 for order in record[name])
            counts['complete_saved_profiles'] += 1
            counts['complete_saved_smith_orders'] += len(record[name])
        counts['reused_occupation_maps' if record['reused_occupation_profile'] else 'new_mixed_maps'] += 1
        if args.recompute_profiles:
            current = profiles(embedding, retraction, 8)
            for name in ('zero_orders_t', 'infinity_orders_t'):
                assert current[name] == record[name]
                counts['complete_profiles_recomputed'] += 1
        proof = certificate(record, winner)
        counts[proof['method']] += 1
        records.append({'kernel': key[0], 'covector': key[1], **proof})
    assert counts['reused_occupation_maps'] == 121 and counts['new_mixed_maps'] == 375
    assert counts['nonnegative_power_basis'] == 326
    assert counts['nonnegative_unnormalized_bernstein_basis'] == 168
    assert counts['identical_product_polynomials'] == 2
    ties = [[record['kernel'], record['covector']] for record in records
            if record['method'] == 'identical_product_polynomials']
    assert sorted(ties) == [[7, 15], [15, 7]]
    names = (Path(__file__).name, 'check_coherent_q32_catalogue.py',
             'check_coherent_q32_pressure.py', 'coherent_hot_mixed.py',
             'coherent_hot_mixed.cpp', 'coherent_hot_compounds.py',
             'coherent_hot_compounds_fast.cpp', 'search_coherent_mixed_boundary.py',
             'selected_gaussian_coefficients.cpp', 'native_rectangular_smith_f4.cpp')
    report = {'status': 'pass', 'scope': __doc__, 'counts': dict(counts),
              'records': records, 'optimal_keys': ties, 'optimality_interval': '0 < r < 1',
              'all_496_maps_compared': True, 'all_scale_catalogue_optimality_proved': True,
              'physical_height_multisets_identical': True,
              'all_saved_profiles_recomputed': args.recompute_profiles,
              'floating_point_pressure_used': False, 'factoring_exponent_improved': False,
              'catalogue_sha256': sha256(args.catalogue.read_bytes()).hexdigest(),
              'source_sha256': {name: sha256((ROOT / 'tools' / name).read_bytes()).hexdigest()
                                for name in names},
              'elapsed_seconds': time.monotonic() - start}
    args.output.write_text(json.dumps(report, separators=(',', ':')) + '\n')
    print(json.dumps({'status': 'pass', 'counts': dict(counts), 'optimal_keys': ties,
                      'elapsed_seconds': report['elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()
