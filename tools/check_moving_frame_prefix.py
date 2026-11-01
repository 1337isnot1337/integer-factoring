#!/usr/bin/env python3
"""Exact all-scale natural-pressure exclusion for a frozen finite prefix.

This checker uses only the complete integer profiles, not their fitted
pressures. It does not exclude other flows, fields, gauges or tensors.
"""

from __future__ import annotations

import argparse
from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def product_profile(sides):
    left, right = map(Counter, sides)
    result = Counter()
    for a, ca in left.items():
        for b, cb in right.items():
            result[a + b] += ca * cb
    return result


def dominance(profile, reference):
    """Nonnegative cumulative coefficients certify a Laurent inequality."""
    cumulative, minimum = 0, None
    for exponent in sorted(profile.keys() | reference.keys()):
        cumulative += profile[exponent] - reference[exponent]
        minimum = cumulative if minimum is None else min(minimum, cumulative)
        if cumulative < 0:
            return False, minimum
    assert cumulative == 64**2 - 49**2
    return True, minimum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest_path = ROOT / 'evidence/moving-frame-height1-prefix16384.json'
    manifest = json.loads(manifest_path.read_text())
    archive = manifest_path.with_name(manifest['archive'])
    compressed = archive.read_bytes()
    assert sha256(compressed).hexdigest() == manifest['archive_sha256']
    raw = gzip.decompress(compressed)
    assert sha256(raw).hexdigest() == manifest['jsonl_sha256']
    core_path = ROOT / 'evidence/moving-boundary-core-floor.json'
    assert sha256(core_path.read_bytes()).hexdigest() == manifest['core_certificate_sha256']
    core = json.loads(core_path.read_text())
    assert core['status'] == 'passed' and core['rho_floor'] == [5, 7]
    reference = product_profile([[2 * value for value in core['core_orders_t']]] * 2)
    assert len(core['core_orders_t']) == 49
    lines = raw.splitlines()
    assert len(lines) == manifest['checked_prefix_cases'] == 16384
    assert len(manifest['a_prefix']) == len(set(manifest['a_prefix'])) == 64
    assert sorted(manifest['b_order']) == list(range(256))
    minimum = None
    for index, line in enumerate(lines):
        record = json.loads(line)
        assert record['a_mask'] == manifest['a_prefix'][index // 256]
        assert record['b_mask'] == manifest['b_order'][index % 256]
        sides = record['orders_infinity_zero']
        assert len(sides) == 2 and all(len(side) == 64 for side in sides)
        assert all(type(value) is int for side in sides for value in side)
        passed, low = dominance(product_profile(sides), reference)
        assert passed, (index, record['a_mask'], record['b_mask'], low)
        minimum = low if minimum is None else min(minimum, low)
    report = {'scope': __doc__, 'status': 'passed', 'checked_cases': len(lines),
              'full_planned_cases': manifest['full_planned_cases'],
              'full_search_exhausted': False, 'rho_natural_floor': [5, 7],
              'minimum_cumulative_coefficient': minimum,
              'final_cumulative_excess_each_case': 64**2 - 49**2,
              'numerical_fits_used': False,
              'all_common_positive_scales_and_physical_types': True,
              'manifest_sha256': sha256(manifest_path.read_bytes()).hexdigest(),
              'archive_sha256': manifest['archive_sha256'],
              'source_commit_of_profile_run': manifest['source_commit'],
              'source_sha256': {str(Path(__file__).relative_to(ROOT)): sha256(Path(__file__).read_bytes()).hexdigest()}}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key: report[key] for key in ('status', 'checked_cases', 'rho_natural_floor', 'minimum_cumulative_coefficient', 'full_search_exhausted')}))


if __name__ == '__main__':
    main()
