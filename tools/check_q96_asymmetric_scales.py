#!/usr/bin/env python3
"""Replay the q96 unequal-cut-scale screen and verify its scope guards."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import math

from screen_q96_asymmetric_scales import screen


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path,
                        default=Path('evidence/q96-asymmetric-scale-screen.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    saved = json.loads(args.input.read_text())
    assert saved['completed'] is True
    fresh = screen()
    assert saved['grid'] == fresh['grid']
    assert saved['interior_best']['name'] == fresh['interior_best']['name']
    for old, new in zip(saved['profiles'], fresh['profiles'], strict=True):
        assert old['name'] == new['name']
        assert old['profile_path'] == new['profile_path']
        assert old['profile_sha256'] == new['profile_sha256']
        assert old['channel_counts'] == new['channel_counts']
        assert old['height_count'] == new['height_count']
        assert old['expected_rho'] == new['expected_rho']
        assert abs(old['balanced']['support'] - new['balanced']['support']) < 1e-10
        assert abs(old['balanced']['rho'] - new['balanced']['rho']) < 1e-10
        assert abs(old['interior']['support'] - new['interior']['support']) < 1e-10
        assert abs(old['interior']['angle'] - new['interior']['angle']) < 1e-10
        assert abs(old['interior']['a'] - new['interior']['a']) < 1e-10
        assert abs(old['interior']['b'] - new['interior']['b']) < 1e-10
        assert abs(old['interior']['scale'] - new['interior']['scale']) < 1e-12
        for old_edge, new_edge in zip(old['boundary'], new['boundary'], strict=True):
            assert old_edge['angle'] == new_edge['angle']
            assert abs(old_edge['support'] - new_edge['support']) < 1e-10
    assert fresh['interior_best']['support'] > 1.0
    assert all(row['interior']['support'] > 1.0 for row in fresh['profiles'])
    assert all(row['boundary'][0]['support'] < 1.001 for row in fresh['profiles'])
    assert all(
        row['calibration_error'] is None or row['calibration_error'] < 1e-9
        for row in fresh['profiles']
    )
    source = Path(__file__).with_name('screen_q96_asymmetric_scales.py')
    receipt = {
        'passed': True,
        'scope': (
            'Replay of the deterministic q96 unequal-cut-scale screen. '
            'The numerical minimum is not a global optimum or an exact '
            'pressure certificate; the factoring exponent is unchanged.'
        ),
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'screen_source_sha256': sha256(source.read_bytes()).hexdigest(),
        'profiles': len(fresh['profiles']),
        'interior_best': fresh['interior_best'],
        'balanced_calibration_errors': [
            row['calibration_error'] for row in fresh['profiles']
            if row['calibration_error'] is not None
        ],
        'all_interior_supports_above_one': True,
        'boundary_layer_approaches_one': True,
        'elapsed_seconds': fresh['elapsed_seconds'],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print('q96 unequal-cut-scale replay passed')


if __name__ == '__main__':
    main()
