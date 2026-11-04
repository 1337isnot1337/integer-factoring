#!/usr/bin/env python3
"""Check the explicit d=8 Bell-diagonal path and its finite rank profile."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import sqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bell_crossing import pencil_profiles, transform

INPUT = Path('evidence/bell-d8-explicit-family.json')
MAX_LENGTH = 100


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=INPUT)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()

    saved = json.loads(args.input.read_text())
    assert saved['completed'] is True
    assert saved['bits'] == 3 and saved['d'] == 8 and saved['q'] == 64
    assert saved['prime'] == 3
    values = tuple(saved['values'])
    assert len(values) == 64
    assert set(values) <= {0, 1, 2}
    assert saved['value_counts'] == {str(k): v for k, v in Counter(values).items()}

    other = transform(values, 3, 3)
    assert tuple(saved['transform_values']) == other
    assert saved['transform_counts'] == {str(k): v for k, v in Counter(other).items()}
    assert transform(other, 3, 3) == values
    one = (1,) * 64
    assert transform(one, 3, 3) == one

    profiles = []
    for length in range(1, MAX_LENGTH + 1):
        row, column = pencil_profiles(values, 3, length)
        maximum = max(max(row), max(column))
        profiles.append({
            't': length,
            'row': list(row),
            'col': list(column),
            'max_rank': maximum,
            'max_ratio': maximum / 64 ** length,
        })
    assert profiles == saved['profiles']
    assert saved['finite_lengths'] == MAX_LENGTH
    maximum = max(item['max_ratio'] for item in profiles)
    assert abs(maximum - saved['max_ratio']) < 1e-18
    assert saved['max_ratio_fraction'] == '43/64'
    assert abs(maximum - 43 / 64) < 1e-18
    assert maximum < 1 / sqrt(2)
    assert all(item['max_ratio'] < 1 / sqrt(2) for item in profiles)

    receipt = {
        'passed': True,
        'd': 8,
        'q': 64,
        'prime': 3,
        'max_ratio': maximum,
        'max_ratio_t': next(item['t'] for item in profiles
                            if item['max_ratio'] == maximum),
        'below_inverse_sqrt_two': True,
        'finite_lengths': MAX_LENGTH,
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'scope': (
            'Exact finite profile of the explicit linear Bell-diagonal path '
            'over F3 at d=8 for tensor lengths 1 through 12. The maximum rank '
            'ratio is 43/64 < 1/sqrt(2). This is not an all-length uniform '
            'two-cut transfer and does not change the factoring exponent.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'd': 8, 'max_ratio': maximum,
                      'max_ratio_t': receipt['max_ratio_t'], 'below': True}))


if __name__ == '__main__':
    main()
