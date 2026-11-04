#!/usr/bin/env python3
"""Check the explicit d=8 Bell-diagonal path and its finite rank profile."""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import log, sqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bell_crossing import linear_profile_ternary, pencil_profiles, transform

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
    assert maximum < 1
    for item in profiles:
        exponent = log(item['max_rank']) / (item['t'] * log(64))
        item['max_rank_exponent'] = exponent
        assert exponent > 1 / sqrt(2)
    exponent_max = max(item['max_rank_exponent'] for item in profiles)
    exponent_min = min(item['max_rank_exponent'] for item in profiles)

    census = []
    for a in range(65):
        for b in range(65 - a):
            c = 64 - a - b
            profile = linear_profile_ternary((a, b, c), 2)
            census.append((max(profile), a, b, c))
    minimum = min(census)
    assert minimum[0] == 1848
    assert sum(1 for row in census if row[0] == 1848) == 2
    assert 1848 > 64 ** sqrt(2)

    monomial = []
    for n0 in range(65):
        for n1 in range(65 - n0):
            n2 = 64 - n0 - n1
            profile = (n0 * n0,
                       2 * n0 * n1,
                       n1 * n1 + 2 * n0 * n2,
                       2 * n1 * n2,
                       n2 * n2)
            monomial.append((max(profile), n0, n1, n2))
    monomial_minimum = min(monomial)
    assert monomial_minimum[0] == 1173
    assert monomial_minimum[1:] in ((13, 17, 34), (34, 17, 13))
    assert 1173 > 64 ** sqrt(2)

    receipt = {
        'passed': True,
        'd': 8,
        'q': 64,
        'prime': 3,
        'max_rank_fraction': maximum,
        'max_rank_fraction_t': next(item['t'] for item in profiles
                                    if item['max_ratio'] == maximum),
        'min_rank_exponent': exponent_min,
        'max_rank_exponent': exponent_max,
        'all_rank_exponents_above_inverse_sqrt_two': True,
        'not_a_width_crossing': True,
        'linear_t2_minimum_rank': 1848,
        'linear_t2_all_excluded': True,
        'monomial_t2_minimum_rank': 1173,
        'monomial_t2_all_excluded': True,
        'finite_lengths': MAX_LENGTH,
        'input_sha256': sha256(args.input.read_bytes()).hexdigest(),
        'scope': (
            'Exact finite profile of the explicit linear Bell-diagonal path '
            'over F3 at d=8 for tensor lengths 1 through 100. The rank '
            'fractions are below 1, but every rank exponent is above '
            '1/sqrt(2). The count census shows every linear ternary path has '
            'length-two rank at least 1848, above 64^sqrt(2); the whole linear '
            'path family is excluded at t=2. No score change.'
        ),
    }
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'passed': True, 'd': 8, 'max_rank_fraction': maximum,
                      'max_rank_exponent': exponent_max, 'linear_t2_min': 1848,
                      'linear_paths_excluded': True}))


if __name__ == '__main__':
    main()
