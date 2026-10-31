#!/usr/bin/env python3
"""Replay every distinct saved F4 operator profile from its Laurent kernel.

This bypasses the affine bit-plane compiler used during enumeration.
It verifies catalogue counts, distinctness, endpoints, and every complete
two-cut profile representative. It does not prove global gauge optimality.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from compiled_f4_identity_operators import combine_kernel
from coupled_identity_operator_filters import identity_kernel
from f4_identity_operator_filters import zero_profiles, infinity_profiles
from search_f4_identity_operators import endpoint_blocks
from search_split_f4_identity_operators import affine_sources, sector_pairs


ROOT = Path(__file__).resolve().parents[1]


def canonical(blocks, infinity):
    return tuple((block['charge'], tuple(block['orders'])) for block in blocks), tuple(infinity)


def two_mode_kernel(run, row):
    middle = ((1, 0), (0, 1))
    constant = identity_kernel(2, 2)
    constant[(0, 0), (0, 0)] = {0: run['vacuum_color']}
    constant[(1, 1), (1, 1)] = {-run['top_exponent']: 1}
    for i, a in enumerate(middle):
        for j, b in enumerate(middle):
            if row['endpoint'][i][j]:
                constant[a, b] = {0: row['endpoint'][i][j]}
            else:
                constant.pop((a, b), None)
    perturbations = [{(a, b): {0: 1, -degree: 1}}
                     for degree in range(1, run['slow_degree']+1) for a, b in product(middle, repeat=2)]
    return combine_kernel(constant, perturbations, row['colors'])


def replay(run, split=False):
    n, ramification, geometry = run['n'], run['ramification'], run['geometry']
    modes = 3 if split else 2
    rows = run['profile_classes']
    expected = len(sector_pairs())*16 if split else len(endpoint_blocks(run['endpoint_policy']))*4**(4*run['slow_degree'])
    assert run['complete'] and run['completed_candidates'] == expected
    assert run['profile_class_count'] == len(rows)
    assert all(isinstance(row['rows'], int) and row['rows'] > 0 for row in rows)
    assert sum(row['rows'] for row in rows) == expected
    keys = set()
    if split:
        constant, perturbations = affine_sources(n, ramification, geometry, run['one_height'], run['two_height'], run['top_height'])
    for row in rows:
        key = canonical(row['zero_blocks'], row['infinity_active'])
        assert key not in keys
        keys.add(key)
        kernel = combine_kernel(constant, perturbations, row['colors']) if split else two_mode_kernel(run, row)
        actual = canonical(zero_profiles(n, ramification, modes, kernel, geometry),
                           infinity_profiles(n, ramification, modes, kernel, geometry))
        assert actual == key
    return {'n': n, 'ramification': ramification, 'geometry': geometry,
            'modes': modes, 'completed_catalogue_candidates': expected,
            'all_distinct_profile_representatives_replayed': len(rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    paths = [ROOT/'evidence'/name for name in ('f4-operator-degree1-mixed-n32.json',
             'f4-operator-degree2-mixed-n32.json', 'split-f4-operator-projectors.json')]
    records = [replay(json.loads(path.read_text())) for path in paths[:2]]
    for run in json.loads(paths[2].read_text())['runs']:
        records.append(replay(run, True))
    result = {'status': 'All distinct saved operator profiles independently reconstructed and replayed.',
              'runs': records, 'total_distinct_profiles_replayed': sum(row['all_distinct_profile_representatives_replayed'] for row in records),
              'source_sha256': {str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
                                for path in paths+[Path(__file__)]}}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
