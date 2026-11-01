#!/usr/bin/env python3
"""Exact affine closure, changed determinants, and finite local profiles."""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from physical_split_crossings import flatten
from reflection_dual_rail_gaussian import tensor
from reflection_laurent_physical_frames import add, scale
from search_affine_coupled_reflection_frames import affine_seed, catalogue
from search_reflection_physical_bases import packed_polynomial
from square_zero_coupled_reflection_frames import construct
from weighted_smith_oracle_f4 import WeightedSmithOracle


def multiply(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return value


def power(a, n):
    out = 1
    for _ in range(n):
        out = multiply(out, a)
    return out


def translate(value, point):
    lo, hi = value
    result = [0] * max(lo.bit_length(), hi.bit_length(), 1)
    for exponent in range(len(result)):
        coefficient = ((lo >> exponent) & 1) + 2 * ((hi >> exponent) & 1)
        subset = exponent
        while True:
            result[subset] ^= multiply(coefficient, power(point, exponent - subset))
            if subset == 0:
                break
            subset = (subset - 1) & exponent
    return (sum((v & 1) << i for i, v in enumerate(result)),
            sum((v >> 1) << i for i, v in enumerate(result)))


def rank(matrix):
    work = [list(row) for row in matrix]
    r = 0
    for j in range(len(work[0])):
        p = next((i for i in range(r, len(work)) if work[i][j]), None)
        if p is None:
            continue
        work[r], work[p] = work[p], work[r]
        inverse = power(work[r][j], 2)
        work[r] = [multiply(x, inverse) for x in work[r]]
        for i in range(r + 1, len(work)):
            c = work[i][j]
            work[i] = [x ^ multiply(c, y) for x, y in zip(work[i], work[r])]
        r += 1
    return r


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    base = tensor(4)
    word = [[56, 2], [58, 2]]
    frames = catalogue(base)
    expected_metadata = [{'target': 2, 'source': 1, 'scalar': c, 'pole': 1} for c in (2, 3)]
    assert [{k: v for k, v in frames[i].items() if k != 'delta'} for i, _ in word] == expected_metadata
    candidate = affine_seed(base, frames, word)
    independent = list(base['values'])
    for scalar in (2, 3):
        g = [[0] * 4 for _ in range(4)]
        g[2][1] = scalar
        transformed = construct(base, g, pole=1)
        independent = [add(current, scale(add(original, changed), 2))
                       for current, original, changed in zip(independent, base['values'], transformed['values'])]
    assert independent == candidate['values']
    assert candidate['degree'] == base['degree'] == 24
    for original, changed in zip(base['values'], candidate['values']):
        assert tuple(x >> 24 for x in original) == tuple(x >> 24 for x in changed)
    controls = []
    for point in (2, 3):
        record = {'point': point, 'base_ranks': [], 'changed_ranks': [], 'changed_profiles': []}
        for seed, field in ((base, 'base_ranks'), (candidate, 'changed_ranks')):
            shifted = dict(seed, values=[translate(v, point) for v in seed['values']])
            for cut in (0, 1):
                matrix, _, _ = flatten(shifted, cut)
                evaluated = [[(lo & 1) + 2 * (hi & 1) for lo, hi in row] for row in matrix]
                observed_rank = rank(evaluated)
                record[field].append(observed_rank)
                if field == 'changed_ranks':
                    oracle = WeightedSmithOracle([[packed_polynomial(v) for v in row] for row in matrix])
                    profile = oracle.profile([0] * 16, [0] * 16, Fraction(1))[0]
                    assert sum(x == 0 for x in profile) == observed_rank
                    assert all(x.denominator == 1 for x in profile)
                    record['changed_profiles'].append([int(x) for x in profile])
        assert record['base_ranks'] == [16, 16]
        assert record['changed_ranks'] == [16, 11]
        controls.append(record)
    names = (Path(__file__).name, 'search_affine_coupled_reflection_frames.py',
             'square_zero_coupled_reflection_frames.py', 'reflection_laurent_physical_frames.py',
             'weighted_smith_oracle_f4.py')
    result = {'status': 'passed', 'scope': 'exact affine strict path and determinant change; no pressure certificate',
              'word': word, 'frame_metadata': expected_metadata, 'degree': 24,
              'monomial_content': candidate['monomial_content'],
              'seed_sha256': sha256(json.dumps(candidate['values'], separators=(',', ':')).encode()).hexdigest(),
              'independent_coefficient_entries': 256, 'unchanged_leading_entries': 256,
              'endpoint_entries': 256, 'complete_finite_profiles': 4, 'controls': controls,
              'source_sha256': {n: sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in names}}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
