#!/usr/bin/env python3
"""Exact four-jet exclusions for every normalized F4 Weyl boundary."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random

from weyl_ordered_crossing import (MULTIPLICATION, ROOT_POWERS, boundary_jets,
                                    canonical_weights, coefficient_cuts,
                                    field_rank, scalar_restriction, toeplitz_rank)

ROOT = Path(__file__).resolve().parents[1]


def dense_product(first, second):
    answer = [[0]*3 for _ in range(3)]
    for i, j, k in product(range(3), repeat=3):
        answer[i][j] ^= MULTIPLICATION[first[i][k]][second[k][j]]
    return answer


def dense_controls():
    identity = [[int(i == j) for j in range(3)] for i in range(3)]
    weyl = {(a, b): [[ROOT_POWERS[b*j % 3] if i == (j+a) % 3 else 0
                      for j in range(3)] for i in range(3)] for a, b in product(range(3), repeat=2)}
    inverse = {}
    for key, matrix in weyl.items():
        a, b = key
        inverse[key] = [[MULTIPLICATION[ROOT_POWERS[a*b % 3]][v] for v in row]
                        for row in weyl[-a % 3, -b % 3]]
        assert dense_product(matrix, inverse[key]) == identity
    rng = Random(347591)
    for trial in range(384):
        length = 1+trial % 6
        physical = [[rng.randrange(3) for _ in range(4)] for _ in range(length)]
        weights = canonical_weights(rng.randrange(4096))
        matrix = identity
        a = b = phase = 0
        for x, y, u, v in physical:
            for factor in (weyl[x, 0], weyl[0, y], weyl[0, -u % 3], weyl[-v % 3, 0]):
                matrix = dense_product(matrix, factor)
            da, db = (x-v) % 3, (y-u) % 3
            phase, a, b = (phase-v*db+b*da) % 3, (a+da) % 3, (b+db) % 3
        value = 0
        for key, weight in enumerate(weights):
            product_matrix = dense_product(inverse[divmod(key, 3)], matrix)
            trace = product_matrix[0][0] ^ product_matrix[1][1] ^ product_matrix[2][2]
            value ^= MULTIPLICATION[weight][trace]
        assert value == MULTIPLICATION[weights[3*a+b]][ROOT_POWERS[phase]]
        # Direct cold specialization is a separated row/column phase.
        cold = [[x, y, x, y] for x, y, _, _ in physical]
        ca = cb = cp = 0
        for x, y, u, v in cold:
            da, db = (x-v) % 3, (y-u) % 3
            cp, ca, cb = (cp-v*db+cb*da) % 3, (ca+da) % 3, (cb+db) % 3
        assert MULTIPLICATION[weights[3*ca+cb]][ROOT_POWERS[cp]] == ROOT_POWERS[
            2*sum(x*x-y*y for x, y, _, _ in physical) % 3]
    return 384


def generate():
    records = []
    for code, jets in boundary_jets():
        r = field_rank(jets[0])
        coefficients = [scalar_restriction(a) for a in jets[1:]]
        r3, r4 = toeplitz_rank(coefficients, 3), toeplitz_rank(coefficients, 4)
        assert r*(r4-r3) >= 738
        records.append([code, r, r3, r4])
    return dict(records=records, minimum_all_scale_floor=738,
                source_sha256={name: sha256((ROOT/'tools'/name).read_bytes()).hexdigest()
                               for name in ('weyl_ordered_crossing.py', 'check_weyl_ordered_boundaries.py',
                                            'check_steiner_cubic_crossings.py')},
                scope='Every F4 boundary of the two-position linear qutrit Weyl construction with rank-one zero cut, up to the stated paired phase equivalence.',
                factoring_exponent='1/2; unchanged')


def verify(data):
    assert len(data['records']) == 4096
    saved = {row[0]: row[1:] for row in data['records']}
    assert set(saved) == set(range(4096))
    minimum = 10**9
    for code, jets in boundary_jets():
        r = field_rank(jets[0], low=True)
        coefficients = [scalar_restriction(a) for a in jets[1:]]
        r3, r4 = toeplitz_rank(coefficients, 3, low=True), toeplitz_rank(coefficients, 4, low=True)
        assert [r, r3, r4] == saved[code]
        minimum = min(minimum, r*(r4-r3))
    assert minimum == data['minimum_all_scale_floor'] == 738 > 9**3
    for name, digest in data['source_sha256'].items():
        assert sha256((ROOT/'tools'/name).read_bytes()).hexdigest() == digest
    return dict(complete_boundaries=4096, exact_leading_and_toeplitz_ranks=3*4096,
                independent_dense_weyl_controls=dense_controls(), all_scale_bare_floor=738,
                exact_comparison='738 > 9^3', strict_bare_rho_lower_bound='3/4',
                factoring_exponent='1/2; unchanged')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    parser.add_argument('--evidence', type=Path, default=ROOT/'evidence/weyl-ordered-boundaries.json')
    args = parser.parse_args()
    if args.generate:
        data = generate()
        args.evidence.write_text(json.dumps(data, separators=(',', ':'))+'\n')
    else:
        data = json.loads(args.evidence.read_text())
    print(json.dumps(verify(data)))


if __name__ == '__main__':
    main()
