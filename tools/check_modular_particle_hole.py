#!/usr/bin/env python3
"""Exact checks of reflected particle-hole compression of modular profiles.

The all-size proof is in research/modular-particle-hole.md. These finite
checks do not establish an improved asymptotic factoring exponent.
"""

import argparse
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random

from check_modular_crossing import regular_partitions
from modular_pfaffian import (
    final_piece_profile, folded_matrix, polynomial_product,
    polynomial_exact_quotient, polynomial_pfaffian,
)


def power(polynomial, exponent):
    answer = 1
    while exponent:
        if exponent & 1:
            answer = polynomial_product(answer, polynomial)
        polynomial = polynomial_product(polynomial, polynomial)
        exponent >>= 1
    return answer


def complement(label, k):
    return tuple(sorted(set(range(1, 2*k, 2)) - {2*k-t for t in label}))


def multiply(left, right):
    return [[sum_product(left[i], [row[j] for row in right])
             for j in range(len(right[0]))] for i in range(len(left))]


def sum_product(left, right):
    value = 0
    for a, b in zip(left, right):
        value ^= polynomial_product(a, b)
    return value


def transpose(matrix):
    return list(map(list, zip(*matrix)))


def multiply_chain(*matrices):
    answer = matrices[0]
    for matrix in matrices[1:]:
        answer = multiply(answer, matrix)
    return answer


def digit_identities():
    # The indeterminate b is encoded by polynomial 2; equality checks every
    # coefficient of the two-by-two identities used in the all-size proof.
    p, r = [[1, 1], [0, 1]], [[0, 1], [1, 0]]
    w, v = [[1, 0], [0, 2]], [[1, 0], [0, 3]]
    j, h, o = [[1, 1], [1, 0]], [[1, 1], [1, 2]], [[1, 2], [2, 2]]
    pt, jt = transpose(p), transpose(j)
    pairs = [((pt, v, r, pt, w), (w, r, h)),
             ((o, r, w, p), (w, r, w, jt, v)),
             ((o, r, h), (pt, v, r, v, j, w)),
             ((o, r, pt, v), (h, r, w, jt, v)),
             ((o, r, w), (w, p, r, v, j, w)),
             ((v, p, r, w), (pt, w, r, h))]
    for left, right in pairs:
        assert multiply_chain(*left) == multiply_chain(*right)
    return len(pairs)*4


def core_case(k):
    core = tuple(range(1, 2*k, 2))
    matrix = folded_matrix(core, core, core, core)
    size = 4*k
    reflection = [(i//k)*k+k-1-i%k for i in range(size)]
    scalar = 1 | (1 << (2*k)) | (1 << (4*k))
    entries = 0
    for i in range(size):
        for j in range(size):
            result = 0
            for t in range(size):
                result ^= polynomial_product(matrix[i][t],
                                             matrix[reflection[t]][reflection[j]])
            assert result == (scalar if i == j else 0)
            entries += 1
    pfaffian_checked = k <= 8
    if pfaffian_checked:
        assert polynomial_pfaffian(matrix) == power(scalar, k)
    return {'core_modes': k, 'inverse_identity_entries': entries,
            'full_pfaffian_checked': pfaffian_checked}


def profile_case(k, weight, cardinality, full):
    labels = tuple(label for label in regular_partitions(2, weight)
                   if len(label) == cardinality and max(label, default=0) < 2*k)
    assert len(labels) >= 2
    labels = labels[:4]
    transformed = tuple(complement(label, k) for label in labels)
    new_weight = k*k-2*k*cardinality+weight
    for label, other in zip(labels, transformed):
        assert complement(other, k) == label
        assert sum(other) == new_weight and len(other) == k-cardinality
    scalar = 1 | (1 << (2*k)) | (1 << (4*k))
    common = power(scalar, abs(k-2*cardinality))
    tuples = list(product(range(len(labels)), repeat=4))
    if k >= 16 and not full:
        random.Random(1729).shuffle(tuples)
        tuples = tuples[:32]
    omitted_factor_fails = False
    missing_reflection_fails = False
    count = 0
    for ix in tuples:
        x, y, u, v = (labels[i] for i in ix)
        cx, cy, cu, cv = (transformed[i] for i in ix)
        old = final_piece_profile(x, v, u, y)
        new = final_piece_profile(cx, cv, cu, cy)
        if k >= 2*cardinality:
            assert new == polynomial_product(common, old)
        else:
            assert old == polynomial_product(common, new)
            assert polynomial_exact_quotient(old, common) == new
        assert old.bit_count() % 2 == new.bit_count() % 2 == int(ix[0] == ix[3]
                                                              and ix[1] == ix[2])
        omitted_factor_fails |= old != new
        if weight != k*cardinality:
            wrong = [tuple(sorted(set(range(1, 2*k, 2))-set(label)))
                     for label in (x, y, u, v)]
            wx, wy, wu, wv = wrong
            wrong_profile = final_piece_profile(wx, wv, wu, wy)
            missing_reflection_fails |= wrong_profile != new
        count += 1
    if k != 2*cardinality:
        assert omitted_factor_fails
    if weight != k*cardinality:
        assert missing_reflection_fails
    polynomial_product.cache_clear()
    return {'core_modes': k, 'weight': weight, 'cardinality': cardinality,
            'selected_labels': labels, 'complement_weight': new_weight,
            'complete_profiles_checked': count,
            'omitted_factor_negative_control': omitted_factor_fails,
            'missing_reflection_negative_control': missing_reflection_fails}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    digit_entries = digit_identities()
    cores = [core_case(k) for k in ((1, 2, 4, 8, 16, 32) if args.full
                                  else (1, 2, 4, 8, 16))]
    settings = [(4, 8, 2), (8, 12, 2), (8, 16, 2), (8, 20, 2),
                (8, 48, 6), (8, 32, 4), (16, 48, 4)]
    if args.full:
        settings += [(16, 96, 6), (16, 128, 6)]
    profiles = [profile_case(*setting, args.full) for setting in settings]
    here = Path(__file__)
    result = {
        'status': 'exact checks of an all-size compression identity; no improved exponent',
        'full': args.full, 'digit_identity_entries': digit_entries,
        'core_cases': cores, 'profile_cases': profiles,
        'sha256': {name: sha256((here.parent/name).read_bytes()).hexdigest()
                   for name in (here.name, 'modular_pfaffian.py')},
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'full': args.full,
                      'digit_identity_entries': digit_entries,
                      'inverse_identity_entries': sum(c['inverse_identity_entries'] for c in cores),
                      'complete_profiles_checked': sum(p['complete_profiles_checked'] for p in profiles)},
                     indent=2))


if __name__ == '__main__':
    main()
