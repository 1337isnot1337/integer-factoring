#!/usr/bin/env python3
"""Exact controls for fixed rational-frame coefficient transfer."""

import argparse
import hashlib
from itertools import product
import json
from math import comb

from rational_frame_cones import FixedRationalCone, multiply


ZERO = (0, 0, 0)
X = (1, 0, 0)
Y = (0, 1, 0)
Z = (0, 0, 1)


def rank(matrix, prime):
    work = [row[:] for row in matrix]
    pivot = 0
    for column in range(len(work[0]) if work else 0):
        row = next((i for i in range(pivot, len(work)) if work[i][column] % prime), None)
        if row is None:
            continue
        work[pivot], work[row] = work[row], work[pivot]
        inverse = pow(work[pivot][column], -1, prime)
        work[pivot] = [value * inverse % prime for value in work[pivot]]
        for i in range(len(work)):
            if i != pivot:
                factor = work[i][column]
                work[i] = [(a - factor * b) % prime for a, b in zip(work[i], work[pivot])]
        pivot += 1
        if pivot == len(work):
            break
    return pivot


def build_family(prime):
    one = {ZERO: 1}
    q = {Y: 1, Z: -1, X: -1}
    functions = [(one, one), ({}, one), ({X: 1}, q), (q, one), ({ZERO: -1}, q),
                 (one, {X: 1, (0, 1, 1): 1})]
    cone = FixedRationalCone(prime, 3, functions)
    # T = L diag(1,x) R, with q=y-z-x.
    left = [[0, 1], [2, 0]]
    right = [[3, 0], [1, 4]]
    tensor = [[q, one], [{X: 1}, {}]]
    return cone, left, right, tensor


def direct_entry(tensor, row, column, prime):
    polynomial = {ZERO: 1}
    for a, b in zip(row, column):
        polynomial = multiply(polynomial, tensor[a][b], prime)
    return polynomial


def coefficient_factorization(cone, left, right, length, target):
    words = list(product(range(2), repeat=length))
    size = len(words)
    result = [[0] * size for _ in range(size)]
    prime = cone.prime
    lower_left = min(cone.lower_weight((index,)) for row in left for index in row
                     if cone.lower_weight((index,)) is not None)
    lower_right = min(cone.lower_weight((index,)) for row in right for index in row
                      if cone.lower_weight((index,)) is not None)
    largest_split = 0
    rank_one_terms = 0
    channel_count = 0
    split_vectors = []
    for channel in words:
        order = sum(channel)
        if order > target[0]:
            continue
        channel_count += 1
        residual = (target[0] - order, target[1], target[2])
        total_weight = cone.weight(residual)
        upper_left = total_weight - length * lower_right
        upper_right = total_weight - length * lower_left
        rows = [cone.product_series([left[a][b] for a, b in zip(row, channel)], upper_left)
                for row in words]
        columns = [cone.product_series([right[a][b] for a, b in zip(channel, column)], upper_right)
                   for column in words]
        exponents = set().union(*(set(row) for row in rows))
        terms = 0
        for exponent in exponents:
            complement = tuple(a - b for a, b in zip(residual, exponent))
            a = [row.get(exponent, 0) for row in rows]
            b = [column.get(complement, 0) for column in columns]
            if not any(a) or not any(b):
                continue
            assert exponent[0] >= 0 and complement[0] >= 0
            terms += 1
            split_vectors.append(exponent)
            for i in range(size):
                for j in range(size):
                    result[i][j] = (result[i][j] + a[i] * b[j]) % prime
        largest_split = max(largest_split, terms)
        rank_one_terms += terms
    expected_channels = sum(comb(length, j) for j in range(min(length, target[0]) + 1)) if target[0] >= 0 else 0
    assert channel_count == expected_channels
    # The proof's uniform coordinate bound with |target_i|<=length.
    shifts = cone.shifts
    shift_size = max((max(map(abs, exponent)) for exponent in shifts), default=0)
    generator_size = max((max(map(abs, exponent)) for exponent in cone.generators), default=0)
    minimum_generator_weight = min((cone.weight(exponent) for exponent in cone.generators), default=1)
    allowance = max(0, sum(cone.weights) - lower_left - lower_right)
    coordinate_bound = length * shift_size + generator_size * (length * allowance // minimum_generator_weight)
    assert all(max(map(abs, exponent)) <= coordinate_bound for exponent in split_vectors)
    assert rank_one_terms <= channel_count * (2 * coordinate_bound + 1) ** 3
    return result, dict(channel_words=channel_count, rank_one_terms=rank_one_terms,
                        largest_splits_for_one_word=largest_split,
                        proved_split_coordinate_bound=coordinate_bound)


def build_receipt():
    controls = []
    scalar_checks = 0
    coefficients_compared = 0
    cone_data = []
    for prime in (2, 5, 7):
        cone, left, right, tensor = build_family(prime)
        assert cone.weights == (4, 1, 2)
        assert all(exponent[0] >= 0 and cone.weight(exponent) > 0 for exponent in cone.generators)
        cone_data.append(dict(prime=prime, weights=list(cone.weights),
                              generators=[list(exponent) for exponent in cone.generators],
                              shifts=[list(exponent) for exponent in cone.shifts]))
        for length in range(1, 5):
            for j in range(2 * length + 1):
                exponent = (j, -length - j, -length - j)
                expected = (-1) ** j * comb(length + j - 1, j) % prime
                assert cone.coefficient([5] * length, exponent) == expected
                scalar_checks += 1
        # The non-Laurent tail is present: -1/(y-z-x) has this coefficient.
        assert cone.coefficient([4], (0, -9, 8)) == -1 % prime
        scalar_checks += 1
        for length in range(1, 5):
            words = list(product(range(2), repeat=length))
            direct = [[direct_entry(tensor, row, column, prime) for column in words] for row in words]
            if length <= 3:
                targets = [exponent for exponent in product(range(length + 1), repeat=3)
                           if sum(exponent) <= length]
            else:
                targets = [(0, 0, 0), (1, 1, 1), (2, 1, 1), (4, 0, 0)]
            targets += [(-1, 0, 0), (1, -1, 0), (length, -length, 0)]
            for target in targets:
                literal = [[entry.get(target, 0) % prime for entry in row] for row in direct]
                reconstructed, counts = coefficient_factorization(cone, left, right, length, target)
                assert reconstructed == literal
                coefficients_compared += len(words) ** 2
                actual_rank = rank(literal, prime)
                assert actual_rank <= counts['rank_one_terms']
                controls.append(dict(prime=prime, length=length, target=list(target),
                                     coefficient_rank=actual_rank, **counts))
    invalid = 0
    try:
        FixedRationalCone(5, 3, [({ZERO: 1}, {X: 1})])
    except ValueError:
        invalid += 1
    assert invalid == 1
    receipt = dict(schema='fixed-rational-frame-transfer-v1',
                   scope='Fixed rational x-adically regular frames, not arbitrary formal series',
                   cone_certificates=cone_data, exact_scalar_coefficient_checks=scalar_checks,
                   coefficient_matrix_controls=controls,
                   literal_matrix_entries_compared=coefficients_compared,
                   rejected_x_pole_controls=invalid, factoring_exponent_improvement=False)
    encoded = json.dumps(receipt, sort_keys=True, separators=(',', ':')).encode()
    receipt['controls_sha256'] = hashlib.sha256(encoded).hexdigest()
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    receipt = build_receipt()
    with open(args.output, 'w', encoding='utf-8') as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps(dict(status='passed', scalar_checks=receipt['exact_scalar_coefficient_checks'],
                          matrix_controls=len(receipt['coefficient_matrix_controls']),
                          entries_compared=receipt['literal_matrix_entries_compared'],
                          controls_sha256=receipt['controls_sha256']), sort_keys=True))


if __name__ == '__main__':
    main()
