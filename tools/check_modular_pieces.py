"""Exact finite checks for the prime-field quotient crossing construction.

This checks identities and finite-field ranks, not an asymptotic saving.
The default suite is small; --full adds the recorded weight-16 F_2 case.
Python's standard library suffices. The shared matrix routine optionally
uses NumPy, with its existing overflow bound.
"""

from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from math import comb, factorial, isqrt, prod
from pathlib import Path
import argparse
import json

from check_half_piece_ranks import (
    afam as source_afam,
    matmul_mod,
    modrank,
    opcoef as source_opcoef,
    partitions,
    remove,
    zeros,
)


@lru_cache(None)
def regular_partitions(weight, prime, maxpart=None):
    """No p-divisible parts, and every multiplicity is less than p."""
    if weight < 0:
        return ()
    if weight == 0:
        return ((),)
    if maxpart is None:
        maxpart = weight
    answer = []
    for part in range(min(weight, maxpart), 0, -1):
        if part % prime == 0:
            continue
        for rest in regular_partitions(weight - part, prime, part):
            if rest.count(part) < prime - 1:
                answer.append((part,) + rest)
    return tuple(answer)


def z_value(partition):
    return prod(k ** count * factorial(count)
                for k, count in Counter(partition).items())


@lru_cache(None)
def quotient_afam(partition, prime):
    """Source a^mu reduced modulo (t_{kp}, t_k^p), coefficientwise."""
    terms = {(): 1}
    for part in partition:
        next_terms = defaultdict(int)
        for monomial, value in terms.items():
            for j in range(1, part + 1):
                if j % prime == 0 or monomial.count(j) == prime - 1:
                    continue
                coefficient = (-1) ** (part - 1) * comb(part, j) % prime
                if coefficient:
                    key = tuple(sorted(monomial + (j,), reverse=True))
                    next_terms[key] += value * coefficient
        terms = {key: value % prime for key, value in next_terms.items()
                 if value % prime}
    return tuple(terms.items())


@lru_cache(None)
def quotient_opcoef(operator, output, state, prime):
    """Source normal-order coefficient recurrence, calculated modulo p."""
    if not operator:
        return int(output == state)
    part, rest = operator[0], operator[1:]
    value = 0
    for k in set(output):
        if k >= part:
            coefficient = comb(k - 1, part - 1) % prime
            if coefficient:
                value += coefficient * quotient_opcoef(
                    rest, remove(output, k), state, prime)
    for k, count in Counter(state).items():
        coefficient = ((-1) ** (part - 1) *
                       comb(k + part - 1, part - 1) * k * count) % prime
        if coefficient:
            value += coefficient * quotient_opcoef(
                rest, output, remove(state, k), prime)
    return value % prime


def transpose_middle(matrix, dimension):
    """Swap the y and u roles; the crossing itself is fixed."""
    d = dimension
    return [[matrix[x*d + y][u*d + v]
             for y in range(d) for v in range(d)]
            for x in range(d) for u in range(d)]


def binary_vector(matrix):
    answer = 0
    for row in matrix:
        for value in row:
            answer = (answer << 1) | (value & 1)
    return answer


def binary_span_contains(vectors, target):
    basis = {}
    for value in vectors:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    while target:
        pivot = target.bit_length() - 1
        if pivot not in basis:
            return len(basis), False
        target ^= basis[pivot]
    return len(basis), True


def finite_recombination(matrices, dimension):
    """Verify the discovered 17-term F_2 crossing on the weight-12 domain."""
    certificate = [
        (-24, (0, 2, 1, 3)),
        (-22, (0, 1, 2, 3)), (-22, (0, 2, 1, 3)),
        (-22, (1, 3, 0, 2)), (-22, (2, 3, 0, 1)),
        (-16, (0, 2, 1, 3)),
        (-14, (0, 1, 2, 3)), (-14, (1, 0, 3, 2)),
        (-14, (1, 3, 0, 2)),
        (-12, (0, 1, 2, 3)), (-12, (1, 0, 3, 2)),
        (-6, (0, 1, 2, 3)), (-6, (1, 0, 3, 2)),
        (8, (0, 1, 2, 3)), (16, (0, 1, 2, 3)),
        (24, (0, 1, 2, 3)), (24, (0, 2, 1, 3)),
    ]
    d = dimension
    total = zeros(d*d, d*d)
    first_piece = None
    record = []
    for shift, permutation in certificate:
        assert {frozenset((permutation[0], permutation[3])),
                frozenset((permutation[1], permutation[2]))} == {
                    frozenset((0, 3)), frozenset((1, 2))}
        original = matrices[shift]
        transformed = zeros(d*d, d*d)
        for x in range(d):
            for y in range(d):
                for u in range(d):
                    for v in range(d):
                        coordinates = (x, y, u, v)
                        a, b, c, e = (coordinates[k] for k in permutation)
                        transformed[x*d+y][u*d+v] = original[a*d+b][c*d+e]
        row_rank = modrank(transformed, 2)
        column_rank = modrank(transpose_middle(transformed, d), 2)
        assert 0 < row_rank * column_rank <= 27
        if first_piece is None:
            first_piece = transformed
        total = [[old ^ new for old, new in zip(oldrow, newrow)]
                 for oldrow, newrow in zip(total, transformed)]
        record.append({'s': shift, 'permutation': list(permutation),
                       'row_rank': row_rank, 'column_rank': column_rank})
    crossing = [[int(x == v and y == u)
                 for u in range(d) for v in range(d)]
                for x in range(d) for y in range(d)]
    assert total == crossing
    omitted = [[old ^ removed for old, removed in zip(oldrow, removedrow)]
               for oldrow, removedrow in zip(total, first_piece)]
    assert omitted != crossing
    return {'pieces': record, 'max_rank_product': 27,
            'crossing_entries_checked': d ** 4,
            'omitted_piece_negative_control_passed': True}


def one_case(weight, prime):
    assert prime >= 2
    assert all(prime % d for d in range(2, isqrt(prime) + 1))
    labels = tuple(sorted(regular_partitions(weight, prime)))
    dimension = len(labels)
    assert dimension > 0
    inverses = [pow(z_value(label), -1, prime) for label in labels]
    duality_count = 0
    for x, first in enumerate(labels):
        for v, second in enumerate(labels):
            result = inverses[v] * sum(
                coefficient * quotient_opcoef(mu, (), second, prime)
                for mu, coefficient in quotient_afam(first, prime))
            assert result % prime == int(x == v)
            duality_count += 1

    d = dimension
    total = zeros(d*d, d*d)
    profiles = []
    binary_pairs = []
    recombination_matrices = {}
    order_entries = 0
    for shift in range(-2 * weight, 2 * weight + 1):
        channels = [(i, partition) for i in range(2 * weight + 1)
                    for partition in regular_partitions(shift + i, prime)]
        left = zeros(d*d, len(channels))
        right = zeros(len(channels), d*d)
        for column, (i, intermediate) in enumerate(channels):
            for x, first in enumerate(labels):
                terms = [(mu, coefficient * (-1) ** i *
                          comb(weight + sum(mu), i) % prime)
                         for mu, coefficient in quotient_afam(first, prime)
                         if weight + sum(mu) >= i]
                terms = [(mu, coefficient) for mu, coefficient in terms
                         if coefficient]
                for y, state in enumerate(labels):
                    left[x*d + y][column] = sum(
                        coefficient * quotient_opcoef(
                            mu, intermediate, state, prime)
                        for mu, coefficient in terms) % prime
            for u, output in enumerate(labels):
                for v, second in enumerate(labels):
                    right[column][u*d + v] = inverses[v] * quotient_opcoef(
                        second, output, intermediate, prime) % prime
        matrix = matmul_mod(left, right, prime)
        if weight == 12 and prime == 2:
            recombination_matrices[shift] = matrix
        total = [[(x+y) % prime for x, y in zip(old, new)]
                 for old, new in zip(total, matrix)]
        column_matrix = transpose_middle(matrix, d)

        # Independently compute the other operator order, at every shift.
        other_channels = [
            (i, partition) for i in range(2 * weight + 1)
            for partition in regular_partitions(2 * weight - shift - i, prime)
        ]
        other_left = zeros(d*d, len(other_channels))
        other_right = zeros(len(other_channels), d*d)
        for column, (i, intermediate) in enumerate(other_channels):
            for x, first in enumerate(labels):
                terms = [(mu, coefficient * (-1) ** i *
                          comb(weight + sum(mu), i) % prime)
                         for mu, coefficient in quotient_afam(first, prime)
                         if weight + sum(mu) >= i]
                terms = [(mu, coefficient) for mu, coefficient in terms
                         if coefficient]
                for u, output in enumerate(labels):
                    other_left[x*d + u][column] = sum(
                        coefficient * quotient_opcoef(
                            mu, output, intermediate, prime)
                        for mu, coefficient in terms) % prime
            for y, state in enumerate(labels):
                for v, second in enumerate(labels):
                    other_right[column][y*d + v] = inverses[v] * quotient_opcoef(
                        second, intermediate, state, prime) % prime
        assert matmul_mod(other_left, other_right, prime) == column_matrix
        order_entries += d ** 4
        row_rank = modrank(matrix, prime)
        column_rank = modrank(column_matrix, prime)
        profiles.append({
            's': shift, 'row_rank': row_rank, 'column_rank': column_rank,
            'row_channels': len(channels),
            'column_channels': len(other_channels),
        })
        if prime == 2:
            binary_pairs.append((binary_vector(matrix),
                                 binary_vector(column_matrix)))

    crossing = [[int(x == v and y == u)
                 for u in range(d) for v in range(d)]
                for x in range(d) for y in range(d)]
    assert total == crossing
    answer = {
        'weight': weight, 'prime': prime, 'labels': [list(p) for p in labels],
        'dimension': dimension, 'normalization': 1, 'ranks': profiles,
        'integer_h_factorial_mod_p': factorial(weight) % prime,
        'duality_entries': duality_count, 'crossing_entries': d ** 4,
        'opposite_order_entries': order_entries,
        'max_rank_product': max(
            row['row_rank'] * row['column_rank'] for row in profiles),
    }
    if prime == 2:
        products = sorted({row['row_rank'] * row['column_rank']
                           for row in profiles})
        target = binary_vector(crossing)
        tested = []
        for cutoff in products:
            vectors = [vector for pair, row in zip(binary_pairs, profiles)
                       if row['row_rank'] * row['column_rank'] <= cutoff
                       for vector in pair]
            rank, contains = binary_span_contains(vectors, target)
            tested.append({'cutoff': cutoff, 'span_rank': rank,
                           'contains_crossing': contains})
            if contains:
                break
        assert tested[-1]['contains_crossing']
        answer['piece_and_exchange_span_screen'] = tested
        answer['minimum_spanning_rank_product_cutoff'] = tested[-1]['cutoff']
    if weight == 12 and prime == 2:
        answer['finite_recombination'] = finite_recombination(
            recombination_matrices, d)
    quotient_opcoef.cache_clear()
    return answer


def source_checks():
    exact_duality = 0
    for weight in range(1, 7):
        for first in partitions(weight):
            for second in partitions(weight):
                result = sum(coefficient * source_opcoef(mu, (), second)
                             for mu, coefficient in source_afam(first))
                assert result == z_value(second) * int(first == second)
                exact_duality += 1

    quotient_coefficients = 0
    null_labels = 0
    for prime in (2, 3, 5):
        states = [partition for weight in range(5)
                  for partition in regular_partitions(weight, prime)]
        outputs = [partition for weight in range(7)
                   for partition in regular_partitions(weight, prime)]
        operators = [partition for weight in range(5)
                     for partition in partitions(weight)]
        for operator in operators:
            for output in outputs:
                for state in states:
                    assert quotient_opcoef(operator, output, state, prime) == (
                        source_opcoef(operator, output, state) % prime)
                    quotient_coefficients += 1
        for operator in ((prime,), (1,) * prime, (prime, 1), (2,) * prime):
            for output in outputs:
                for state in states:
                    assert source_opcoef(operator, output, state) % prime == 0
                    null_labels += 1
        for weight in range(1, 7):
            for partition in regular_partitions(weight, prime):
                expected = {
                    mu: coefficient % prime
                    for mu, coefficient in source_afam(partition)
                    if coefficient % prime and all(k % prime for k in mu)
                    and all(count < prime for count in Counter(mu).values())
                }
                assert dict(quotient_afam(partition, prime)) == expected

    # The old normalized b has zero coefficient for every regular label
    # when h >= p; it cannot give the new scale-one identity.
    bad_normalization = []
    for prime, weight in ((2, 8), (3, 3), (5, 5)):
        labels = regular_partitions(weight, prime)
        assert labels
        assert all((factorial(weight) // z_value(mu)) % prime == 0
                   for mu in labels)
        bad_normalization.append({'prime': prime, 'weight': weight,
                                  'old_scale': 0, 'required_scale': 1})
    # The state ideal is not stable over Q: d(t_1^p)=p t_1^(p-1).
    characteristic_zero_failure = source_opcoef((1,), (1,), (1, 1))
    assert characteristic_zero_failure == 2
    assert characteristic_zero_failure % 2 == 0
    quotient_opcoef.cache_clear()
    source_opcoef.cache_clear()
    return {
        'integer_unscaled_duality_entries': exact_duality,
        'source_mod_quotient_coefficient_comparisons': quotient_coefficients,
        'null_operator_label_entries': null_labels,
        'zero_factorial_negative_controls': bad_normalization,
        'rational_ideal_failure_derivative_coefficient': characteristic_zero_failure,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    cases = [(8, 2), (12, 2), (3, 3), (4, 3), (5, 5)]
    if args.full:
        cases += [(16, 2), (6, 3)]
    checks = source_checks()
    records = [one_case(weight, prime) for weight, prime in cases]
    here = Path(__file__).resolve()
    record = {
        'claim': 'Exact modular crossing construction and finite ranks; no stronger rate',
        'source_revision': '16a5ce41304fcdaea4e0be7c132e1c361dc26da5',
        'checker_sha256': sha256(here.read_bytes()).hexdigest(),
        'shared_operator_sha256': sha256(
            here.with_name('check_half_piece_ranks.py').read_bytes()).hexdigest(),
        'full_suite': args.full, 'source_checks': checks, 'cases': records,
        'totals': {
            key: sum(record[key] for record in records)
            for key in ('duality_entries', 'crossing_entries', 'opposite_order_entries')
        },
        'all_checks_passed': True,
    }
    output = json.dumps(record, indent=2, sort_keys=True) + '\n'
    if args.output:
        args.output.write_text(output)
        print(json.dumps({'output': str(args.output), 'cases': len(records),
                          'totals': record['totals'], 'all_checks_passed': True}))
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
