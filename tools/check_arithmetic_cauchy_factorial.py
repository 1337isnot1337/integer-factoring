#!/usr/bin/env python3
"""Exact finite controls for the all-size arithmetic Cauchy factorial law.

The ordinary proof is in research/arithmetic-cauchy-factorial-law.md.
These controls are not an asymptotic factoring or individual Smith certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
import json
from math import factorial, prod
from pathlib import Path
from random import Random

from boolean_cauchy_operator import smith_orders
from walsh_cauchy_permanents import permanent, principal_permanent_polynomials


def valuation(value, prime):
    value = Fraction(value)
    if not value:
        return None
    numerator, denominator = abs(value.numerator), value.denominator
    answer = 0
    while numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator % prime == 0:
        denominator //= prime
        answer -= 1
    return answer


def falling(number, count):
    return prod(number - index for index in range(count))


@lru_cache(maxsize=None)
def positive_partitions(total, least=1):
    if total == 0:
        return ((),)
    return tuple((part,) + tail for part in range(least, total + 1)
                 for tail in positive_partitions(total - part, part))


def sampling_moment(size, exponents):
    """Independent inclusion-exclusion recurrence for the distinct sum."""
    power_sums = [sum(value ** power for value in range(size))
                  for power in range(sum(exponents) + 1)]

    @lru_cache(maxsize=None)
    def distinct_sum(powers):
        if not powers:
            return 1
        first, rest = powers[0], list(powers[1:])
        answer = power_sums[first] * distinct_sum(tuple(rest))
        for index in range(len(rest)):
            merged = rest[:]
            merged[index] += first
            answer -= distinct_sum(tuple(sorted(merged)))
        return answer

    return Fraction(distinct_sum(tuple(sorted(exponents))),
                    falling(size, len(exponents)))


def poset_data(vertices, edges):
    predecessors = [0] * vertices
    neighbors = [set() for _ in range(vertices)]
    for lower, upper in edges:
        if lower == upper:
            return None
        predecessors[upper] |= 1 << lower
        neighbors[lower].add(upper)
        neighbors[upper].add(lower)
    orders = [0] * (1 << vertices)
    orders[0] = 1
    for state, count in enumerate(orders):
        if not count:
            continue
        for vertex in range(vertices):
            if not (state >> vertex) & 1 and predecessors[vertex] & ~state == 0:
                orders[state | (1 << vertex)] += count
    if not orders[-1]:
        return None
    unseen = set(range(vertices))
    components = []
    while unseen:
        seed = unseen.pop()
        stack, count = [seed], 0
        while stack:
            vertex = stack.pop()
            count += 1
            following = neighbors[vertex] & unseen
            unseen.difference_update(following)
            stack.extend(following)
        components.append(count)
    return orders[-1], components


def equality_partition_moments(exponents):
    """Enumerate leaf equalities, retaining exact linear-extension weights."""
    roots = len(exponents)
    attached_to = tuple(root for root, count in enumerate(exponents)
                        for _ in range(count))
    records = []
    discarded_cycles = 0

    def extend(index, blocks, edges):
        nonlocal discarded_cycles
        if index == len(attached_to):
            data = poset_data(blocks, set(edges))
            if data is None:
                discarded_cycles += 1
                return
            extensions, components = data
            weight = Fraction(extensions, factorial(blocks))
            assert sum(component - 1 for component in components) <= len(set(edges))
            for prime in (2, 3, 5):
                assert valuation(weight, prime) >= -(len(attached_to) // (prime - 1))
            records.append((blocks, weight))
            return
        root = attached_to[index]
        for block in range(blocks + 1):
            if block == root:
                continue
            extend(index + 1, max(blocks, block + 1), edges + [(block, root)])

    extend(0, roots, [])
    return records, discarded_cycles


def quadratic_add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def quadratic_multiply(left, right):
    a, b = left
    c, d = right
    # omega^2 = -1 - omega.
    return a * c - b * d, a * d + b * c - b * d


def quadratic_inverse(value):
    a, b = value
    norm = a * a - a * b + b * b
    return Fraction(a - b, norm), Fraction(-b, norm)


def local_smith_orders(matrix, prime, precision):
    """Minimum-valuation elimination, independently of moment determinants."""
    modulus = prime ** precision
    work = [[value % modulus for value in row] for row in matrix]
    size = len(work)
    result = []
    for index in range(size):
        best = None
        for row in range(index, size):
            for column in range(index, size):
                if work[row][column]:
                    order = valuation(work[row][column], prime)
                    if best is None or order < best[0]:
                        best = order, row, column
        if best is None:
            return result + [None] * (size - index)
        order, row, column = best
        work[index], work[row] = work[row], work[index]
        for values in work:
            values[index], values[column] = values[column], values[index]
        divisor = prime ** order
        inverse = pow(work[index][index] // divisor, -1, modulus)
        work[index] = [value * inverse % modulus for value in work[index]]
        for row in range(index + 1, size):
            factor = work[row][index] // divisor
            if factor:
                work[row][index:] = [(a - factor * b) % modulus
                                    for a, b in zip(work[row][index:], work[index][index:])]
        work[index][index + 1:] = [0] * (size - index - 1)
        result.append(order)
    return result


def orthogonal_basis_controls(prime, step, start, size):
    moments = [Fraction(1, (start + step * index) ** 2) for index in range(2 * size + 1)]

    def pairing(left, right):
        return sum((a * b * moments[i + j] for i, a in enumerate(left)
                    for j, b in enumerate(right)), Fraction())

    polynomials, norms = [], []
    coefficient_checks = 0
    for degree in range(size + 1):
        values = [Fraction(0)] * degree + [Fraction(1)]
        for previous, norm in zip(polynomials, norms):
            factor = pairing(values, previous) / norm
            for index, value in enumerate(previous):
                values[index] -= factor * value
        norm = pairing(values, values)
        expected = (2 * degree * valuation(step, prime)
                    + 2 * valuation(factorial(degree), prime)
                    + valuation(degree + 1, prime))
        assert valuation(norm, prime) == expected
        assert valuation(values[0], prime) == 0
        assert all(not value or valuation(value, prime) >= 0 for value in values)
        assert all(pairing(values, previous) == 0 for previous in polynomials)
        coefficient_checks += len(values)
        polynomials.append(values)
        norms.append(norm)

    for degree in range(size):
        beta = norms[degree] / norms[degree - 1] if degree else Fraction(0)
        alpha = -(polynomials[degree + 1][0]
                  + (beta * polynomials[degree - 1][0] if degree else 0)) / polynomials[degree][0]
        assert valuation(alpha, prime) == 0
        assert not degree or valuation(beta, prime) > 0
        reconstructed = [Fraction(0)] + polynomials[degree]
        for index, value in enumerate(polynomials[degree]):
            reconstructed[index] -= alpha * value
        if degree:
            for index, value in enumerate(polynomials[degree - 1]):
                reconstructed[index] -= beta * value
        assert reconstructed == polynomials[degree + 1]
    return coefficient_checks, [valuation(norm, prime) for norm in norms]


def verify(max_size=14):
    if not 1 <= max_size <= 18:
        raise ValueError("Use a maximum direct matrix size from 1 through 18.")
    moment_controls = direct_moment_controls = quotient_controls = 0
    partition_records = []
    for exponents in ((1,), (2,), (1, 1), (3,), (2, 1), (1, 1, 1),
                      (4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)):
        records, discarded = equality_partition_moments(exponents)
        quotient_controls += len(records)
        for size in range(len(exponents), 9):
            reconstructed = sum((weight * falling(size - len(exponents), blocks - len(exponents))
                                 for blocks, weight in records), Fraction())
            expected = sampling_moment(size, exponents)
            assert reconstructed == expected
            if size <= 6:
                literal = sum(prod(value ** power for value, power in zip(values, exponents))
                              for values in permutations(range(size), len(exponents)))
                assert expected == Fraction(literal, falling(size, len(exponents)))
                direct_moment_controls += 1
            moment_controls += 1
        partition_records.append({'exponents': list(exponents), 'acyclic_partitions': len(records),
                                  'discarded_cyclic_partitions': discarded})

    denominator_controls = 0
    for size in (2, 3, 4, 6, 8, 12, 16, 24, 32, 48):
        for degree in range(1, 11):
            for exponents in positive_partitions(degree):
                if len(exponents) > size:
                    continue
                value = sampling_moment(size, exponents)
                for prime in (2, 3, 5, 7):
                    actual = valuation(value, prime)
                    assert actual is None or actual >= -(degree // (prime - 1))
                    denominator_controls += 1

    random = Random(2026103101)
    permanent_records = []
    for prime, steps in ((2, (2, 4, 6, 10)), (3, (3, 6)), (5, (5, 10))):
        for size in range(1, max_size + 1):
            loss = valuation(factorial(size), prime)
            modulus = prime ** (loss + 3)
            target_modulus = prime ** (loss + 1)
            for step in steps:
                for trial in range(3):
                    columns = [prime * random.randrange(-1000, 1000) + random.randrange(1, prime)
                               for _ in range(size)]
                    matrix = [[pow(column + step * row, -1, modulus) for column in columns]
                              for row in range(size)]
                    actual = permanent(matrix, modulus)
                    expected = factorial(size) * prod(pow(column, -1, modulus) for column in columns)
                    assert (actual - expected) % target_modulus == 0
                    assert valuation(actual, prime) == loss
                    permanent_records.append({'prime': prime, 'size': size, 'spacing': step,
                                              'trial': trial, 'columns': columns,
                                              'permanent_residue': actual, 'precision': loss + 3,
                                              'valuation': loss})

    # Full polynomial checks exercise arbitrary-size principal subsets.
    polynomial_records = []
    for length in (4, 8, 16):
        exponent = length.bit_length() - 1
        selected = {}
        for spacing in (2, 4):
            for start in (0, 1):
                for size in range(1, min(5, (length - 1 - start) // spacing + 1) + 1):
                    subset = sum(1 << (start + spacing * index) for index in range(size))
                    selected[subset] = (start, spacing, size)
        for center in sorted(set((1 - length, -1, 1, length - 1))):
            for subset, (start, spacing, size) in selected.items():
                content = size * (exponent + 1) + valuation(factorial(size), 2)
                polynomial = principal_permanent_polynomials(
                    length, center, content + 1, [subset])[subset]
                expected = [0] * len(polynomial)
                expected[size] = 1 << content
                assert polynomial == expected
                polynomial_records.append({'packet_length': length, 'center': center, 'start': start,
                                           'spacing': spacing, 'size': size, 'content': content,
                                           'checked_coefficients': len(polynomial)})

    positions = (0, 2, 4, 8)
    bad = sum((prod((Fraction(1, 1 + x - positions[column])
                     for x, column in zip(positions, permutation)), start=Fraction(1))
               for permutation in permutations(range(4))), Fraction())
    assert bad == Fraction(63488, 165375) and valuation(bad, 2) == 11
    assert valuation(bad, 2) != valuation(factorial(4), 2)
    omega, omega_squared = (0, 1), (-1, -1)
    shift = (2, 0)
    mixed = quadratic_add(
        quadratic_multiply(quadratic_inverse(omega),
                           quadratic_inverse(quadratic_add(omega_squared, shift))),
        quadratic_multiply(quadratic_inverse(omega_squared),
                           quadratic_inverse(quadratic_add(omega, shift))))
    assert mixed == (0, 0)

    smith_records = []
    for prime, step in ((2, 2), (2, 4), (2, 6), (3, 3), (3, 9), (5, 5)):
        for size in (1, 2, 3, 5, 8, 12, 17, 24, 33):
            expected = [2 * degree * valuation(step, prime)
                        + 2 * valuation(factorial(degree), prime)
                        + valuation(degree + 1, prime) for degree in range(size)]
            precision = expected[-1] + 4
            modulus = prime ** precision
            for start in (1, prime + 1):
                matrix = [[pow(start + step * (row + column), -2, modulus)
                           for column in range(size)] for row in range(size)]
                actual = (smith_orders(matrix, precision) if prime == 2
                          else local_smith_orders(matrix, prime, precision))
                assert actual == expected
                smith_records.append({'prime': prime, 'spacing': step, 'A': start,
                                      'size': size, 'orders': actual, 'precision': precision})
    orthogonal_records = []
    for prime, step, start in ((2, 2, 1), (2, 4, -7), (2, 6, 3),
                               (3, 3, 1), (3, 9, 4), (5, 5, 6)):
        checks, orders = orthogonal_basis_controls(prime, step, start, 6)
        orthogonal_records.append({'prime': prime, 'spacing': step, 'A': start,
                                   'maximum_degree': 6, 'coefficient_checks': checks,
                                   'norm_orders': orders})

    root = Path(__file__).resolve().parent.parent
    sources = ['research/arithmetic-cauchy-factorial-law.md',
               'tools/check_arithmetic_cauchy_factorial.py',
               'tools/walsh_cauchy_permanents.py', 'tools/cauchy_packet_polynomials.py',
               'tools/boolean_cauchy_operator.py']
    return {
        'status': 'pass',
        'claim': 'Ordinary all-size arithmetic Cauchy permanent factorial valuation and '
                 'complete squared-matrix Smith profile; no improved factoring exponent.',
        'moment_controls': {'equality_partition_reconstructions': moment_controls,
                            'literal_injective_sums': direct_moment_controls,
                            'acyclic_quotient_weights': quotient_controls,
                            'prime_denominator_checks': denominator_controls,
                            'partition_records': partition_records},
        'permanent_controls': {'count': len(permanent_records), 'maximum_size': max_size,
                               'records': permanent_records},
        'complete_polynomial_controls': {'count': len(polynomial_records),
                                         'records': polynomial_records},
        'smith_controls': {'count': len(smith_records), 'records': smith_records},
        'orthogonal_basis_controls': {'count': len(orthogonal_records), 'records': orthogonal_records},
        'negative_controls': {
            'nonconsecutive_even_positions': {'positions': list(positions), 'T': 1,
                                               'permanent': str(bad), 'valuation': 11,
                                               'factorial_valuation': 3},
            'different_residues_in_unramified_quadratic_extension': {
                'relation': 'omega^2+omega+1=0', 'spacing': 2,
                'column_shifts': ['omega', 'omega^2'], 'permanent': 0}},
        'source_sha256': {source: sha256((root / source).read_bytes()).hexdigest()
                          for source in sources},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-size', type=int, default=14)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify(args.max_size)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': result['status'],
                      'moment_reconstructions': result['moment_controls']['equality_partition_reconstructions'],
                      'prime_denominator_checks': result['moment_controls']['prime_denominator_checks'],
                      'permanents': result['permanent_controls']['count'],
                      'complete_packet_polynomials': result['complete_polynomial_controls']['count'],
                      'smith_profiles': result['smith_controls']['count'],
                      'orthogonal_bases': result['orthogonal_basis_controls']['count'],
                      'negative_controls': list(result['negative_controls'])}, indent=2))


if __name__ == '__main__':
    main()
