#!/usr/bin/env python3
"""Polynomial-space quotient factors over each fixed prime field.

The production routines use no memoization.  Reference recursion and
partition enumeration are imported only by the finite verification code.
"""

from collections import Counter
from itertools import product
from math import comb, factorial
import argparse
import json


def polynomial_product(left, right, degree, prime):
    result = [0] * (degree + 1)
    for i, c in enumerate(left):
        if not c or i > degree:
            continue
        for j, d in enumerate(right[:degree - i + 1]):
            result[i + j] = (result[i + j] + c * d) % prime
    return result


def streamed_permanent(matrix, degree, prime):
    size = len(matrix)
    if not size:
        return [1] + [0] * degree
    assert all(len(row) == size for row in matrix)
    result = [0] * (degree + 1)
    for subset in range(1 << size):
        term = [1]
        for row in matrix:
            row_sum = [0] * (degree + 1)
            for column in range(size):
                if subset >> column & 1:
                    for i, c in enumerate(row[column][:degree + 1]):
                        row_sum[i] = (row_sum[i] + c) % prime
            term = polynomial_product(term, row_sum, degree, prime)
            if not any(term):
                break
        sign = (-1) ** (size - subset.bit_count())
        for i, c in enumerate(term):
            result[i] = (result[i] + sign * c) % prime
    return result


def regular(parts, prime):
    return all(k > 0 and k % prime for k in parts) and all(
        count < prime for count in Counter(parts).values())


def unaffected_data(prime, operator, output, input_):
    """Yield signed repeated columns and the exact unit multiplier."""
    remaining = len(output) + len(input_) - len(operator)
    if remaining < 0 or remaining % 2:
        return
    remaining //= 2
    incoming, outgoing = Counter(input_), Counter(output)
    common = sorted(incoming.keys() & outgoing.keys())
    modes = sorted(incoming.keys() | outgoing.keys())
    ranges = [range(min(incoming[k], outgoing[k]) + 1) for k in common]
    for choices in product(*ranges):
        if sum(choices) != remaining:
            continue
        unchanged = dict(zip(common, choices))
        columns = []
        weight = 1
        for k in modes:
            eta = unchanged.get(k, 0)
            annihilated, created = incoming[k] - eta, outgoing[k] - eta
            denominator = factorial(eta) * factorial(annihilated) * factorial(created)
            assert denominator % prime
            weight = (weight * factorial(incoming[k]) * pow(k, annihilated, prime)
                      * pow(denominator, -1, prime)) % prime
            columns.extend([(True, k)] * created)
            columns.extend([(False, k)] * annihilated)
        assert len(columns) == len(operator)
        yield columns, weight


def scalar_entry(part, creation, mode):
    if creation:
        return comb(mode - 1, part - 1) if mode >= part else 0
    return (-1) ** (part - 1) * comb(mode + part - 1, part - 1)


def streamed_operator(prime, label, output, input_, normalized=False):
    if not all(regular(parts, prime) for parts in (label, output, input_)):
        raise ValueError("the production interface requires retained state partitions")
    result = 0
    for columns, weight in unaffected_data(prime, label, output, input_):
        matrix = [[[scalar_entry(n, creation, k) % prime] for creation, k in columns]
                  for n in label]
        result = (result + weight * streamed_permanent(matrix, 0, prime)[0]) % prime
    if normalized:
        denominator = 1
        for mode, multiplicity in Counter(label).items():
            denominator = denominator * pow(mode, multiplicity, prime) * factorial(multiplicity) % prime
        result = result * pow(denominator, -1, prime) % prime
    return result


def folded_entry(prime, part, creation, mode):
    result = [0] * (part + 1)
    for lowered in range(1, part + 1):
        coefficient = ((-1) ** (part - 1) * comb(part, lowered)
                       * scalar_entry(lowered, creation, mode)) % prime
        if coefficient:
            for index in range(lowered + 1):
                result[index] = (result[index] + coefficient * comb(lowered, index)) % prime
    return result


def streamed_top_profile(prime, label, output, input_):
    if not all(regular(parts, prime) for parts in (label, output, input_)):
        raise ValueError("the production interface requires retained state partitions")
    weight = sum(label)
    result = [0] * (weight + 1)
    for columns, multiplier in unaffected_data(prime, label, output, input_):
        matrix = [[folded_entry(prime, n, creation, k) for creation, k in columns]
                  for n in label]
        permanent = streamed_permanent(matrix, weight, prime)
        for index, value in enumerate(permanent):
            result[index] = (result[index] + multiplier * value) % prime
    return polynomial_product(result, [comb(weight, i) % prime for i in range(weight + 1)],
                              2 * weight, prime)


def verify():
    from check_modular_crossing import modular_angle, modular_operator, normalized_piece, regular_partitions
    from check_partition_channels import partition_denominator
    top_entries = bottom_entries = 0
    largest_particles = 0
    for prime, maximum in ((3, 7), (5, 5)):
        states = tuple(mu for energy in range(maximum + 1) for mu in regular_partitions(prime, energy))
        for weight in range(1, maximum + 1):
            for label in regular_partitions(prime, weight):
                largest_particles = max(largest_particles, len(label))
                normalization = pow(partition_denominator(label), -1, prime)
                for output in states:
                    for input_ in states:
                        profile = streamed_top_profile(prime, label, output, input_)
                        assert len(profile) == 2 * weight + 1
                        for index, value in enumerate(profile):
                            assert value == modular_angle(prime, label, index, output, input_)
                            top_entries += 1
                        assert streamed_operator(prime, label, output, input_) == modular_operator(
                            prime, label, output, input_)
                        assert streamed_operator(prime, label, output, input_, normalized=True) == (
                            normalization * modular_operator(prime, label, output, input_)) % prime
                        bottom_entries += 2
    complete_entries = 0
    for prime, weight in ((3, 2), (3, 3), (5, 2)):
        labels = regular_partitions(prime, weight)
        dimension = len(labels)
        total = [[0] * dimension ** 2 for _ in range(dimension ** 2)]
        for shift in range(-2 * weight, 2 * weight + 1):
            matrix = [[0] * dimension ** 2 for _ in range(dimension ** 2)]
            for index in range(2 * weight + 1):
                for intermediate in regular_partitions(prime, shift + index):
                    left = [(-1) ** index * streamed_top_profile(prime, x, intermediate, y)[index]
                            for x in labels for y in labels]
                    right = [streamed_operator(prime, v, u, intermediate, normalized=True)
                             for v in labels for u in labels]
                    for r, c in enumerate(left):
                        for j, d in enumerate(right):
                            matrix[r][j] = (matrix[r][j] + c * d) % prime
            assert matrix == normalized_piece(prime, weight, shift, labels)
            for r in range(dimension ** 2):
                for c in range(dimension ** 2):
                    total[r][c] = (total[r][c] + matrix[r][c]) % prime
            complete_entries += dimension ** 4
        assert total == [[int(r == c) for c in range(dimension ** 2)]
                         for r in range(dimension ** 2)]

    # These entries separately require duplicate-column division and falling factors.
    assert streamed_operator(3, (1, 1), (1, 1), ()) == 1
    assert streamed_permanent([[[1], [1]], [[1], [1]]], 0, 3)[0] == 2
    assert streamed_operator(3, (1,), (1,), (1, 1)) == 2
    assert streamed_operator(3, (2,), (), (1,)) == 1  # Annihilator sign.
    assert streamed_top_profile(3, (2,), (1,), ()) == [1, 0, 0, 1, 0]
    assert [(-x) % 3 for x in streamed_top_profile(3, (2,), (1,), ())] != [1, 0, 0, 1, 0]
    return {
        "status": "pass",
        "scope": "finite source comparisons for the proved polynomial-space fixed-prime factor interface",
        "primes": [3, 5],
        "top_angle_entries": top_entries,
        "bottom_factor_entries_including_normalization": bottom_entries,
        "complete_piece_entries": complete_entries,
        "largest_tested_operator_particle_count": largest_particles,
        "negative_controls": ["repeated-column division", "derivative multiplicity",
                              "annihilator sign", "dual sign"],
        "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    encoded = json.dumps(verify(), indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end="")
