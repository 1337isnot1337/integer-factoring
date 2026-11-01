#!/usr/bin/env python3
"""Exact Cauchy-minor witnesses for the proposed exterior-power deformation."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import permutations
import json
from math import comb, prod
from pathlib import Path


def valuation(value):
    value = abs(value)
    return (value & -value).bit_length() - 1 if value else 10**9


def majority_subset(length, bits):
    assert bits % 2 and 1 <= bits <= length.bit_length() - 1
    mask = (1 << bits) - 1
    return tuple(index for index in range(length)
                 if (index & mask).bit_count() > bits // 2)


def residue_data(length, subset):
    layers = []
    for level in range(1, length.bit_length() + 1):
        counts = [0] * (1 << level)
        for value in subset:
            counts[value % len(counts)] += 1
        layers.append(counts)
    return layers


def correlation(counts, shift):
    return sum(value * counts[(index + shift) % len(counts)]
               for index, value in enumerate(counts))


def content_by_residues(length, subset, layers, center):
    content = len(subset) - correlation(layers[-1], center)
    assert content >= 0
    for counts in layers[:-1]:
        difference = sum(value * value for value in counts) - correlation(counts, center)
        assert difference >= 0
        assert 2 * difference == sum((value - counts[(index + center) % len(counts)]) ** 2
                                    for index, value in enumerate(counts))
        content += difference
    return content


def content_by_cauchy(length, subset, center):
    a = length.bit_length() - 1
    vandermonde_value = sum(valuation(y - x) for index, x in enumerate(subset)
                            for y in subset[index + 1:])
    denominator_value = sum(min(valuation(center + x - y), a + 1)
                            for x in subset for y in subset)
    return len(subset) * (a + 1) + 2 * vandermonde_value - denominator_value


def polynomial_multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        if x:
            for j, y in enumerate(right):
                result[i + j] += x * y
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def polynomial_roots(roots):
    result = [1]
    for root in roots:
        result = polynomial_multiply(result, [-root, 1])
    return result


def minor_factorization(length, subset):
    size = len(subset)
    counts = Counter(y - x for x in subset for y in subset)
    roots = [root for root in range(-length, length)
             for _ in range(size - counts[root])]
    numerator = (-1) ** (size * (size - 1) // 2)
    numerator *= prod((y - x) ** 2 for index, x in enumerate(subset) for y in subset[index + 1:])
    derivative = prod(-t for t in range(-length, length) if t)
    denominator = derivative ** size
    polynomial = polynomial_roots(roots)
    assert len(polynomial) - 1 == size * (2 * length - size)
    assert numerator * polynomial[0] == denominator
    return polynomial, numerator, denominator


def local_polynomial_content(coefficients, numerator, denominator, center, scale):
    degree = len(coefficients) - 1
    local = [sum(coefficients[k] * comb(k, j) * center ** (k - j)
                 for k in range(j, degree + 1)) * scale ** j
             for j in range(degree + 1)]
    return min(valuation(value) for value in local) + valuation(numerator) - valuation(denominator)


def determinant_polynomial(length, subset):
    size = len(subset)
    entries = {}
    for difference in {y - x for x in subset for y in subset}:
        entries[difference] = polynomial_roots(root for root in range(-length, length)
                                               if root != difference)
    result = [0] * (size * (2 * length - 1) + 1)
    for permutation in permutations(range(size)):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(size) for j in range(i + 1, size))
        term = [(-1) ** inversions]
        for row, column in enumerate(permutation):
            term = polynomial_multiply(term, entries[subset[column] - subset[row]])
        for degree, coefficient in enumerate(term):
            result[degree] += coefficient
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    statistics = Counter()
    rows = []
    for a in range(1, 9):
        length = 1 << a
        for cutoff in range(1, a + 1):
            bits = cutoff if cutoff % 2 else cutoff + 1
            if bits > a:
                continue
            subset = majority_subset(length, bits)
            assert len(subset) == length // 2
            layers = residue_data(length, subset)
            middle = (bits - 1) // 2
            assert (4 * middle + 1) * comb(2 * middle, middle) ** 2 >= 16 ** middle
            for level in range(cutoff):
                counts = layers[level]
                half = 1 << level
                differences = [counts[index] - counts[index + half] for index in range(half)]
                discrepancy = sum(value * value for value in differences)
                correlation_value = sum(differences)
                expected = -comb(bits - 1, (bits - 1) // 2) * (length >> bits)
                assert correlation_value == expected
                assert half * discrepancy >= correlation_value ** 2
                assert 8 * bits * correlation_value ** 2 >= length ** 2
                assert 8 * bits * half * discrepancy >= length ** 2
                statistics['majority_discrepancy_controls'] += 1
            minimum_ratio = None
            contents = []
            for center in range(-length, length):
                if center == 0 or valuation(center) >= cutoff:
                    continue
                level = valuation(center)
                content = content_by_residues(length, subset, layers, center)
                counts = layers[level]
                discrepancy = sum((counts[index] - counts[index + (1 << level)]) ** 2
                                  for index in range(1 << level))
                assert content >= discrepancy
                assert 8 * bits * (1 << level) * content >= length ** 2
                assert 16 * cutoff * (1 << level) * content >= length ** 2
                ratio = Fraction(content * cutoff * (1 << level), length ** 2)
                minimum_ratio = ratio if minimum_ratio is None else min(minimum_ratio, ratio)
                if length <= 32:
                    assert content == content_by_cauchy(length, subset, center)
                    statistics['direct_cauchy_content_controls'] += 1
                contents.append((center, content))
                statistics['active_jet_content_controls'] += 1
            if length <= 8:
                polynomial, numerator, denominator = minor_factorization(length, subset)
                assert determinant_polynomial(length, subset) == [numerator * value for value in polynomial]
                statistics['permutation_determinant_controls'] += 1
                for center, content in contents:
                    measured = local_polynomial_content(polynomial, numerator, denominator, center, 2 * length)
                    assert measured == content
                    statistics['expanded_local_polynomial_controls'] += 1
                for copies in (2, 3):
                    if length <= 4:
                        repeated = [1]
                        for _ in range(copies):
                            repeated = polynomial_multiply(repeated, polynomial)
                        for center, content in contents:
                            assert local_polynomial_content(repeated, numerator ** copies, denominator ** copies,
                                                            center, 2 * length) == copies * content
                            statistics['repeated_packet_controls'] += 1
            rows.append({'packet_length': length, 'active_cutoff': cutoff,
                         'majority_bits': bits, 'occupied_slots': len(subset),
                         'polynomial_degree': 3 * length ** 2 // 4,
                         'minimum_normalized_content': {
                             'numerator': minimum_ratio.numerator,
                             'denominator': minimum_ratio.denominator},
                         'smallest_content': min(content for _, content in contents)})
    result = {'scope': 'Exact 0/1 rank-one state witnesses for the proposed exterior lift; no claim for the original Boolean lift or divisor predicates.',
              'counts': dict(statistics), 'cases': rows,
              'source_revision': '4eaee87548250b2682affda7d05d947c05554bae',
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    text = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
