#!/usr/bin/env python3
"""Exact checks for the regular fusion rank and subset-sum interval proof."""

import argparse
from itertools import product
import json
from pathlib import Path
import random

from check_modular_crossing import (
    column_bits, determinant_operator_two, homogeneous_family_two,
    homogeneous_piece_two, inverse_two, modular_operator, rank_two,
    regular_partitions,
)


def subset_sums(parts):
    bits = 1
    for part in parts:
        bits |= bits << part
    return bits


def interval_checks():
    cases = extensions = 0
    rng = random.Random(301026)
    for count in range(1, 13):
        for odd in range(1, 2 * count, 2):
            pairs = [(101 + 4 * i, 103 + 4 * i) for i in range(count)]
            parts = [value for pair in pairs for value in pair] + [odd]
            lower = sum(pair[0] for pair in pairs) + odd
            upper = sum(pair[0] for pair in pairs) + 2 * count
            bits = subset_sums(parts)
            assert lower + upper == sum(parts)
            assert ((bits >> lower) & ((1 << (upper - lower + 1)) - 1)) == (
                (1 << (upper - lower + 1)) - 1)
            cases += 1
            for _ in range(30):
                choices = [j for j in range(1, min(81, upper - lower + 1) + 1, 2)
                           if j not in parts]
                if not choices:
                    break
                value = rng.choice(choices)
                parts.append(value)
                bits |= bits << value
                upper += value
                assert lower + upper == sum(parts)
                assert ((bits >> lower) & ((1 << (upper - lower + 1)) - 1)) == (
                    (1 << (upper - lower + 1)) - 1)
                extensions += 1
    return {"initial_pair_intervals": cases, "successive_interval_extensions": extensions}


def partition_profiles():
    maximum = 128
    counts = [1] + [0] * maximum
    for part in range(1, maximum + 1, 2):
        for total in range(maximum, part - 1, -1):
            counts[total] += counts[total - part]
    profiles = []
    for weight in (4, 8, 16, 24, 32, 48, 64):
        labels = regular_partitions(2, weight)
        outputs = regular_partitions(2, 2 * weight)
        assert len(labels) == counts[weight] and len(outputs) == counts[2 * weight]
        splittable = sum((subset_sums(parts) >> weight) & 1 for parts in outputs)
        profiles.append({"weight": weight, "labels": len(labels),
                         "output_dimension": len(outputs), "splittable_outputs": splittable})
    return profiles


def matrix_square(rows):
    result = []
    for row in rows:
        value = 0
        for index in range(len(rows)):
            if row >> index & 1:
                value ^= rows[index]
        result.append(value)
    return result


def fusion_profile(weight):
    labels, duals = homogeneous_family_two(weight)
    output = regular_partitions(2, 2 * weight)
    dimension = len(labels)
    columns = []
    filtration_entries = independent_entries = 0
    rng = random.Random(73510 + weight)
    for first, second in product(labels, repeat=2):
        vector = 0
        top_length = len(first) + len(second)
        top = tuple(sorted(first + second)) if set(first).isdisjoint(second) else None
        for index, target in enumerate(output):
            value = determinant_operator_two(first, target, second)
            if len(target) >= top_length:
                assert value == int(target == top)
                filtration_entries += 1
            if weight <= 16 or rng.randrange(257) == 0:
                assert value == modular_operator(2, first, target, second)
                independent_entries += 1
            vector |= value << index
        columns.append(vector)
    fusion_rank = rank_two(columns)
    balanced = sum((subset_sums(parts) >> weight) & 1 for parts in output)
    assert fusion_rank >= balanced

    gram = [sum(determinant_operator_two(first, (), second) << j
                for j, second in enumerate(labels)) for first in labels]
    assert matrix_square(gram) == [1 << i for i in range(dimension)]
    assert inverse_two(gram, dimension) == gram
    transformed = []
    for first in range(dimension):
        for second in range(dimension):
            vector = 0
            for term in range(dimension):
                if gram[first] >> term & 1:
                    vector ^= columns[term * dimension + second]
            transformed.append(vector)
    predicted = []
    for i, vector in enumerate(transformed):
        first, second = divmod(i, dimension)
        row = 0
        for j, other in enumerate(transformed):
            third, fourth = divmod(j, dimension)
            value = ((vector & other).bit_count() & 1) ^ (
                first == second and third == fourth)
            row |= int(value) << j
        predicted.append(row)
    actual = homogeneous_piece_two(weight, 0, labels, duals)
    assert actual == predicted
    row_rank = rank_two(actual)
    column_rank = rank_two(column_bits(actual, dimension))
    assert row_rank == column_rank
    assert max(0, 2 * fusion_rank - len(output) - 1) <= row_rank <= fusion_rank + 1
    return {"weight": weight, "labels": dimension, "output_dimension": len(output),
            "splittable_outputs": balanced, "fusion_rank": fusion_rank,
            "central_row_rank": row_rank, "central_column_rank": column_rank,
            "top_length_entries": filtration_entries,
            "independent_source_recursion_entries": independent_entries,
            "central_gram_entries": dimension ** 4}


def ideal_checks():
    states = tuple(part for weight in range(13) for part in regular_partitions(2, weight))
    parameters = [(m,) for m in range(2, 11, 2)] + [(m, m) for m in (1, 3, 5)]
    count = 0
    for parameter, first, second in product(parameters, states, states):
        assert modular_operator(2, parameter, first, second) == 0
        count += 1
    return count


def verify(large=False):
    intervals = interval_checks()
    profiles = [fusion_profile(h) for h in ((4, 8, 16, 32) if large else (4, 8, 16))]
    return {"scope": "exact local identities supporting an ordinary asymptotic proof",
            "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
            **intervals, "parameter_ideal_entries": ideal_checks(),
            "partition_profiles": partition_profiles(), "fusion_profiles": profiles,
            "asymptotic_statement": "full regular homogeneous central ranks have rho=1/sqrt(2)",
            "limitations": ["full regular label domain", "no arbitrary-prefix claim",
                            "no improved factoring exponent"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight-thirty-two", action="store_true")
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    result = verify(options.weight_thirty_two)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if options.output:
        options.output.write_text(encoded)
    print(encoded, end="")
