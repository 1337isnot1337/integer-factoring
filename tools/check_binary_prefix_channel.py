#!/usr/bin/env python3
"""Finite exact support for the actual binary-prefix channel-rank theorem.

The asymptotic proof is in research/binary-prefix-top-channel.md.
This checker does not establish the rank of the complete central piece.
"""

from collections import defaultdict
import json
from math import isqrt

from check_partition_channels import operator, partitions, rank_mod


def bounded(weight: int, maximum: int) -> tuple[tuple[int, ...], ...]:
    return tuple(part for part in partitions(weight)
                 if not part or part[-1] <= maximum)


def replace_large(part: tuple[int, ...], maximum: int) -> tuple[int, ...]:
    if maximum < 1:
        raise ValueError("positive part bound required")
    small = tuple(value for value in part if value <= maximum)
    removed_weight = sum(value for value in part if value > maximum)
    return tuple(sorted((1,) * removed_weight + small))


def split_with_ones(theta: tuple[int, ...], half: int, reserved: int
                    ) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Append `reserved` ones and split into two weight-`half` partitions."""
    if not 1 <= reserved <= half or sum(theta) != 2 * half - reserved:
        raise ValueError("invalid reserved-one parameters")
    if any(value < 1 or value > reserved for value in theta):
        raise ValueError("input part exceeds reserve")
    left, total, cut = [], 0, 0
    for value in theta:
        if total + value > half:
            break
        left.append(value)
        total += value
        cut += 1
    gap = half - total
    assert 0 <= gap <= reserved
    first = tuple(sorted(tuple(left) + (1,) * gap))
    second = tuple(sorted(theta[cut:] + (1,) * (reserved - gap)))
    assert sum(first) == sum(second) == half
    return first, second


def verify() -> dict:
    for prime in (65521, 1000003):
        assert all(prime % divisor for divisor in range(2, isqrt(prime) + 1))
    fiber_cases = fiber_inputs = split_cases = 0
    for weight in range(1, 19):
        for maximum in range(1, weight + 1):
            fibers = defaultdict(list)
            for part in partitions(weight):
                output = replace_large(part, maximum)
                removed = tuple(value for value in part if value > maximum)
                assert output[-1] <= maximum and sum(output) == weight
                ones = list(output)
                for _ in range(sum(removed)):
                    ones.remove(1)
                assert tuple(sorted(tuple(ones) + removed)) == part
                fibers[output].append(removed)
                fiber_inputs += 1
            limit = weight // (maximum + 1)
            budget = sum(len(tuple(part for part in partitions(total)
                                   if len(part) <= limit))
                         for total in range(weight + 1))
            assert max(map(len, fibers.values())) <= budget
            assert len(partitions(weight)) <= len(bounded(weight, maximum)) * budget
            fiber_cases += 1

    for half in range(1, 13):
        for reserved in range(1, half + 1):
            seen = set()
            for theta in bounded(2 * half - reserved, reserved):
                first, second = split_with_ones(theta, half, reserved)
                combined = tuple(sorted(first + second))
                assert combined == tuple(sorted(theta + (1,) * reserved))
                assert combined not in seen
                assert max(first) <= reserved and max(second) <= reserved
                seen.add(combined)
                split_cases += 1
            assert len(seen) == len(bounded(2 * half - reserved, reserved))

    largest_part_cases = no_ones_inputs = omission_bound_cases = 0
    for weight in range(1, 19):
        for maximum in range(weight + 1):
            tail = len(partitions(weight)) - len(bounded(weight, maximum))
            classes = sum(len(bounded(weight - largest, largest))
                          for largest in range(maximum + 1, weight + 1))
            assert tail == classes
            assert tail <= sum(len(partitions(weight - largest))
                               for largest in range(maximum + 1, weight + 1))
            largest_part_cases += 1
        if weight < 2:
            continue
        without_ones = tuple(part for part in partitions(weight) if 1 not in part)
        assert len(without_ones) == len(partitions(weight)) - len(partitions(weight - 1))
        images = set()
        for part in without_ones:
            image = part[:-1] + (part[-1] + 1,)
            assert image[-1] > (image[-2] if len(image) > 1 else 0)
            assert image in partitions(weight + 1) and 1 not in image
            assert image[:-1] + (image[-1] - 1,) == part
            images.add(image)
            no_ones_inputs += 1
        assert len(images) == len(without_ones)
        for step in range(6):
            later_differences = [
                len(partitions(weight + offset)) - len(partitions(weight + offset - 1))
                for offset in range(step + 1)
            ]
            assert all(value >= len(without_ones) for value in later_differences)
            assert (step + 1) * len(without_ones) <= len(partitions(weight + step))
            omission_bound_cases += 1

    records = []
    prefix_cases = triangular_cases = 0
    for width in range(1, 5):
        half = 2
        while len(partitions(half)) < 1 << width:
            half += 1
        assert len(partitions(half - 1)) < 1 << width
        labels = tuple(sorted(partitions(half), key=lambda part: tuple(reversed(part))))[:1 << width]
        label_set = set(labels)
        unions = {tuple(sorted(first + second)) for first in labels for second in labels}
        for maximum in range(1, half + 1):
            permitted = bounded(half, maximum)
            if len(permitted) > len(labels):
                continue
            assert set(permitted) <= label_set
            for theta in bounded(2 * half - maximum, maximum):
                first, second = split_with_ones(theta, half, maximum)
                assert first in label_set and second in label_set
                assert tuple(sorted(first + second)) in unions
            prefix_cases += 1
        outputs = partitions(2 * half)
        matrix = []
        for first in labels:
            for second in labels:
                expected = tuple(sorted(first + second))
                row = [operator(first, output, second) for output in outputs]
                for output, value in zip(outputs, row):
                    if len(output) >= len(expected):
                        assert value == int(output == expected)
                        triangular_cases += 1
                matrix.append(row)
        ranks = [rank_mod(matrix, prime) for prime in (65521, 1000003)]
        assert all(len(unions) <= rank <= len(outputs) for rank in ranks)
        records.append({"label_bits": width, "weight": half,
                        "binary_labels": len(labels), "distinct_unions": len(unions),
                        "top_channel_dimension": len(outputs),
                        "rank_mod_65521": ranks[0], "rank_mod_1000003": ranks[1]})
        operator.cache_clear()

    # A suffix of the same size need not contain the bounded initial segment.
    ordered = tuple(sorted(partitions(4), key=lambda part: tuple(reversed(part))))
    assert set(bounded(4, 3)) <= set(ordered[:4])
    assert not set(bounded(4, 3)) <= set(ordered[-4:])

    # The reserve is required: (2,2,2) alone cannot split into two weight-three parts.
    reachable = {0}
    for value in (2, 2, 2):
        reachable |= {total + value for total in tuple(reachable)}
    assert 3 not in reachable

    return {
        "scope": "finite combinatorial maps and actual binary-prefix channel ranks; no complete-piece rank claim",
        "replacement_fiber_parameter_cases": fiber_cases,
        "replacement_map_inputs": fiber_inputs,
        "reserved_one_splits": split_cases,
        "largest_part_class_cases": largest_part_cases,
        "no_ones_injection_inputs": no_ones_inputs,
        "prefix_omission_bound_cases": omission_bound_cases,
        "bounded_prefix_cases": prefix_cases,
        "triangular_coefficient_cases": triangular_cases,
        "binary_channel_cases": records,
        "negative_controls": ["lexicographic suffix instead of prefix", "omitted reserved ones"],
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
