#!/usr/bin/env python3
"""Exact checks for partition-defect labels and their existing top channel.

The proof and limitations are in research/near-ones-labels.md.
Small full ranks do not prove an asymptotic rank bound for multiplicity boxes.
The operator is the pinned Appendix B construction; it is not a new algorithm.
"""

import json
from math import factorial, isqrt

from check_partition_channels import operator, partitions, rank_mod


def embed(defect_partition: tuple[int, ...], weight: int) -> tuple[int, ...]:
    defect = sum(defect_partition)
    ones = weight - defect - len(defect_partition)
    if ones < 0:
        raise ValueError("weight too small for this defect partition")
    return (1,) * ones + tuple(value + 1 for value in defect_partition)


def verify() -> dict:
    for prime in (65521, 1000003):
        assert all(prime % divisor for divisor in range(2, isqrt(prime) + 1))

    embedding_inputs = union_cases = 0
    for defect in range(1, 11):
        for weight in (2 * defect, 3 * defect + 1):
            images = [embed(part, weight) for part in partitions(defect)]
            assert len(set(images)) == len(partitions(defect))
            for part, image in zip(partitions(defect), images):
                assert sum(image) == weight
                assert weight - len(image) == defect
                assert tuple(value - 1 for value in image if value > 1) == part
                embedding_inputs += 1
            for first in partitions(defect):
                for second in partitions(defect):
                    original_union = tuple(sorted(first + second))
                    image_union = tuple(sorted(embed(first, weight) + embed(second, weight)))
                    assert image_union == embed(original_union, 2 * weight)
                    union_cases += 1

    escaping_coefficients = []
    for weight in range(2, 13):
        value = operator((1,) * weight, (1, 2 * weight - 1), (1,) * weight)
        assert value == weight * factorial(weight)
        escaping_coefficients.append({
            "weight": weight,
            "input_defects": [0, 0],
            "output_defect": 2 * weight - 2,
            "coefficient": value,
        })
    operator.cache_clear()

    records = []
    support_cases = triangular_cases = 0
    for weight in range(2, 9):
        labels = [(1,) * (weight - 2 * count) + (2,) * count
                  for count in range(weight // 2 + 1)]
        outputs = partitions(2 * weight)
        unions = {tuple(sorted(first + second)) for first in labels for second in labels}
        assert len(unions) == 2 * len(labels) - 1
        matrix = []
        for first in labels:
            for second in labels:
                first_defect = weight - len(first)
                second_defect = weight - len(second)
                leading = tuple(sorted(first + second))
                row = []
                for output in outputs:
                    value = operator(first, output, second)
                    row.append(value)
                    if value:
                        output_defect = 2 * weight - len(output)
                        extra = output_defect - first_defect - second_defect
                        assert extra >= 0 and extra % 2 == 0
                        annihilations = extra // 2
                        assert annihilations < len(first)
                        assert annihilations <= len(second)
                        support_cases += 1
                    if len(output) >= len(leading):
                        assert value == int(output == leading)
                        triangular_cases += 1
                matrix.append(row)
        ranks = [rank_mod(matrix, prime) for prime in (65521, 1000003)]
        assert ranks == [len(labels) ** 2] * 2
        records.append({
            "weight": weight,
            "labels": len(labels),
            "distinct_unions": len(unions),
            "pair_dimension": len(labels) ** 2,
            "rank_mod_65521": ranks[0],
            "rank_mod_1000003": ranks[1],
        })
        operator.cache_clear()

    return {
        "source_commit": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "scope": "defect embeddings, exact support identities, and finite top-channel ranks; no new factoring rate",
        "embedding_inputs": embedding_inputs,
        "union_embedding_cases": union_cases,
        "escaping_all_ones_coefficients": escaping_coefficients,
        "nonzero_defect_support_cases": support_cases,
        "triangular_coefficient_cases": triangular_cases,
        "one_two_label_cases": records,
        "negative_control": "cutting output defect at the sum of the two input defects removes a nonzero all-ones coefficient",
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
