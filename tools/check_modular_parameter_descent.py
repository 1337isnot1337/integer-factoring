#!/usr/bin/env python3
"""Exact finite controls for descent in the source operator parameter.

Source recursion: anthropics/formal-math, Appendix B at 16a5ce4.
The proof and scope are in research/modular-parameter-descent.md.
No asymptotic rank or factoring improvement is asserted.
"""

from collections import Counter
from functools import lru_cache
from math import comb
import json

from check_modular_crossing import (
    modular_angle, modular_operator, normalized_piece, regular,
    regular_partitions,
)
from check_partition_channels import dual_polynomial, operator, partitions


@lru_cache(None)
def reduced_dual(prime, parts):
    """Compile directly in the quotient; never build forbidden monomials."""
    result = {(): 1}
    for part in parts:
        next_result = {}
        sign = (-1) ** (part - 1)
        for monomial, value in result.items():
            counts = Counter(monomial)
            for index in range(1, part + 1):
                coefficient = sign * comb(part, index) % prime
                if not coefficient or index % prime == 0:
                    continue
                if counts[index] == prime - 1:
                    continue
                target = tuple(sorted(monomial + (index,)))
                next_result[target] = (
                    next_result.get(target, 0) + value * coefficient
                ) % prime
        result = {m: c for m, c in next_result.items() if c}
    return tuple(sorted(result.items()))


def reduced_angle(prime, parts, index, output, input_):
    weight = sum(parts)
    return sum(coefficient * comb(weight + sum(term), index)
               * modular_operator(prime, term, output, input_)
               for term, coefficient in reduced_dual(prime, parts)
               if weight + sum(term) >= index) % prime


def reduced_piece(prime, weight, shift, labels):
    """Use the existing factorization, replacing only its dual compiler."""
    from check_partition_channels import partition_denominator
    dimension = len(labels)
    size = dimension ** 2
    result = [[0] * size for _ in range(size)]
    scales = [pow(partition_denominator(v) % prime, -1, prime)
              for v in labels]
    for index in range(2 * weight + 1):
        for intermediate in regular_partitions(prime, shift + index):
            left = [(-1) ** index * reduced_angle(
                prime, x, index, intermediate, y) % prime
                for x in labels for y in labels]
            right = [scale * modular_operator(prime, v, u, intermediate) % prime
                     for v, scale in zip(labels, scales) for u in labels]
            for row, value in zip(result, left):
                if value:
                    for column, other in enumerate(right):
                        row[column] = (row[column] + value * other) % prime
    return result


def verify():
    parameter_checks = 0
    parameter_records = []
    for prime in (2, 3, 5, 7):
        states = tuple(x for h in range(7)
                       for x in regular_partitions(prime, h))
        bad = tuple(x for h in range(1, 9) for x in partitions(h)
                    if not regular(x, prime))
        for parts in bad:
            for output in states:
                for input_ in states:
                    assert modular_operator(prime, parts, output, input_) == 0
                    parameter_checks += 1
        parameter_records.append({
            "prime": prime, "max_parameter_weight": 8,
            "max_state_weight": 6, "nonregular_parameters": len(bad),
            "regular_states": len(states),
        })

    angle_checks = 0
    dual_records = []
    for prime in (2, 3, 5):
        states = tuple(x for h in range(5)
                       for x in regular_partitions(prime, h))
        labels = tuple(x for h in range(1, 10)
                       for x in regular_partitions(prime, h))
        full_terms = reduced_terms = 0
        for parts in labels:
            expected = tuple(sorted((m, c % prime)
                             for m, c in dual_polynomial(parts)
                             if regular(m, prime) and c % prime))
            assert reduced_dual(prime, parts) == expected
            full_terms += len(dual_polynomial(parts))
            reduced_terms += len(expected)
            for output in states:
                for input_ in states:
                    for index in range(2 * sum(parts) + 1):
                        assert reduced_angle(prime, parts, index, output, input_) == (
                            modular_angle(prime, parts, index, output, input_))
                        angle_checks += 1
        dual_records.append({
            "prime": prime, "max_label_weight": 9,
            "max_state_weight": 4, "labels": len(labels),
            "full_integer_dual_terms": full_terms,
            "nonzero_reduced_dual_terms": reduced_terms,
        })

    crossing_entries = piece_entries = 0
    crossing_records = []
    for prime in (2, 3, 5):
        weight = next(h for h in range(1, 12)
                      if len(regular_partitions(prime, h)) >= 2)
        labels = regular_partitions(prime, weight)[:2]
        size = 4
        total = [[0] * size for _ in range(size)]
        for shift in range(-2 * weight, 2 * weight + 1):
            piece = reduced_piece(prime, weight, shift, labels)
            assert piece == normalized_piece(prime, weight, shift, labels)
            piece_entries += size ** 2
            for row in range(size):
                for col in range(size):
                    total[row][col] = (total[row][col] + piece[row][col]) % prime
        assert total == [[int(i == j) for j in range(size)]
                         for i in range(size)]
        crossing_entries += size ** 2
        crossing_records.append({"prime": prime, "weight": weight,
                                 "labels": labels})

    # The same deletion is false before reducing the coefficients modulo p.
    # t_2 annihilates t_1 with coefficient -2 over the integers.
    assert operator((2,), (), (1,)) == -2
    assert modular_operator(2, (2,), (), (1,)) == 0
    # A surviving parameter cannot be erased: t_1[-1]t_1 = 1.
    assert modular_operator(2, (1,), (), (1,)) == 1
    return {
        "status": "exact finite controls passed; ordinary proof in note",
        "parameter_ideal_entries": parameter_checks,
        "parameter_ranges": parameter_records,
        "dual_angle_entries": angle_checks,
        "dual_ranges_and_term_counts": dual_records,
        "unchanged_piece_entries": piece_entries,
        "complete_crossing_entries": crossing_entries,
        "crossing_examples": crossing_records,
        "negative_controls": 2,
        "factoring_score_improvement": False,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
