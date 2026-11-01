#!/usr/bin/env python3
"""Check the Wick engine and the four-piece, three-label certificate."""

from math import comb
import argparse
import json

from check_partition_channels import dual_polynomial
from check_modular_crossing import (
    column_bits, modular_operator, normalized_piece, rank_two,
    regular, regular_partitions,
)
from fermionic_piece_polynomial import piece_matrices, reflected_duals


def direct_generalized_piece(labels, dual_monomials, shift, exchanged=True):
    """Independent literal source recursion with parameter h+weight(a)."""
    h = sum(labels[0])
    a_weight = sum(dual_monomials[0])
    assert all(sum(a) == a_weight for a in dual_monomials)
    locality = a_weight + h
    dimension = len(labels)
    result = [[0] * (dimension * dimension) for _ in range(dimension * dimension)]
    for index in range(locality + 1):
        if comb(locality, index) % 2 == 0:
            continue
        energy = shift + index if exchanged else 2 * h - shift - index
        for intermediate in regular_partitions(2, energy):
            if exchanged:
                left = [modular_operator(2, a, intermediate, y)
                        for a in dual_monomials for y in labels]
                right = [modular_operator(2, v, u, intermediate)
                         for v in labels for u in labels]
                nonzero = [j for j, value in enumerate(right) if value]
                for row, value in zip(result, left):
                    if value:
                        for j in nonzero:
                            row[j] ^= 1
            else:
                for x, a in enumerate(dual_monomials):
                    for u_index, u in enumerate(labels):
                        if modular_operator(2, a, u, intermediate):
                            for v_index, v in enumerate(labels):
                                for y_index, y in enumerate(labels):
                                    result[x * dimension + y_index][v_index * dimension + u_index] ^= (
                                        modular_operator(2, v, intermediate, y))
    return [sum(value << j for j, value in enumerate(row)) for row in result]


def verify():
    source_entries = 0
    for h in (1, 4, 8):
        labels = regular_partitions(2, h)
        duals = [tuple((a, c) for a, c in dual_polynomial(label)
                       if regular(a, 2) and c % 2) for label in labels]
        matrices = piece_matrices(labels, duals)
        dimension = len(labels)
        for shift in range(-2 * h, 2 * h + 1):
            source = normalized_piece(2, h, shift, labels)
            rows = [sum(c << j for j, c in enumerate(row)) for row in source]
            assert matrices.get(shift, [0] * (dimension * dimension)) == rows
            source_entries += dimension ** 4

    labels = ((1, 11), (3, 9), (5, 7))
    duals = reflected_duals(labels, 8)
    monomials = tuple(terms[0][0] for terms in duals)
    assert monomials == ((5, 15), (7, 13), (9, 11))
    dual_entries = 0
    for x, a in enumerate(monomials):
        for v, b in enumerate(labels):
            assert modular_operator(2, a, (), b) == int(x == v)
            dual_entries += 1
    matrices = piece_matrices(labels, duals)
    expected = {
        -32: [273, 0, 0, 0, 273, 0, 0, 0, 273],
        -16: [256, 0, 0, 0, 256, 0, 0, 0, 17],
        0: [272, 2, 0, 8, 257, 0, 0, 0, 17],
        16: [256, 0, 4, 0, 256, 32, 64, 128, 17],
    }
    assert matrices == expected
    generalized_entries = 0
    for shift in range(-32, 25):
        row_order = direct_generalized_piece(labels, monomials, shift)
        column_order = direct_generalized_piece(labels, monomials, shift, exchanged=False)
        assert row_order == column_order == matrices.get(shift, [0] * 9)
        generalized_entries += 81
    total = [0] * 9
    profile = []
    for shift, rows in sorted(matrices.items()):
        total = [a ^ b for a, b in zip(total, rows)]
        profile.append([shift, rank_two(rows), rank_two(column_bits(rows, 3))])
    assert total == [1 << i for i in range(9)]
    assert profile == [[-32, 1, 9], [-16, 2, 4], [0, 4, 5], [16, 6, 2]]
    assert max(row * column for _, row, column in profile) == 20
    # The certificate is a literal sum; dropping its highest-shift piece fails.
    assert [a ^ b for a, b in zip(total, matrices[16])] != total

    return {
        "status": "pass",
        "scope": "exact finite tensor improvement; no asymptotic contraction claim",
        "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5",
        "source_entry_comparisons": source_entries,
        "reflected_duality_entries": dual_entries,
        "generalized_both_order_entry_comparisons": generalized_entries,
        "summed_crossing_entries": 81,
        "certificate": {
            "characteristic": 2,
            "labels": labels,
            "dual_monomials": monomials,
            "input_weight": 12,
            "dual_weight": 20,
            "locality_parameter": 32,
            "row_bitsets": {str(shift): rows for shift, rows in sorted(matrices.items())},
            "row_column_rank_profile": profile,
            "maximum_rank_product": 20,
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output")
    args = parser.parse_args()
    result = json.dumps(verify(), indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(result)
    print(result, end="")
