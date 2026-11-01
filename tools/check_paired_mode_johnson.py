#!/usr/bin/env python3
"""Cross-check the paired-mode reduction and report exact finite rank profiles."""

from itertools import product
from math import log2
import argparse
import json

from check_fermionic_piece_polynomial import direct_generalized_piece
from check_modular_crossing import column_bits, rank_two
from fermionic_piece_polynomial import monomial_piece_polynomial
from paired_mode_johnson import complete_entry, expand_power, paired_labels, rank_profile, subsets


def verify(max_pairs=10):
    wick_entries = source_entries = 0
    for pair_count, q in ((2, 1), (3, 1), (3, 2), (4, 1), (4, 2)):
        box = 1 << (2 * pair_count - 1).bit_length()
        labels = paired_labels(pair_count, q, box)
        codes = subsets(pair_count, q)
        dimension = len(labels)
        matrices = {}
        for x, y, v, u in product(range(dimension), repeat=4):
            numerator, offset = complete_entry(q, codes[x], codes[y], codes[u], codes[v])
            expected = expand_power(numerator, 4 * box, offset)
            actual, actual_offset = monomial_piece_polynomial(labels[x], labels[v], labels[u], labels[y])
            actual = expand_power(actual, 1, actual_offset)
            assert actual == expected
            wick_entries += 1
            for shift in expected:
                matrices.setdefault(shift, [0] * dimension ** 2)[x * dimension + y] |= (
                    1 << (v * dimension + u))
        total = [0] * dimension ** 2
        actual_profile = []
        for shift, rows in sorted(matrices.items()):
            total = [a ^ b for a, b in zip(total, rows)]
            actual_profile.append([shift // (4 * box), rank_two(rows),
                                   rank_two(column_bits(rows, dimension))])
        assert total == [1 << i for i in range(dimension ** 2)]
        expected_profile = [row for row in rank_profile(pair_count, q)["profile"] if row[1]]
        assert actual_profile == expected_profile
        if q == 1 and pair_count <= 3:
            weight = sum(labels[0])
            for shift in range(-2 * weight, 2 * weight + 1):
                expected = matrices.get(shift, [0] * dimension ** 2)
                assert direct_generalized_piece(labels, labels, shift) == expected
                assert direct_generalized_piece(labels, labels, shift, exchanged=False) == expected
                source_entries += dimension ** 4
    reports = []
    for pairs in range(2, max_pairs + 1):
        report = rank_profile(pairs, pairs // 2)
        report["finite_log_rank_ratio"] = (
            log2(report["maximum_rank_product"]) / (4 * log2(report["label_dimension"])))
        reports.append(report)
    # Replacing the reciprocal profile by an identical profile changes an axis.
    assert reports[0]["profile"][0][1:] == [1, 4]
    return {"status": "pass", "wick_entry_comparisons": wick_entries,
            "literal_source_both_order_entries": source_entries,
            "profiles": reports,
            "scope": "exact all-size block identity and finite ranks; no exponent improvement established",
            "source_revision": "16a5ce41304fcdaea4e0be7c132e1c361dc26da5"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-pairs", type=int, default=10)
    parser.add_argument("--output")
    args = parser.parse_args()
    encoded = json.dumps(verify(args.max_pairs), indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path
        Path(args.output).write_text(encoded)
    print(encoded, end="")
